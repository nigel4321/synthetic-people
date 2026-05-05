#!/usr/bin/env python3
"""QC + validation for a single VCF before downstream processing.

Two tiers of check:

Hard (errors) — these mean the file is unreadable or unusable downstream:
  - file exists and is non-empty
  - tabix index is present alongside the VCF
  - bcftools can parse the header
  - header has a #CHROM line with sample columns
  - index reports at least one contig
  - file contains at least one variant record

Soft (warnings) — the file is readable but its attributes do not match
what a human-genome VCF typically looks like. These never abort the
pipeline; they surface in qc_report.md for the user to eyeball:
  - reference field names a recognised human build
    (hg19/GRCh37/b37/hs37d5 or hg38/GRCh38/b38)
  - indexed contigs are among the expected human chromosomes
  - FORMAT=GT is declared
  - INFO declares AF, or both AC and AN (required by scan_variant)
  - variant positions fit inside plausible human chromosome lengths

Always writes the JSON output file. In --strict mode, also exits 1 when
any hard check fails, which aborts the Nextflow workflow.

Optional --mqc-out writes a MultiQC custom-content sidecar (a
``*_mqc.json`` file) containing the same per-file findings reshaped as a
table section — one section per VCF, merged into the cohort MultiQC
report when the pipeline runs the MULTIQC stage.
"""

import argparse
import json
import os
import re
import subprocess
import sys


# Canonical (non-prefixed) human chromosome names, plus chr-prefixed forms.
_HUMAN_CANONICAL = set(str(i) for i in range(1, 23)) | {"X", "Y", "M", "MT"}
HUMAN_CHROMS = _HUMAN_CANONICAL | {f"chr{c}" for c in _HUMAN_CANONICAL}

# Lowercased substrings that indicate a recognised human reference build.
HUMAN_REF_TOKENS = {
    "GRCh37": ("grch37", "b37", "hs37d5", "hg19"),
    "GRCh38": ("grch38", "b38", "hg38", "hs38"),
}

# Upper bound for plausible positions on any human chromosome (GRCh37 chr1
# is ~249 Mb — anything past ~300 Mb indicates wrong genome / scaffold).
MAX_HUMAN_CHROM_LENGTH = 300_000_000


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


def _detect_build(reference: str) -> str | None:
    ref = reference.lower()
    for build, tokens in HUMAN_REF_TOKENS.items():
        if any(tok in ref for tok in tokens):
            return build
    return None


def _parse_header(header_text: str) -> dict:
    chrom_line = None
    info_ids: set[str] = set()
    reference = ""
    has_gt = False
    for ln in header_text.splitlines():
        if ln.startswith("#CHROM"):
            chrom_line = ln
        elif ln.startswith("##reference"):
            reference = ln.split("=", 1)[1].strip() if "=" in ln else ""
        elif ln.startswith("##INFO=<ID="):
            m = re.match(r"##INFO=<ID=([A-Za-z0-9_]+)", ln)
            if m:
                info_ids.add(m.group(1))
        elif ln.startswith("##FORMAT=<ID=GT"):
            has_gt = True
    samples: list[str] = []
    if chrom_line:
        fields = chrom_line.split("\t")
        samples = fields[9:] if len(fields) > 9 else []
    return {
        "chrom_line_present": chrom_line is not None,
        "samples": samples,
        "reference": reference,
        "info_ids": sorted(info_ids),
        "has_format_gt": has_gt,
    }


def _parse_index(idx_stdout: str) -> tuple[list[str], int | None]:
    """Extract (contigs, total_records) from `bcftools index -s` output."""
    contigs = []
    total = 0
    have_counts = False
    for line in idx_stdout.strip().splitlines():
        parts = line.split()
        if not parts:
            continue
        contigs.append(parts[0])
        if len(parts) >= 3:
            try:
                total += int(parts[2])
                have_counts = True
            except ValueError:
                pass
    return contigs, (total if have_counts else None)


def _has_any_record(vcf: str, contig: str) -> bool:
    """True if the VCF has at least one variant record on `contig`.

    Implemented via Popen + readline() so the check is O(1) even on very
    large contigs (chr20 in 1000G Phase 3 is 1.7M records × 2504 samples).
    """
    proc = subprocess.Popen(
        ["bcftools", "view", "-H", "-r", contig, vcf],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    assert proc.stdout is not None
    first = proc.stdout.readline()
    proc.stdout.close()
    proc.wait()
    return bool(first.strip())


def _contigs_with_overlong_positions(vcf: str, contigs: list[str],
                                     threshold: int) -> list[str]:
    """Return contigs that have any variant past `threshold` (via tabix seek).

    Cheap: a single region query `contig:threshold+1-999999999`. tabix returns
    zero rows if no such variant exists, so this is O(1) per contig even on
    a 300 MB chromosome, vs. O(n_variants) for a full-scan max-position.
    """
    suspect: list[str] = []
    for c in contigs:
        region = f"{c}:{threshold + 1}-999999999"
        # -H: headers off. Capture only the first line of output via Popen so
        # we don't buffer the whole matching region into memory.
        proc = subprocess.Popen(
            ["bcftools", "view", "-H", "-r", region, vcf],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
        )
        assert proc.stdout is not None
        first = proc.stdout.readline()
        proc.stdout.close()
        proc.wait()
        if first.strip():
            suspect.append(c)
    return suspect


def run_qc(vcf_path: str, name: str) -> dict:
    errors: list[str] = []
    warnings: list[str] = []
    checks: dict = {}

    # --- Hard: file present ---
    if not os.path.isfile(vcf_path):
        errors.append(f"file does not exist: {vcf_path}")
        return _result(name, vcf_path, errors, warnings, checks)
    if os.path.getsize(vcf_path) == 0:
        errors.append(f"file is empty: {vcf_path}")
        return _result(name, vcf_path, errors, warnings, checks)
    checks["size_bytes"] = os.path.getsize(vcf_path)

    # --- Hard: tabix index ---
    tbi = vcf_path + ".tbi"
    csi = vcf_path + ".csi"
    if not (os.path.isfile(tbi) or os.path.isfile(csi)):
        errors.append(f"missing tabix index: expected {tbi} (or .csi)")
    checks["tabix_index_present"] = os.path.isfile(tbi) or os.path.isfile(csi)

    # --- Hard: header parseable ---
    header_proc = _run(["bcftools", "view", "-h", vcf_path])
    if header_proc.returncode != 0:
        errors.append(
            f"bcftools cannot read header: {header_proc.stderr.strip()[:200]}"
        )
        return _result(name, vcf_path, errors, warnings, checks)
    parsed = _parse_header(header_proc.stdout)
    checks["header_parseable"] = True
    checks["reference"] = parsed["reference"]
    checks["info_fields"] = parsed["info_ids"]
    checks["has_format_gt"] = parsed["has_format_gt"]
    checks["sample_count"] = len(parsed["samples"])

    # --- Hard: #CHROM line + samples ---
    if not parsed["chrom_line_present"]:
        errors.append("header has no #CHROM line")
    elif not parsed["samples"]:
        errors.append("#CHROM line has no sample columns")

    # --- Hard: indexed contigs + record count ---
    idx_proc = _run(["bcftools", "index", "-s", vcf_path])
    contigs, total_records = _parse_index(idx_proc.stdout)
    # Older tabix indices (e.g. 1000G Phase 3) lack the count metadata that
    # `bcftools index -s` needs — it then exits 1 with empty stdout. Fall back
    # to `tabix -l`, which works for any tabix-format index.
    if not contigs:
        tbx = _run(["tabix", "-l", vcf_path])
        contigs = [c.strip() for c in tbx.stdout.splitlines() if c.strip()]
    checks["contigs"] = contigs
    checks["n_variants"] = total_records
    if not contigs:
        errors.append("no contigs found in index "
                      "(`bcftools index -s` and `tabix -l` both empty)")
    if total_records is not None and total_records == 0:
        errors.append("file contains zero variant records")
    elif total_records is None and contigs:
        # No per-contig counts available — peek at the first record via Popen
        # so we don't buffer the whole contig into memory (chr20 in the 1000G
        # release is 1.7M records × 2.5k samples — several GB).
        if not _has_any_record(vcf_path, contigs[0]):
            errors.append("file contains no variant records in first contig")

    # --- Soft: reference build ---
    build = _detect_build(parsed["reference"])
    checks["reference_build"] = build
    if not build:
        warnings.append(
            f"reference '{parsed['reference'] or '(none)'}' does not match "
            "any recognised human build (hg19/GRCh37/hg38/GRCh38)"
        )

    # --- Soft: human contigs ---
    human_contigs = [c for c in contigs if c in HUMAN_CHROMS]
    non_human = [c for c in contigs if c not in HUMAN_CHROMS]
    checks["human_contigs"] = human_contigs
    checks["non_human_contigs"] = non_human
    if contigs and not human_contigs:
        warnings.append(
            f"no indexed contigs match human chromosome names; saw {contigs}"
        )
    elif non_human:
        warnings.append(
            f"indexed contigs include names that are not standard human "
            f"chromosomes: {non_human}"
        )

    # --- Soft: FORMAT=GT ---
    if not parsed["has_format_gt"]:
        warnings.append("header does not declare FORMAT=GT")

    # --- Soft: INFO fields the downstream pipeline relies on ---
    has_af = "AF" in parsed["info_ids"]
    has_ac_an = "AC" in parsed["info_ids"] and "AN" in parsed["info_ids"]
    checks["info_has_af"] = has_af
    checks["info_has_ac_an"] = has_ac_an
    if not has_af and not has_ac_an:
        warnings.append(
            "header declares neither AF nor AC+AN — scan_variant will not be "
            "able to classify allele frequency for this file"
        )

    # --- Soft: variant position sanity ---
    if contigs and not errors:
        # Only when the rest of the file checks out — otherwise don't spend
        # the I/O.
        overlong = _contigs_with_overlong_positions(
            vcf_path, contigs, MAX_HUMAN_CHROM_LENGTH,
        )
        checks["contigs_with_overlong_positions"] = overlong
        if overlong:
            warnings.append(
                f"contigs have variants past plausible human chromosome "
                f"length ({MAX_HUMAN_CHROM_LENGTH:,} bp): {overlong}"
            )

    return _result(name, vcf_path, errors, warnings, checks)


def _result(name: str, vcf_path: str, errors: list[str],
            warnings: list[str], checks: dict) -> dict:
    return {
        "name": name,
        "file": os.path.basename(vcf_path),
        "pass": not errors,
        "errors": errors,
        "warnings": warnings,
        "checks": checks,
    }


def build_mqc_payload(result: dict) -> dict:
    """Reshape a qc_validate result into a MultiQC custom-content section.

    MultiQC merges files that share the top-level ``id`` by concatenating
    their ``data`` dicts, so each invocation emits one section with one
    sample and the cohort report stitches itself together at MultiQC time.

    The headers list defines the column order, scaling, and tooltips; we
    keep it identical across files so the merged table stays consistent.
    """
    c = result.get("checks", {})
    if c.get("info_has_af"):
        af_acan = "AF"
    elif c.get("info_has_ac_an"):
        af_acan = "AC+AN"
    else:
        af_acan = "-"
    sample = result["name"]
    row = {
        "qc_status": "PASS" if result.get("pass") else "FAIL",
        "n_errors": len(result.get("errors", [])),
        "n_warnings": len(result.get("warnings", [])),
        "reference_build": c.get("reference_build") or "-",
        "n_variants": c.get("n_variants") if c.get("n_variants") is not None else 0,
        "sample_count": c.get("sample_count", 0),
        "has_gt": "yes" if c.get("has_format_gt") else "no",
        "has_af_or_acan": af_acan,
        "non_human_contigs": len(c.get("non_human_contigs", []) or []),
    }
    return {
        "id": "vcf_qc",
        "section_name": "VCF QC checks",
        "section_anchor": "vcf-qc",
        "description": (
            "Custom QC from qc_validate.py — header completeness, contig "
            "sanity, INFO/FORMAT presence, position sanity."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "vcf_qc_table",
            "title": "VCF QC",
            "save_file": True,
            "no_violin": True,
        },
        "headers": {
            "qc_status":         {"title": "QC",              "description": "PASS or FAIL — hard checks",                "scale": False},
            "n_errors":          {"title": "Errors",          "description": "Hard-check error count",                    "scale": "Reds",    "format": "{:,d}"},
            "n_warnings":        {"title": "Warns",           "description": "Soft-check warning count",                  "scale": "Oranges", "format": "{:,d}"},
            "reference_build":   {"title": "Build",           "description": "Detected reference build (GRCh37/GRCh38)",  "scale": False},
            "n_variants":        {"title": "# variants",      "description": "Total variant records (from index)",        "format": "{:,d}"},
            "sample_count":      {"title": "# samples",       "description": "Sample columns in #CHROM",                  "format": "{:,d}"},
            "has_gt":            {"title": "GT",              "description": "FORMAT=GT declared",                        "scale": False},
            "has_af_or_acan":    {"title": "AF/AC+AN",        "description": "Frequency-classifiable INFO available",     "scale": False},
            "non_human_contigs": {"title": "Non-human ctgs",  "description": "Indexed contigs not on standard human list", "scale": "Oranges", "format": "{:,d}"},
        },
        "data": {sample: row},
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--vcf", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--out", required=True, help="Output JSON path")
    p.add_argument("--mqc-out",
                   help="Optional path for a MultiQC custom-content "
                        "_mqc.json sidecar (one section per file, merged "
                        "by sample at report time)")
    p.add_argument("--strict", action="store_true",
                   help="Exit 1 if any hard check fails")
    args = p.parse_args()

    result = run_qc(args.vcf, args.name)

    with open(args.out, "w") as fh:
        json.dump(result, fh, indent=2)

    if args.mqc_out:
        with open(args.mqc_out, "w") as fh:
            json.dump(build_mqc_payload(result), fh, indent=2)

    if result["errors"]:
        print(f"[qc_validate] HARD FAILURES for {args.name}:", file=sys.stderr)
        for e in result["errors"]:
            print(f"  - {e}", file=sys.stderr)
    if result["warnings"]:
        print(f"[qc_validate] WARNINGS for {args.name}:", file=sys.stderr)
        for w in result["warnings"]:
            print(f"  - {w}", file=sys.stderr)

    if args.strict and result["errors"]:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
