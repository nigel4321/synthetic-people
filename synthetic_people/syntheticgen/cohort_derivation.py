"""Per-person record derivation from cohort BCFs (Phase 5b2).

Replaces ``person_records_from_cohort`` for callers that source their
cohort sites from disk-backed per-chromosome BCFs (Phase 5b1) rather
than from an in-memory list. The output dict shape is identical so
``write_person_vcf`` can consume either.

Each call to :func:`derive_person_records` spawns one
``bcftools view`` subprocess per cohort chromosome BCF, filtered to a
single sample and to records where that sample carries at least one
alt allele. Records are parsed from the resulting VCF text and
materialised as a list of per-person record dicts in genome order
(``bcftools view`` already emits sorted output per file; the caller
iterates chromosomes in the requested order).

Memory is bounded by the largest per-person record list — at the
realistic ~10% non-hom-ref-per-sample fraction, that's a few
thousand records per person × a few hundred bytes ≈ low single-digit
MB even at n=100 000. The N-people fan-out runs the per-person
derivation in workers; since the cohort BCFs are read-only, multiple
``bcftools view`` subprocesses against the same BCF run cleanly in
parallel.
"""

from __future__ import annotations

import subprocess
from pathlib import Path


# Fields the per-person record dict carries when written by
# write_person_vcf. The cohort_derivation parser pulls them out of
# the INFO column on each BCF row and surfaces only the ones that
# are actually present (so a record without CLNSIG won't carry an
# empty `clnsig` key — same pattern as person_records_from_cohort).
_INFO_FIELDS_TO_CARRY = (
    "CLNSIG", "CLNDN",
    "COSMIC_ID", "COSMIC_GENE",
    "SVTYPE", "SVLEN", "END", "CIPOS",
)
_INFO_KEY_REMAP = {
    "CLNSIG": "clnsig",
    "CLNDN": "clndn",
    "COSMIC_ID": "cosmic_id",
    "COSMIC_GENE": "cosmic_gene",
    "SVTYPE": "svtype",
    "SVLEN": "svlen",
    "END": "end",
    "CIPOS": "cipos",
}


def _parse_info(info: str) -> dict:
    """Split an INFO field into a key→value dict, parsing only the
    fields ``write_person_vcf`` cares about."""
    out: dict = {}
    if not info or info == ".":
        return out
    for kv in info.split(";"):
        if "=" not in kv:
            # Flag-only INFO entries (e.g. HIGHLIGHT) — not in the
            # set we surface to per-person records.
            continue
        key, value = kv.split("=", 1)
        if key not in _INFO_FIELDS_TO_CARRY:
            continue
        dest = _INFO_KEY_REMAP[key]
        if key == "SVLEN" or key == "END":
            try:
                out[dest] = int(value)
            except ValueError:
                continue
        elif key == "CIPOS":
            try:
                lo, hi = value.split(",", 1)
                out[dest] = (int(lo), int(hi))
            except (ValueError, TypeError):
                continue
        else:
            out[dest] = value
    return out


def _parse_record(line: str) -> dict | None:
    """Parse a single VCF data row into a per-person record dict.

    Returns ``None`` for malformed rows (missing FORMAT block, sample
    column with the wrong number of fields, etc.) — the caller skips
    those rather than crashing the whole derivation.
    """
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 10:
        return None
    chrom, pos_s, vid, ref, alt, _qual, _filt, info, fmt = fields[:9]
    sample = fields[9]
    try:
        pos = int(pos_s)
    except ValueError:
        return None
    alts = alt.split(",") if alt and alt != "." else []
    if not alts:
        # Reference-only row — the bcftools `-i 'GT="alt"'` filter
        # should have dropped these, but guard against the edge case.
        return None
    fmt_keys = fmt.split(":")
    sample_vals = sample.split(":")
    try:
        gt_idx = fmt_keys.index("GT")
    except ValueError:
        return None
    if gt_idx >= len(sample_vals):
        return None
    gt = sample_vals[gt_idx]
    rec: dict = {
        "chrom": chrom,
        "pos": pos,
        "id": vid if vid else ".",
        "ref": ref,
        "alts": alts,
        # afs aren't carried back into the per-person record — the
        # writer recomputes per-record AC/AN/AF from the per-sample
        # genotype anyway, so any cohort-level AF value would just
        # be dropped on the floor.
        "afs": [None] * len(alts),
        "gt": gt,
    }
    rec.update(_parse_info(info))
    return rec


def derive_person_records(cohort_bcf_paths: list,
                          sample_id: str) -> list:
    """Stream per-person non-hom-ref records out of cohort BCFs.

    ``cohort_bcf_paths`` should be in genome order — typically the
    same list the manifest's ``cohort_bcfs[]`` field carries. Each
    BCF is queried via a two-step bcftools pipeline::

        bcftools view -s <sample_id> <bcf> | bcftools view -e 'GT="ref"'

    The ``-s`` first restricts to a single sample column; the second
    pass filters out records where that sample's now-single GT is
    homozygous-reference. ``GT="ref"`` matches ``0|0`` / ``0/0`` only
    — missing (``./.``) survives the filter, matching the keep-missing
    semantics of ``person_records_from_cohort``.

    Why pipelined rather than ``bcftools view -s SAMPLE -e 'GT="ref"'``:
    the bcftools filter expression evaluates *before* the sample
    subset is applied, against the multi-sample GT. ``GT="ref"`` on a
    multi-sample record returns true if *any* sample is hom-ref, so a
    single-call form drops every record that any sample has 0|0 at.
    The pipelined form filters on the post-subset single-sample GT,
    which is what we want.

    Returns a list of per-person record dicts shaped exactly like
    ``person_records_from_cohort``'s output, including any
    overlay-supplied INFO metadata (CLNSIG / COSMIC / SVTYPE).
    """
    records: list = []
    for bcf_path in cohort_bcf_paths:
        records.extend(_derive_from_one_bcf(Path(bcf_path), sample_id))
    return records


def _derive_from_one_bcf(bcf_path: Path, sample_id: str) -> list:
    """Spawn the ``bcftools view -s … | bcftools view -e 'GT="ref"'``
    pipeline and parse its output.

    Errors from bcftools surface as ``RuntimeError`` with the
    captured stderr — the caller's per-person worker should let
    these propagate so the pipeline fails fast on a missing or
    corrupted cohort BCF rather than silently producing an
    incomplete per-person VCF.
    """
    subset = subprocess.Popen(
        ["bcftools", "view", "-s", sample_id, str(bcf_path)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    drop_ref = subprocess.Popen(
        ["bcftools", "view", "-e", 'GT="ref"'],
        stdin=subset.stdout,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True,
    )
    # Closing the parent's read end of the subset pipe lets bcftools
    # propagate SIGPIPE if drop_ref exits early (e.g. on malformed
    # input).
    assert subset.stdout is not None
    subset.stdout.close()
    out: list = []
    assert drop_ref.stdout is not None
    for line in drop_ref.stdout:
        if not line or line.startswith("#"):
            continue
        rec = _parse_record(line)
        if rec is not None:
            out.append(rec)
    drop_ref.stdout.close()
    drop_ref_stderr = ""
    if drop_ref.stderr is not None:
        drop_ref_stderr = drop_ref.stderr.read()
        drop_ref.stderr.close()
    drop_ref.wait()
    subset_stderr = ""
    if subset.stderr is not None:
        subset_stderr = subset.stderr.read().decode(
            "utf-8", errors="replace")
        subset.stderr.close()
    subset.wait()
    if subset.returncode != 0:
        raise RuntimeError(
            f"bcftools view -s failed on {bcf_path} "
            f"(sample={sample_id!r}, exit {subset.returncode}): "
            f"{subset_stderr.strip()[:500] or '(no stderr)'}"
        )
    if drop_ref.returncode != 0:
        raise RuntimeError(
            f"bcftools view -e 'GT=ref' failed on {bcf_path} "
            f"(sample={sample_id!r}, exit {drop_ref.returncode}): "
            f"{drop_ref_stderr.strip()[:500] or '(no stderr)'}"
        )
    return out
