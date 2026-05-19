"""Single-pass cohort-to-person dispatch — Phase 5g.4.

Replaces the per-batch ``derive_persons_batch`` scan loop for callers
that need to fan out across many persons. The current per-batch
algorithm scans every cohort BCF once per batch; at large ``n``
that's quadratic in ``n`` because cohort BCF size is itself linear
in ``n`` (cohort BCF total scanned ≈ ``n / B × cohort_bcf_size``).

This module's :func:`dispatch_cohort_to_staging` scans each cohort
BCF exactly **once**, writing per-person records to per-sample
staging files as the scan proceeds. A downstream consumer reads its
own sample's staging file back into the same record-dict shape
``derive_persons_batch`` produces (see :func:`read_person_staging`).

The contract is **byte-equivalence**: feeding the dispatch+read path
the same cohort BCFs and sample IDs as ``derive_persons_batch`` must
produce identical per-person record lists. This is the gating
property for the Phase 5g.4 rollout; the validation tests in
``tests/test_cohort_dispatch.py`` enforce it.

Implementation hazards (each is a quality-degradation surface if
mishandled — see DATA_QUALITY_ASSESSMENT.md §6.2):

* **M13.5 MT carve-out.** Every MT record must reach every sample's
  staging file regardless of original simulator GT, so the lineage
  clonality override in ``write_person_vcf`` sees the same record
  set for every member of a lineage. Implemented by skipping the
  hom-ref GT filter on ``MT`` / ``M`` / ``chrMT`` / ``chrM`` chroms,
  matching ``cohort_derivation.derive_persons_batch``.
* **``afs=[None]`` shape.** The streamed coalescent path's
  ``write_person_vcf`` MT lineage-carrier fallback (AF=0.1) depends
  on the ``afs=[None] * len(alts)`` shape. Preserved on read.
* **Full ``_INFO_FIELDS_TO_CARRY`` set.** Round-tripped via JSON in
  the staging row's INFO column. ``cipos`` is re-tupled on read so
  the writer sees ``tuple[int, int]`` (matching the in-memory
  derivation shape).
* **Record ordering.** Within each cohort BCF rows arrive in
  bcftools-sort order (chrom + pos). Across BCFs the order follows
  ``cohort_bcf_paths``. Both are preserved through staging because
  records are appended in the order they're scanned and the writer
  consumer reads sequentially.

Out of scope for this module (tracked separately in
PERFORMANCE_PLAN.md §5g.4):

* Per-chromosome chunked staging (reduces peak disk at large ``n``).
* FD-limit window rotation for ``n > ulimit -n``.
* Per-chrom resume mid-dispatch state machine.
* Compression of staging fragments.

Those land in PR-A.2 / PR-B once the byte-equivalence contract has
bake time at n ≤ 100.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path


# Same set + remap as ``cohort_derivation`` — these MUST stay in
# lockstep so the dispatch staging round-trip preserves every INFO
# field the in-memory derivation does.
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

# JSON round-tripping turns tuples into lists; re-tuple these keys on
# read so the shape matches ``cohort_derivation``'s record dicts.
_TUPLE_KEYS = ("cipos",)


def _parse_info(info: str) -> dict:
    """Extract the carry-forward keys from a VCF INFO column.

    Mirrors ``cohort_derivation._parse_info`` exactly. Kept in this
    module rather than re-imported because the byte-equivalence
    contract requires the two parsers to behave identically — a
    refactor that drifts the dispatch parser away from the
    in-memory parser would silently break the contract.
    """
    out: dict = {}
    if not info or info == ".":
        return out
    for kv in info.split(";"):
        if "=" not in kv:
            continue
        key, value = kv.split("=", 1)
        if key not in _INFO_FIELDS_TO_CARRY:
            continue
        dest = _INFO_KEY_REMAP[key]
        if key in ("SVLEN", "END"):
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


def _staging_path(staging_dir: Path, sample_id: str) -> Path:
    # Sample IDs in this codebase come from the cohort meta and are
    # generator-controlled (alphanumeric NA/HG-prefixed identifiers
    # like ``NA12345``). Scrub the two characters that would cause
    # filesystem trouble defensively.
    safe = sample_id.replace("/", "_").replace("\0", "_")
    return staging_dir / f"person_{safe}.tsv"


def dispatch_cohort_to_staging(
    cohort_bcf_paths: list,
    sample_ids: list,
    staging_dir: Path,
) -> dict:
    """Scan cohort BCFs once, write each person's records to its own
    staging file.

    Returns ``{sid: staging_path}`` for every sample_id in
    ``sample_ids``. The staging files contain one row per per-person
    record in the TSV shape :func:`read_person_staging` parses back
    to the record-dict format ``derive_persons_batch`` produces.

    The cohort BCFs are scanned in the order given. Within each BCF,
    rows are processed in bcftools-iteration order (chrom + pos).
    Per-person record ordering therefore matches what
    ``derive_persons_batch`` produces when passed the same list.

    ``staging_dir`` is created if it doesn't exist. Per-sample
    staging files are opened in ``w`` mode (overwrite). The caller
    is responsible for cleaning up the staging dir after consumers
    finish; this function only writes.
    """
    staging_dir.mkdir(parents=True, exist_ok=True)
    paths = {sid: _staging_path(staging_dir, sid) for sid in sample_ids}
    if not sample_ids:
        return paths

    # Open every staging writer up front. The caller is expected to
    # provide a sample_ids list that fits the host's FD ulimit; the
    # chunked window-rotation variant for ``n > ulimit -n`` is a
    # follow-up (see PERFORMANCE_PLAN §5g.4).
    writers = {
        sid: paths[sid].open("w", encoding="utf-8")
        for sid in sample_ids
    }
    sids_arg = ",".join(sample_ids)
    fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO" + "[\t%GT]\n"
    n = len(sample_ids)
    expected_field_count = 6 + n

    try:
        for bcf_path in cohort_bcf_paths:
            cmd = ["bcftools", "query", "-s", sids_arg,
                   "-f", fmt, str(bcf_path)]
            proc = subprocess.Popen(
                cmd, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, text=True,
            )
            assert proc.stdout is not None
            for line in proc.stdout:
                if not line:
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < expected_field_count:
                    # Malformed row — skip so a single bad record
                    # doesn't poison the whole dispatch pass.
                    continue
                chrom = fields[0]
                try:
                    pos = int(fields[1])
                except ValueError:
                    continue
                vid = fields[2] if fields[2] else "."
                ref = fields[3]
                alt = fields[4]
                info = fields[5]
                gts = fields[6:6 + n]
                alts = alt.split(",") if alt and alt != "." else []
                if not alts:
                    # Reference-only row — drop outright. Matches
                    # derive_persons_batch.
                    continue
                info_extras = _parse_info(info)
                info_json = (
                    json.dumps(info_extras, sort_keys=True)
                    if info_extras else ""
                )
                alts_csv = ",".join(alts)
                is_mt = chrom in ("MT", "M", "chrMT", "chrM")
                row_prefix = (
                    f"{chrom}\t{pos}\t{vid}\t{ref}\t{alts_csv}\t"
                )
                row_suffix = f"\t{info_json}\n"
                for sid, gt in zip(sample_ids, gts):
                    # M13.5 carve-out: every MT record reaches every
                    # sample's staging regardless of GT, so the
                    # write-time lineage override sees the same MT
                    # record set for every member of a lineage.
                    if not is_mt and (
                        not gt or gt == "0|0" or gt == "0/0"
                    ):
                        continue
                    writers[sid].write(row_prefix + gt + row_suffix)
            proc.stdout.close()
            stderr_buf = ""
            if proc.stderr is not None:
                stderr_buf = proc.stderr.read()
                proc.stderr.close()
            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError(
                    f"bcftools query -s failed on {bcf_path} "
                    f"(N={n}, exit {proc.returncode}): "
                    f"{stderr_buf.strip()[:500] or '(no stderr)'}"
                )
    finally:
        for w in writers.values():
            w.close()

    return paths


def read_person_staging(staging_path: Path) -> list:
    """Parse a per-person staging file back to the record-dict list
    ``write_person_vcf`` consumes.

    The returned shape matches ``derive_persons_batch``'s per-person
    record list exactly: same keys, same value types, same
    ``afs=[None] * len(alts)`` shape, same ``cipos`` tuple type.
    """
    records: list = []
    if not staging_path.exists():
        return records
    with staging_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 7:
                continue
            chrom, pos_s, vid, ref, alts_csv, gt, info_json = parts
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            alts = alts_csv.split(",") if alts_csv else []
            if not alts:
                continue
            rec: dict = {
                "chrom": chrom,
                "pos": pos,
                "id": vid if vid else ".",
                "ref": ref,
                "alts": alts,
                "afs": [None] * len(alts),
                "gt": gt,
            }
            if info_json:
                extras = json.loads(info_json)
                for k in _TUPLE_KEYS:
                    if k in extras:
                        extras[k] = tuple(extras[k])
                rec.update(extras)
            records.append(rec)
    return records
