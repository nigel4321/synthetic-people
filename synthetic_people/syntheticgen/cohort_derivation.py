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


def derive_persons_batch(cohort_bcf_paths: list,
                         sample_ids: list) -> dict:
    """Stream per-person records for a batch of sample_ids in one
    bcftools subprocess per chromosome.

    Phase 5g: replaces the per-person ``derive_person_records``
    fan-out for callers that have all their target sample IDs up
    front (the streamed cohort path). The single-person form pays
    a full multi-sample BCF decode on every (person, chrom) pair —
    on a measured run at ``n=3000 × 22 chroms`` that was 66,000
    bcftools subprocesses each scanning a few hundred MB of cohort
    BCF, taking ~45 s/person and projecting to ~38 hours.

    The batched form runs::

        bcftools query -s s1,s2,...,sB -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO[\\t%GT]\\n'

    once per chrom — the ``[\\t%GT]`` template expands to one
    tab-separated GT per sample in ``-s`` order, so a single decode
    pass yields all B columns. The parser dispatches each row's
    GTs into ``per_person[sample_id]`` lists, dropping hom-ref
    (``0|0`` / ``0/0``) the same way ``derive_person_records``
    does (missing ``./.`` / ``.|.`` is kept). At ``B=20`` that
    drops the per-chrom invocation count by ~20×; at ``B=100`` by
    ~100×.

    Returns ``{sample_id: list_of_record_dicts}`` for every sid in
    ``sample_ids``. Each record dict has the same shape as
    ``derive_person_records``'s output, so downstream callers
    (``write_person_vcf``, etc.) consume either source identically.

    Memory bound: the returned dict holds every record for every
    listed sample for every chrom, simultaneously. Callers tuning
    against host RAM should batch their sample list so that
    ``B × E[per-person record list size]`` fits comfortably; at
    ``n=3000`` × full-22-chrom the per-person list is a few hundred
    MB to ~1 GB of parsed dicts, so ``B`` should typically stay
    under ~50 on 32 GB hosts.
    """
    if not sample_ids:
        return {}

    per_person: dict = {sid: [] for sid in sample_ids}
    sids_arg = ",".join(sample_ids)
    fmt = ("%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO"
           + "[\t%GT]\n")
    n = len(sample_ids)
    expected_field_count = 6 + n

    for bcf_path in cohort_bcf_paths:
        cmd = ["bcftools", "query", "-s", sids_arg,
               "-f", fmt, str(bcf_path)]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            if not line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < expected_field_count:
                # Malformed row (truncated, missing GT columns) —
                # skip so a single bad record doesn't poison the
                # entire batch.
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
                # Reference-only row — `bcftools query` doesn't filter
                # GT="ref" the way the per-person pipeline does, so we
                # check here instead. (It does happen for sites where
                # every sample is hom-ref; we drop them outright.)
                continue

            info_extras = _parse_info(info)
            # Build immutable shared template fields once per row;
            # per-person dicts share the alts/afs lists. The
            # downstream writer reads but never mutates them, so the
            # sharing is safe — and at n=3000 it cuts row-template
            # allocation cost dramatically.
            afs = [None] * len(alts)

            # M13.5: MT records must reach every per-person record
            # list regardless of the original simulator GT so the
            # lineage-clonality override in ``write_person_vcf`` can
            # rewrite each person's MT GT. Without this carve-out a
            # lineage carrier whose original simulator GT was hom-
            # ref would miss the MT record entirely, breaking the
            # "same-lineage → same MT record set" contract.
            is_mt = chrom in ("MT", "M", "chrMT", "chrM")
            for sid, gt in zip(sample_ids, gts):
                # `bcftools view -e 'GT="ref"'` matches 0/0 and 0|0
                # only; missing genotypes (./. / .|.) survive. The
                # per-person fan-out preserved that semantics, so
                # mirror it here.
                if not is_mt and (
                    not gt or gt == "0|0" or gt == "0/0"
                ):
                    continue
                rec = {
                    "chrom": chrom,
                    "pos": pos,
                    "id": vid,
                    "ref": ref,
                    "alts": alts,
                    "afs": afs,
                    "gt": gt,
                }
                if info_extras:
                    rec.update(info_extras)
                per_person[sid].append(rec)

        # Drain stderr and check returncode after stdout has been
        # consumed; bcftools occasionally writes warnings even on
        # success (mismatched mu in stdpopsim runs has surfaced in
        # this codebase before), so we keep the buffer for the
        # error path.
        proc.stdout.close()
        stderr_buf = ""
        if proc.stderr is not None:
            stderr_buf = proc.stderr.read()
            proc.stderr.close()
        proc.wait()
        if proc.returncode != 0:
            raise RuntimeError(
                f"bcftools query -s failed on {bcf_path} "
                f"(batch size={n}, exit {proc.returncode}): "
                f"{stderr_buf.strip()[:500] or '(no stderr)'}"
            )

    return per_person


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
    # M13.5: don't drop hom-ref MT records — the lineage-clonality
    # override in ``write_person_vcf`` needs every MT site to be
    # present in every person's per-person record list, regardless
    # of the original simulator GT. Without this carve-out, a
    # lineage carrier whose original simulator GT happened to be
    # hom-ref would never see the MT record, breaking the
    # "same-lineage → same MT record set" contract.
    drop_ref = subprocess.Popen(
        ["bcftools", "view", "-e",
         'GT="ref" && CHROM!="MT" && CHROM!="M"'
         ' && CHROM!="chrMT" && CHROM!="chrM"'],
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
