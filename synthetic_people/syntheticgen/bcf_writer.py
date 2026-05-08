"""Multi-sample cohort BCF writer (Phase 5).

Streams a cohort's per-chromosome sites into a bgzipped BCF on disk so
the cohort genotype state never has to fit in RAM end-to-end. Pairs
with the chromosome-streaming variant of ``simulate_cohort``: each
chromosome's sites are passed through ``write_chrom_bcf`` and freed
before the next chromosome is simulated.

Why BCF (binary VCF) and not text VCF: ``bcftools view -s SAMPLE``
extracts a single person from a cohort BCF an order of magnitude
faster than the equivalent text-VCF query, and that extraction is the
per-person derivation step the rest of Phase 5 relies on.

Why subprocess into ``bcftools view -O b`` rather than pysam: pysam
isn't on the dependency tree (msprime / tskit / stdpopsim don't pull
it in), and adding it is a meaningful install-time cost — pysam ships
its own htslib build with C extensions. The existing writer.py
already streams text into ``bgzip -c`` via ``Popen``; this module
follows the same pattern, just with ``bcftools view -O b`` swapped
in. A future commit can switch to pysam if the dep tree changes.
"""

from __future__ import annotations

import io
import multiprocessing as mp
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from .builds import BUILDS
from .cohort_sites import dense_gts_from_carriers, dense_gts_from_carriers_slice
from .header import _ALT_SV, _FORMAT, _INFO_CORE, _INFO_SV


def build_cohort_header(build: str, sample_ids: list[str]) -> str:
    """Return the VCF header for a multi-sample cohort BCF.

    Same INFO / FORMAT / ALT declarations as the per-person header in
    ``header.build_header`` — keeping them in lockstep means a per-
    person VCF derived via ``bcftools view -s SAMPLE`` is
    indistinguishable from one written directly by ``writer.py``
    (modulo per-record FORMAT-tag ordering, which bcftools preserves).
    """
    info = BUILDS[build]
    assembly = info["assembly"]
    reference = info["reference"]
    contigs = info["contigs"]

    lines: list[str] = [
        "##fileformat=VCFv4.2",
        "##source=synthetic_people/generate_people.py",
        f"##reference={reference}",
    ]
    for chrom, length in contigs.items():
        lines.append(
            f"##contig=<ID={chrom},length={length},assembly={assembly}>"
        )
    lines.extend(_INFO_CORE)
    lines.extend(_INFO_SV)
    lines.extend(_ALT_SV)
    lines.extend(_FORMAT)
    chrom_line = (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        + "".join(f"\t{s}" for s in sample_ids)
    )
    lines.append(chrom_line)
    return "\n".join(lines) + "\n"


def _format_info(site: dict, n_people: int) -> str:
    """Build the per-record INFO field for a cohort BCF row.

    Cohort-level AC / AN / AF reflect the *true* genotypes across all
    samples, not what any single per-person VCF will eventually emit
    after sequencing-noise injection. The per-person derivation step
    re-computes these from each person's view of the record.
    """
    acs = site.get("acs") or [0] * len(site["alts"])
    an = 2 * n_people  # diploid cohort — N people × 2 haplotypes
    afs = [(ac / an) if an else 0.0 for ac in acs]
    parts: list[str] = [
        f"AC={','.join(str(c) for c in acs)}",
        f"AN={an}",
        f"AF={','.join(f'{f:.6f}' for f in afs)}",
    ]
    if site.get("clnsig") and site["clnsig"] != ".":
        parts.append(f"CLNSIG={site['clnsig']}")
    if site.get("clndn") and site["clndn"] != ".":
        parts.append(f"CLNDN={site['clndn']}")
    if site.get("cosmic_id"):
        parts.append(f"COSMIC_ID={site['cosmic_id']}")
    if site.get("cosmic_gene"):
        parts.append(f"COSMIC_GENE={site['cosmic_gene']}")
    if site.get("svtype"):
        parts.append(f"SVTYPE={site['svtype']}")
        if site.get("svlen") is not None:
            parts.append(f"SVLEN={site['svlen']}")
        if site.get("end") is not None:
            parts.append(f"END={site['end']}")
        if site.get("cipos"):
            lo, hi = site["cipos"]
            parts.append(f"CIPOS={lo},{hi}")
    return ";".join(parts)


class CohortBcfWriter:
    """Streaming writer for a single cohort BCF.

    Spawns a ``bcftools view -O b`` subprocess and pipes a text VCF
    into its stdin. The subprocess writes the bgzipped binary form to
    ``out_path`` directly, so we never materialise an intermediate
    text VCF on disk.

    Use as a context manager::

        with CohortBcfWriter(out_path, build, sample_ids) as w:
            for site in chrom_sites:
                w.write_site(site)
        # On context exit, bcftools is allowed to flush, the BCF is
        # closed, and we run `bcftools index` to drop a CSI index next
        # to it.

    The CSI index (rather than tabix's TBI) is required because BCF
    files cannot be tabix-indexed — only the sister text-bgzip format.
    CSI is htslib-native and ``bcftools view -s SAMPLE -r REGION`` works
    against it identically.
    """

    def __init__(self, out_path: Path, build: str,
                 sample_ids: list[str],
                 sample_slice: tuple[int, int] | None = None,
                 cohort_size: int | None = None):
        """Construct a cohort BCF writer.

        ``sample_ids`` is the *header* sample list — what shows up
        as the ``#CHROM`` header columns in the output BCF. By
        default, ``write_site`` expands carriers across all
        ``len(sample_ids)`` people (the original full-cohort write).

        Phase 5e Phase A — sample-slice mode:

        - ``sample_slice`` is an optional ``(slice_lo, slice_hi)``
          person-index tuple. When set, ``sample_ids`` should be
          the slice's sample names, and ``cohort_size`` should be
          the *full* cohort size (so the per-site INFO ``AN`` /
          ``AF`` reflect the cohort, not just this slice). At write
          time, only the slice's GTs are formatted into the sample
          block.
        - When ``sample_slice`` is unset, ``cohort_size`` defaults
          to ``len(sample_ids)`` and the full-cohort path runs
          unchanged.

        AC counts come from ``site["acs"]`` already, which the
        cohort simulator computes across the whole cohort, so they
        don't need slice adjustment — the slice-write just trims
        the per-sample columns.
        """
        self.out_path = Path(out_path)
        self.build = build
        self.sample_ids = list(sample_ids)
        self.n_header_samples = len(self.sample_ids)
        self._sample_slice = sample_slice
        if sample_slice is not None:
            slice_lo, slice_hi = sample_slice
            if slice_lo < 0 or slice_hi < slice_lo:
                raise ValueError(
                    f"invalid sample_slice {sample_slice!r}: "
                    f"requires 0 <= lo <= hi"
                )
            slice_n = slice_hi - slice_lo
            if slice_n != self.n_header_samples:
                raise ValueError(
                    f"sample_slice covers {slice_n} persons but "
                    f"len(sample_ids)={self.n_header_samples}; pass "
                    f"the slice's sample_ids list as sample_ids"
                )
            if cohort_size is None:
                raise ValueError(
                    "cohort_size is required when sample_slice is "
                    "set (so per-site AN/AF reflect the full cohort)"
                )
            self.n_people = cohort_size
        else:
            if cohort_size is not None and \
                    cohort_size != self.n_header_samples:
                raise ValueError(
                    f"cohort_size={cohort_size} but no sample_slice "
                    f"and len(sample_ids)={self.n_header_samples}; "
                    f"either drop cohort_size or pass sample_slice"
                )
            self.n_people = self.n_header_samples
        self._proc: subprocess.Popen | None = None
        self._fh: io.TextIOWrapper | None = None
        self._out_fh = None

    def __enter__(self) -> "CohortBcfWriter":
        self.out_path.parent.mkdir(parents=True, exist_ok=True)
        # Open the destination ourselves and hand it to bcftools as
        # stdout, mirroring the writer.py bgzip pattern.
        self._out_fh = open(self.out_path, "wb")
        self._proc = subprocess.Popen(
            ["bcftools", "view", "-O", "b", "-"],
            stdin=subprocess.PIPE,
            stdout=self._out_fh,
            stderr=subprocess.PIPE,
        )
        assert self._proc.stdin is not None
        self._fh = io.TextIOWrapper(
            self._proc.stdin, encoding="utf-8",
            write_through=True, line_buffering=False,
        )
        self._fh.write(build_cohort_header(self.build, self.sample_ids))
        return self

    def write_site(self, site: dict) -> None:
        """Write one cohort-site row.

        Phase 5c: accepts either ``site["carriers"]`` (sparse
        ``[(haplotype_idx, allele_idx), ...]`` — the canonical form
        emitted by simulate_cohort_iter) or ``site["gts"]`` (legacy
        dense list of ``"a|b"`` strings — kept as a fallback for
        tests and any caller that hasn't migrated). Carriers expand
        to dense GT strings at write time.

        Per-call DP/GQ/AD are layered in at per-person derivation
        time, identical to today's writer.py behaviour.
        """
        if self._fh is None:
            raise RuntimeError("write_site called outside `with` block")

        if "carriers" in site:
            if self._sample_slice is not None:
                slice_lo, slice_hi = self._sample_slice
                gts = dense_gts_from_carriers_slice(
                    site["carriers"], slice_lo, slice_hi)
            else:
                gts = dense_gts_from_carriers(
                    site["carriers"], self.n_people)
        elif site.get("gts") is not None:
            gts = site["gts"]
            # Dense-GT input is taken at face value; the caller is
            # expected to slice the list themselves before passing
            # it in (this branch is a fallback for legacy callers
            # and tests that already work in dense form).
            if len(gts) != self.n_header_samples:
                raise ValueError(
                    f"site at {site.get('chrom')}:{site.get('pos')} has "
                    f"{len(gts)} GTs; expected {self.n_header_samples}"
                )
        else:
            raise ValueError(
                f"site at {site.get('chrom')}:{site.get('pos')} has "
                f"neither `carriers` nor `gts`; expected one"
            )

        # FORMAT carries GT only at the cohort level — DP/GQ/AD are
        # per-person and get drawn during per-person derivation.
        sample_block = "\t".join(gts)
        info = _format_info(site, self.n_people)
        line = "\t".join([
            site["chrom"], str(site["pos"]),
            site.get("id") or ".",
            site["ref"], ",".join(site["alts"]),
            "100", "PASS", info,
            "GT", sample_block,
        ]) + "\n"
        self._fh.write(line)

    def write_sites(self, sites) -> int:
        """Bulk-write an iterable of sites; returns the count written."""
        n = 0
        for s in sites:
            self.write_site(s)
            n += 1
        return n

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        if self._fh is not None:
            try:
                self._fh.flush()
                self._fh.detach()
            except Exception:
                pass
            self._fh = None
        if self._proc is not None and self._proc.stdin is not None:
            if not self._proc.stdin.closed:
                self._proc.stdin.close()
        rc = self._proc.wait() if self._proc else 0
        stderr_tail = ""
        if self._proc is not None and self._proc.stderr is not None:
            stderr_tail = self._proc.stderr.read().decode(
                "utf-8", errors="replace").strip()[-500:]
            self._proc.stderr.close()
        if self._out_fh is not None:
            self._out_fh.close()
            self._out_fh = None
        if rc != 0:
            # Surface the subprocess's own diagnostic so a malformed
            # site (e.g. wrong sample count) doesn't look like a phantom
            # writer failure.
            raise RuntimeError(
                f"bcftools view -O b exited {rc} writing {self.out_path}: "
                f"{stderr_tail or '(no stderr captured)'}"
            )
        # CSI index for region + sample queries downstream.
        subprocess.run(
            ["bcftools", "index", "-f", str(self.out_path)],
            check=True, capture_output=True,
        )


# ---------------------------------------------------------------------------
# Phase 5e Phase A — sample-slice parallel cohort BCF write
# ---------------------------------------------------------------------------
#
# The serial CohortBcfWriter above pipes a single 100s-of-MB text VCF
# through one ``bcftools view -O b`` subprocess. At ``n=3000 × 70 Mb``
# the per-chrom write was ~960 s — Python text formatting + bgzip
# encode are both single-threaded inside that subprocess and dominate
# the cohort-phase wall time.
#
# Parallel-write splits the per-site sample block across W workers
# along contiguous person slices. Each worker writes its own partial
# BCF (same site set, different sample columns); the parent then
# runs ``bcftools merge`` to join them into the final cohort BCF.
# Workers fork-COW share the parent's ``sites`` list — no per-worker
# copy of the records — by reading from a module-level
# ``_COHORT_WRITE_STATE`` dict the parent populates before forking.
# Per-worker private memory stays small (just the worker's bgzip
# pipe buffer + the slice's GT-string list per site).
#
# bcftools merge (not concat): partials are sample-disjoint with
# the same variant set, so merge collapses to a sample-column join
# at each (chrom, pos, ref, alt) key. Concat is the wrong tool —
# that's for region-disjoint inputs.

_COHORT_WRITE_STATE: dict = {}


def _partial_writer_target(slice_idx: int) -> None:
    """Worker entry point — writes one sample-slice partial BCF.

    Reads shared inputs (sites list, full sample list, build name,
    output paths, slice ranges) from ``_COHORT_WRITE_STATE`` so they
    are inherited via fork rather than pickled per task. Only the
    slice index is passed across the task boundary.
    """
    state = _COHORT_WRITE_STATE
    out_path = state["partial_paths"][slice_idx]
    slice_lo, slice_hi = state["slices"][slice_idx]
    full_sample_ids = state["sample_ids"]
    sites = state["sites"]
    slice_sample_ids = full_sample_ids[slice_lo:slice_hi]
    with CohortBcfWriter(
        out_path, state["build"], slice_sample_ids,
        sample_slice=(slice_lo, slice_hi),
        cohort_size=len(full_sample_ids),
    ) as w:
        w.write_sites(sites)


def _split_into_slices(n: int, workers: int) -> list[tuple[int, int]]:
    """Split ``n`` people into ``workers`` contiguous slices.

    Returns a list of ``(lo, hi)`` pairs covering ``[0, n)`` with
    ``hi - lo`` differing by at most 1 across slices. Empty slices
    (when ``workers > n``) are dropped — fewer partials get written
    than workers requested. Deterministic and stable across runs.
    """
    if n <= 0 or workers <= 0:
        return []
    base, rem = divmod(n, workers)
    slices: list[tuple[int, int]] = []
    cursor = 0
    for i in range(workers):
        size = base + (1 if i < rem else 0)
        if size == 0:
            continue
        slices.append((cursor, cursor + size))
        cursor += size
    return slices


def write_cohort_bcf_parallel(out_path: Path, build: str,
                              sample_ids: list[str],
                              sites: list,
                              workers: int) -> None:
    """Write a cohort BCF using ``workers`` parallel sample-slice
    writers + a final ``bcftools merge``.

    Falls back to the serial :class:`CohortBcfWriter` when
    ``workers <= 1`` or ``len(sample_ids) <= 1`` — at that scale
    the parallel orchestration costs more than it saves.

    Output equivalence: at any ``--workers``, the merged BCF
    contains the same sites in the same order with the same
    per-sample columns as a serial-write would. The merge step
    is deterministic given the partial inputs.

    Resume: the ``.partials/<bcf-stem>/`` directory only exists
    *during* the per-chrom write. On success it's cleaned up before
    the function returns, so a successfully-completed chromosome
    never has stale partials on disk for the next run to trip over.
    On failure the partials are left in place for postmortem.
    """
    n_samples = len(sample_ids)
    if workers <= 1 or n_samples <= 1:
        with CohortBcfWriter(out_path, build, sample_ids) as w:
            w.write_sites(sites)
        return

    slices = _split_into_slices(n_samples, workers)
    if len(slices) <= 1:
        # workers > n_samples — only one non-empty slice would
        # exist; serial write is simpler.
        with CohortBcfWriter(out_path, build, sample_ids) as w:
            w.write_sites(sites)
        return

    out_path = Path(out_path)
    partials_dir = out_path.parent / ".partials" / out_path.stem
    partials_dir.mkdir(parents=True, exist_ok=True)
    partial_paths = [
        partials_dir / f"slice_{i:03d}.bcf"
        for i in range(len(slices))
    ]

    # Stage shared inputs in module-level state so fork-spawned
    # workers see them via COW. Cleared in ``finally`` so a
    # workers-failed-mid-run doesn't leave the state populated.
    _COHORT_WRITE_STATE.update({
        "partial_paths": partial_paths,
        "slices": slices,
        "sample_ids": sample_ids,
        "sites": sites,
        "build": build,
    })
    try:
        ctx = mp.get_context("fork")
        with ProcessPoolExecutor(max_workers=len(slices),
                                 mp_context=ctx) as ex:
            futures = [
                ex.submit(_partial_writer_target, i)
                for i in range(len(slices))
            ]
            # ``fut.result()`` re-raises any worker exception so a
            # malformed site or bcftools-write failure surfaces with
            # its original traceback rather than as a silent partial-
            # missing failure during the merge step.
            for fut in futures:
                fut.result()
    finally:
        _COHORT_WRITE_STATE.clear()

    # Merge partials into the final cohort BCF. ``bcftools merge``
    # joins by ``(chrom, pos, ref, alt)`` key — every partial has the
    # same variants in the same order, so the join collapses to a
    # sample-column concatenation in the order the partials appear
    # on the command line. Pass them in slice order so the merged
    # sample columns align with ``sample_ids``.
    merge_cmd = [
        "bcftools", "merge", "-O", "b",
        "-o", str(out_path),
        *[str(p) for p in partial_paths],
    ]
    merge_proc = subprocess.run(merge_cmd, capture_output=True)
    if merge_proc.returncode != 0:
        stderr = merge_proc.stderr.decode("utf-8",
                                          errors="replace").strip()
        raise RuntimeError(
            f"bcftools merge exited {merge_proc.returncode} writing "
            f"{out_path}: {stderr[-500:] or '(no stderr)'}"
        )
    subprocess.run(
        ["bcftools", "index", "-f", str(out_path)],
        check=True, capture_output=True,
    )

    # Clean up partials + their indexes. The .partials/<stem>/
    # directory only existed for the duration of this write; remove
    # it (and the parent .partials/ if it's now empty) so a resume
    # doesn't see stale state.
    for p in partial_paths:
        for path in (p, Path(str(p) + ".csi")):
            if path.exists():
                path.unlink()
    try:
        partials_dir.rmdir()
        # Best-effort cleanup of the .partials/ container; only
        # remove if it's empty (other chroms may still be writing).
        if not any(partials_dir.parent.iterdir()):
            partials_dir.parent.rmdir()
    except OSError:
        pass
