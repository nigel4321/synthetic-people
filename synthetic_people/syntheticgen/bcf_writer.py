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
import subprocess
from pathlib import Path

from .builds import BUILDS
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
                 sample_ids: list[str]):
        self.out_path = Path(out_path)
        self.build = build
        self.sample_ids = list(sample_ids)
        self.n_people = len(sample_ids)
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

        The site dict is the same shape produced by ``simulate_cohort``
        (chrom / pos / id / ref / alts / acs / gts plus optional
        overlay metadata). Only the cohort-level GT block is emitted
        here; per-call DP/GQ/AD are layered in at per-person derivation
        time, identical to today's writer.py behaviour.
        """
        if self._fh is None:
            raise RuntimeError("write_site called outside `with` block")
        gts = site.get("gts")
        if gts is None or len(gts) != self.n_people:
            raise ValueError(
                f"site at {site.get('chrom')}:{site.get('pos')} has "
                f"{len(gts) if gts is not None else 'no'} GTs; expected "
                f"{self.n_people}"
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
