"""Reference-FASTA loading for M12.

Until M12 the cli drew ``REF`` from ``rng.choice("ACGT")`` —
fabricated bases that didn't match the GRCh38 primary assembly.
Anything downstream that re-validated the VCF against the
reference (``bcftools norm --check-ref``, aligners, validators)
would fail on every record. Tier 1's REF-check gate was
designed precisely so that a passing run *after* M12 lands is
empirical evidence the wiring works.

This module lets the cli load a real FASTA and emit real REF
bases. The four producers (materialised path in
``coalescent._tree_sequence_to_sites``, streaming pass 1 in
``coalescent._tree_sequence_to_sites_meta``, streaming pass 2
in ``coalescent._stream_cohort_pass2``, admixture in
``admixture.py``) look up the REF at simulation time via
:func:`fetch_ref_base` instead of drawing from the cli's rng.

Falls back gracefully when no FASTA is provided — the old
fabricated-REF path is preserved for tests, smoke runs, and
quick development that doesn't need real reference content.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def have_pysam() -> bool:
    """Probe whether pysam is importable.

    Pysam is required only when ``--reference-fasta`` is used.
    The rest of the pipeline (and the test suite) doesn't depend
    on it; treat it as an optional install on the same footing
    as pyarrow.
    """
    try:
        import pysam  # noqa: F401
        return True
    except ImportError:
        return False


def load_fasta(path: Path) -> Any:
    """Open a FASTA via ``pysam.FastaFile``.

    Memory-mapped: resident-set is tiny (~50 MB for GRCh38 primary)
    regardless of FASTA size, because the kernel pages bases on
    demand. The returned handle is shareable across workers via
    the standard ``ProcessPoolExecutor`` fork-mmap path.

    Raises ``ImportError`` (clear cli message) if pysam isn't
    installed, ``FileNotFoundError`` if the path is missing, or
    ``ValueError`` if the FASTA isn't indexed (``.fai`` missing —
    pysam auto-creates one on open, but read-only filesystems
    can break that; the message points at ``samtools faidx``).
    """
    if not have_pysam():
        raise ImportError(
            "pysam is required for --reference-fasta. Install with "
            "`pip install pysam` (optional dep alongside pyarrow)."
        )
    import pysam
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            f"--reference-fasta: file not found: {path}",
        )
    try:
        return pysam.FastaFile(str(path))
    except (OSError, ValueError) as exc:
        raise ValueError(
            f"--reference-fasta: could not open {path} — {exc}. "
            f"If the .fai index is missing, run `samtools faidx "
            f"{path}` first.",
        ) from exc


def fetch_ref_base(fasta: Any, chrom: str, pos: int) -> str:
    """Look up the reference base at 1-based ``pos`` on ``chrom``.

    Returns the uppercase base (``A``/``C``/``G``/``T``/``N``).

    Handles the common chr-prefix mismatch between the cli's
    ``BUILDS`` naming (``"22"``) and UCSC/Ensembl FASTA naming
    (``"chr22"``) by trying both. The cli's ``BUILDS`` keys are
    unprefixed; many widely-distributed GRCh38 FASTAs use the
    ``chr`` prefix. We don't want users to have to massage either
    side to make them line up.

    Returns ``"N"`` if the chromosome doesn't exist in the FASTA
    or the position is out of range — defensive for malformed
    input rather than a hard error per-variant, since per-variant
    raises would tank a multi-hour run on one bad coordinate.
    Pre-flight validation in :func:`validate_fasta` catches the
    typical "wrong FASTA" cases at startup instead.
    """
    if fasta is None:
        return "N"
    # Try exact + cross-prefix to bridge BUILDS-vs-UCSC convention.
    names = [chrom]
    if chrom.startswith("chr"):
        names.append(chrom[3:])
    else:
        names.append(f"chr{chrom}")
    for name in names:
        try:
            base = fasta.fetch(name, pos - 1, pos)
        except (KeyError, IndexError, ValueError):
            continue
        if base:
            return base.upper()
    return "N"


def resolve_chrom_name(fasta: Any, chrom: str) -> str | None:
    """Return the FASTA-side name for ``chrom``, or ``None`` if
    neither prefix variant exists in the FASTA.

    Used by :func:`validate_fasta` to surface a clear error at
    startup if the user passes a FASTA whose chromosome names
    don't match the requested ``--chromosomes`` set, instead of
    silently producing ``N``-only output.
    """
    if fasta is None:
        return None
    available = set(fasta.references)
    if chrom in available:
        return chrom
    if chrom.startswith("chr"):
        alt = chrom[3:]
    else:
        alt = f"chr{chrom}"
    if alt in available:
        return alt
    return None


def validate_fasta(
    fasta: Any, chromosomes: list, chr_length_mb: float,
) -> None:
    """Pre-flight check: every requested chrom exists in the FASTA
    AND is long enough for the configured ``chr_length_mb``.

    Raises ``ValueError`` with a clear message on mismatch so the
    user can spot a "wrong FASTA" error in seconds instead of
    seeing a stream of N's hours into the run.
    """
    if fasta is None:
        return
    missing = []
    too_short = []
    for chrom in chromosomes:
        resolved = resolve_chrom_name(fasta, chrom)
        if resolved is None:
            missing.append(chrom)
            continue
        if chr_length_mb > 0:
            try:
                length_bp = fasta.get_reference_length(resolved)
            except (KeyError, ValueError):
                length_bp = 0
            if length_bp < chr_length_mb * 1_000_000:
                too_short.append(
                    f"{chrom} (FASTA has "
                    f"{length_bp / 1_000_000:.1f} Mb, "
                    f"need {chr_length_mb} Mb)",
                )
    errors = []
    if missing:
        errors.append(
            "missing chromosomes: " + ", ".join(missing)
            + f". Available: {sorted(fasta.references)[:10]}"
            + ("…" if len(fasta.references) > 10 else ""),
        )
    if too_short:
        errors.append(
            "chromosomes shorter than --chr-length-mb: "
            + ", ".join(too_short),
        )
    if errors:
        raise ValueError(
            "--reference-fasta validation failed:\n  - "
            + "\n  - ".join(errors),
        )
