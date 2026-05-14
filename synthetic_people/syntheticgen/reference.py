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

import gzip
import shutil
import sys
import urllib.request
from pathlib import Path
from typing import Any

from .builds import BUILDS


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


def fetch_reference_fasta(cache_dir: Path, build: str) -> Path:
    """Ensure the reference FASTA + ``.fai`` index are cached on disk.

    Mirrors the ``clinvar.fetch_clinvar`` pattern: first run downloads
    the FASTA, decompresses it, and runs ``pysam.faidx``; subsequent
    runs find the cached file and return immediately. Returns the
    resolved ``.fa`` path.

    Layout under ``cache_dir``:

    - ``reference/<build>.fa``      — uncompressed FASTA
    - ``reference/<build>.fa.fai``  — pysam index

    Cache contract: a present-and-indexed ``.fa`` is treated as
    complete and is NOT re-downloaded. Partial downloads land at
    ``.fa.part`` and are renamed only on completion, so an
    interrupted run never leaves a half-written FASTA that the
    next run mistakes for valid. A missing ``.fai`` re-indexes
    in place without re-downloading.

    Disk footprint at GRCh38: ~900 MB peak during download
    (``.fa.gz.part`` + decompressing target) → ~3.1 GB steady
    state (``.fa`` + 50 KB ``.fai``). The ``.gz`` is deleted
    after decompression so only the indexed FASTA persists.

    Raises ``ImportError`` if pysam isn't installed (the index
    step needs it). Raises ``ValueError`` if the build isn't in
    the ``BUILDS`` table.
    """
    if build not in BUILDS:
        raise ValueError(
            f"unknown build {build!r}; supported: {sorted(BUILDS)}",
        )
    url = BUILDS[build].get("reference_fasta_url")
    if not url:
        raise ValueError(
            f"build {build!r} has no reference_fasta_url in BUILDS",
        )

    ref_dir = Path(cache_dir) / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)
    fa_path = ref_dir / f"{build}.fa"
    fai_path = fa_path.with_suffix(fa_path.suffix + ".fai")

    if fa_path.is_file() and fai_path.is_file():
        return fa_path

    if not fa_path.is_file():
        # Download .fa.gz, decompress to .fa.part, atomically rename.
        # Two-stage rename avoids a half-decompressed file masquerading
        # as a valid one if the process is killed mid-decompress.
        if not have_pysam():
            raise ImportError(
                "pysam is required to index the downloaded FASTA. "
                "Install with `pip install pysam` (optional dep).",
            )
        gz_part = ref_dir / f"{build}.fa.gz.part"
        gz_final = ref_dir / f"{build}.fa.gz"
        fa_part = ref_dir / f"{build}.fa.part"
        # PR #86 review (Copilot): if a previous run died mid-
        # decompress we may have a stale ``.gz`` on disk that
        # ``Path.rename`` won't overwrite on Windows. Clear stale
        # sidecars before the download starts so the rename below
        # is unambiguous on every platform.
        gz_final.unlink(missing_ok=True)
        fa_part.unlink(missing_ok=True)
        _download(url, gz_part)
        # Use ``.replace`` rather than ``.rename`` for atomic-publish
        # semantics on Windows + a race-safe overwrite on POSIX.
        gz_part.replace(gz_final)
        print(
            f"  decompressing {gz_final.name} → {fa_path.name} "
            f"(~3 GB; ~30 s on a fast SSD)",
            file=sys.stderr,
        )
        with gzip.open(gz_final, "rb") as src, open(fa_part, "wb") as dst:
            shutil.copyfileobj(src, dst, length=1 << 20)
        fa_part.replace(fa_path)
        # Remove the .gz now that the decompressed copy is on disk.
        gz_final.unlink(missing_ok=True)

    if not fai_path.is_file():
        print(
            f"  indexing {fa_path.name} (~60 s)",
            file=sys.stderr,
        )
        if not have_pysam():
            raise ImportError(
                "pysam is required for FASTA indexing. Install with "
                "`pip install pysam`.",
            )
        import pysam
        pysam.faidx(str(fa_path))

    return fa_path


def _download(url: str, dest: Path) -> None:
    """Stream-download ``url`` to ``dest`` with a progress log.

    Mirrors ``clinvar._download`` so the cli's two cache-fetch
    surfaces (ClinVar + reference FASTA) print consistently.
    """
    print(f"  downloading {url}", file=sys.stderr)
    with urllib.request.urlopen(url) as resp:
        total = int(resp.getheader("Content-Length") or 0)
        downloaded = 0
        last_report = 0
        with open(dest, "wb") as fh:
            while True:
                chunk = resp.read(1 << 20)  # 1 MiB
                if not chunk:
                    break
                fh.write(chunk)
                downloaded += len(chunk)
                if total and downloaded - last_report >= 50 * (1 << 20):
                    pct = downloaded * 100 / total
                    print(
                        f"    {downloaded / 1e6:7.1f} / "
                        f"{total / 1e6:.1f} MB ({pct:.0f}%)",
                        file=sys.stderr,
                    )
                    last_report = downloaded


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

    Returns one of ``A``/``C``/``G``/``T``/``N`` — non-canonical
    bases (IUPAC ambiguity codes ``R``/``Y``/``W``/``S``/``K``/``M``
    /``B``/``D``/``H``/``V``, masked-region ``N``/``n``, or any
    other byte the FASTA might contain) are normalised to ``N`` so
    callers can rely on a 5-character return alphabet.

    Why: ``choose_alt`` in the cli's Ti/Tv calibrator only knows how
    to weight transitions/transversions for the four canonical
    bases — any non-canonical input makes it return ``None``, which
    trips the producer's ``assert alt is not None``. Returning ``N``
    here funnels every non-canonical base through ``_pick_ref``'s
    existing rng fallback, preserving the producer's downstream
    invariants.

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
            upper = base.upper()
            # Canonical-only output; any IUPAC ambiguity or masked
            # base becomes ``N`` for callers' alphabet contract.
            return upper if upper in _CANONICAL_BASES else "N"
    return "N"


_CANONICAL_BASES = frozenset("ACGT")


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
