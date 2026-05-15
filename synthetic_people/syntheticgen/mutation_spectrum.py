"""96-channel mutation spectrum diagnostic (DATA_QUALITY_ASSESSMENT.md
Tier 2 #5).

A SNV's mutation channel is its substitution type plus the trinucleotide
context (flanking base on each side). After folding purine-context calls
to pyrimidine-context via reverse-complement, there are exactly 6
substitution types {C>A, C>G, C>T, T>A, T>C, T>G} × 16 contexts
(4 left flanks × 4 right flanks) = **96 channels**.

The spectrum is the relative density across those 96 channels. Real
human germline / somatic spectra are dominated by clock-like SBS1
(C>T at NpCpG contexts) plus tissue-specific signatures. Today's
``BinaryMutationModel`` in ``coalescent.py`` produces a degenerate
spectrum — the simulator pins Ti/Tv globally but doesn't condition on
trinucleotide context, so post-M12 the cohort's spectrum is far from
SBS1. This module is the empirical gate for M14, which adds
context-aware μ via stdpopsim's mutation-rate tables. Pre-M14: spectrum
is flat-ish (Ti/Tv-pinned). Post-M14: should match SBS1.

Scope of this module: compute and emit the 96-channel counts +
fractions. **Comparison against a published reference signature (SBS1
cosine similarity) is deferred to a follow-up PR** so the reference
vector can be sourced carefully (see
``DATA_QUALITY_ASSESSMENT.md`` §A.4 Tier 2 #5 deferral note).
"""

from __future__ import annotations

import collections
import subprocess
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .reference import fetch_ref_base


# Bounded stderr tail kept from ``bcftools query`` for diagnostics on
# non-zero exit. Same convention as ``validate.py``'s ref-check gate.
_BCFTOOLS_STDERR_TAIL_LINES = 200


# ---------------------------------------------------------------------------
# Channel definition
# ---------------------------------------------------------------------------
# Canonical SigProfiler ordering: 6 substitutions, each crossed with 16
# trinucleotide contexts (4 left flanks × 4 right flanks). Channel
# index = substitution_idx * 16 + left_idx * 4 + right_idx. Locking the
# order means the JSON output is comparable across versions and against
# COSMIC's published signature vectors (which use the same convention).

_SUBSTITUTIONS: tuple[str, ...] = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
_FLANKS: tuple[str, ...] = ("A", "C", "G", "T")
_COMPLEMENT: dict[str, str] = {"A": "T", "C": "G", "G": "C", "T": "A"}
_PYRIMIDINES: frozenset[str] = frozenset({"C", "T"})
N_CHANNELS = 96


def channel_index(ref_pyr: str, alt_pyr: str,
                  left: str, right: str) -> int:
    """Return 0..95 channel index for a pyrimidine-context SNV.

    Inputs must already be pyrimidine-normalised; pass purine-context
    bases through :func:`pyrimidine_normalize` first. Raises
    ``ValueError`` on any out-of-alphabet input so a silent mis-binning
    of an ``N`` flank or an unexpected ALT can't poison the spectrum.
    """
    sub = f"{ref_pyr}>{alt_pyr}"
    try:
        sub_idx = _SUBSTITUTIONS.index(sub)
        left_idx = _FLANKS.index(left)
        right_idx = _FLANKS.index(right)
    except ValueError as exc:
        raise ValueError(
            f"channel_index: out-of-alphabet input "
            f"(ref={ref_pyr!r}, alt={alt_pyr!r}, "
            f"left={left!r}, right={right!r}): {exc}"
        ) from exc
    return sub_idx * 16 + left_idx * 4 + right_idx


def channel_label(idx: int) -> str:
    """Render a channel index as ``"<left>[<ref>><alt>]<right>"``.

    Matches the canonical SigProfiler / COSMIC label format so the
    JSON output's keys line up with published signature tables.
    """
    if not 0 <= idx < N_CHANNELS:
        raise ValueError(f"channel index out of range: {idx}")
    sub_idx, ctx_idx = divmod(idx, 16)
    left_idx, right_idx = divmod(ctx_idx, 4)
    return (
        f"{_FLANKS[left_idx]}"
        f"[{_SUBSTITUTIONS[sub_idx]}]"
        f"{_FLANKS[right_idx]}"
    )


def all_channel_labels() -> list[str]:
    """All 96 channel labels in canonical (index) order."""
    return [channel_label(i) for i in range(N_CHANNELS)]


def pyrimidine_normalize(
    ref: str, alt: str, left: str, right: str,
) -> tuple[str, str, str, str] | None:
    """Fold a SNV's REF/ALT/flank tuple to pyrimidine-context form.

    Returns ``(ref_pyr, alt_pyr, left_pyr, right_pyr)`` or ``None`` when
    any input is outside ``{A, C, G, T}`` — N-flanked or otherwise
    non-canonical SNVs are excluded from the spectrum rather than
    silently funneled into a wrong channel.

    Reverse-complement semantics: if ref is a purine (A or G), the
    three flanking bases get complemented AND the flank order
    reverses (the original 3' flank becomes the new 5' flank), so
    e.g. A[G>T]C becomes G[C>A]T after RC.
    """
    if any(b not in _COMPLEMENT for b in (ref, alt, left, right)):
        return None
    if ref in _PYRIMIDINES:
        return ref, alt, left, right
    return (
        _COMPLEMENT[ref],
        _COMPLEMENT[alt],
        _COMPLEMENT[right],  # new left = complemented OLD right
        _COMPLEMENT[left],   # new right = complemented OLD left
    )


# ---------------------------------------------------------------------------
# VCF/BCF iteration
# ---------------------------------------------------------------------------
# We pull (chrom, pos, ref, alt) via ``bcftools query`` rather than a
# full VCF parse — the spectrum only cares about those four fields and
# bcftools handles BCF + tabix-style streaming without us writing a
# parser. Matches the dependency pattern used by ``validate.py``'s
# ``check_ref_against_fasta``.


def _iter_snv_loci(
    vcf_path: Path,
) -> Iterator[tuple[str, int, str, str]]:
    """Yield ``(chrom, pos, ref, alt)`` for every biallelic SNV.

    Filters at the bcftools-query layer with an include expression:
    ``-i 'TYPE="snp" && N_ALT=1'`` keeps only biallelic SNVs. A single
    subprocess (rather than a ``bcftools view | bcftools query`` pipe)
    means one stderr stream + one return code to manage — eliminates
    the deadlock risk of a two-process pipeline where the upstream
    proc blocks on a full stderr buffer while we're consuming
    downstream stdout.

    Indels, SVs, and multi-allelic SNVs are filtered by bcftools, but
    we also defensively skip any row whose REF or ALT isn't a single
    base (e.g. a passed-through deletion that the filter let slip).

    Raises ``RuntimeError`` if ``bcftools query`` exits non-zero
    (e.g. corrupt BCF, missing index, bad input format) with a
    bounded stderr tail in the message — fail-fast rather than
    silently producing an empty spectrum.
    """
    cmd = [
        "bcftools", "query",
        "-i", 'TYPE="snp" && N_ALT=1',
        "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
        str(vcf_path),
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True,
    )
    stderr_tail: collections.deque[str] = collections.deque(
        maxlen=_BCFTOOLS_STDERR_TAIL_LINES,
    )
    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                continue
            chrom, pos_str, ref, alt = parts
            if len(ref) != 1 or len(alt) != 1:
                # Defensive: bcftools should have filtered these out,
                # but a future schema change shouldn't silently
                # contaminate the spectrum.
                continue
            try:
                pos = int(pos_str)
            except ValueError:
                continue
            yield chrom, pos, ref, alt
    finally:
        # Drain stderr after stdout is exhausted (or the iterator was
        # abandoned). Under normal operation bcftools query writes
        # nothing to stderr, so this loop is empty and no buffer-
        # fill deadlock is possible in practice. On error, stderr
        # holds the diagnostic message that became the RuntimeError
        # below.
        assert proc.stderr is not None
        for line in proc.stderr:
            stderr_tail.append(line)
        proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(
            f"bcftools query exited {proc.returncode} on {vcf_path}: "
            f"{''.join(stderr_tail).strip()[-1500:] or '(no stderr)'}"
        )


# ---------------------------------------------------------------------------
# Spectrum computation
# ---------------------------------------------------------------------------


@dataclass
class MutationSpectrum:
    """96-channel mutation-spectrum result for a single VCF/BCF.

    ``counts[i]`` is the raw SNV count in channel ``i`` (see
    :func:`channel_label` for the i → label mapping). ``n_total`` is
    the count of biallelic SNVs read from the source; ``n_excluded`` is
    the count of those that couldn't be binned (N flank, off-chrom-end
    flank, IUPAC ambiguity). The sum of ``counts`` plus ``n_excluded``
    equals ``n_total``.
    """
    counts: list[int] = field(
        default_factory=lambda: [0] * N_CHANNELS,
    )
    n_total: int = 0
    n_excluded: int = 0

    def fractions(self) -> list[float]:
        """``counts`` normalised to sum to 1.0 (zero vector when no
        binnable SNVs were seen)."""
        binned = sum(self.counts)
        if binned == 0:
            return [0.0] * N_CHANNELS
        return [c / binned for c in self.counts]

    def add(self, other: "MutationSpectrum") -> None:
        """Accumulate another spectrum's counts into this one."""
        for i in range(N_CHANNELS):
            self.counts[i] += other.counts[i]
        self.n_total += other.n_total
        self.n_excluded += other.n_excluded


def compute_spectrum(
    vcf_path: Path, fasta: Any,
) -> MutationSpectrum:
    """Walk a VCF/BCF's biallelic SNVs and bin them into 96 channels.

    ``fasta`` is a ``pysam.FastaFile`` (or compatible) — typically the
    same handle the cli already loaded for M12. The trinucleotide
    context at each SNV is ``[fetch_ref_base(pos-1), REF, fetch_ref_base(pos+1)]``.
    Records whose flank lookup returns ``N`` (off the chrom end, IUPAC
    ambiguity, or non-matching chrom name) are counted in
    ``n_excluded`` rather than silently dropped.
    """
    result = MutationSpectrum()
    for chrom, pos, ref, alt in _iter_snv_loci(vcf_path):
        result.n_total += 1
        # Flank lookups use ``fetch_ref_base`` so the chr-prefix
        # convention (e.g. "22" vs "chr22") is handled the same way as
        # M12's REF resolution.
        left = fetch_ref_base(fasta, chrom, pos - 1)
        right = fetch_ref_base(fasta, chrom, pos + 1)
        normalised = pyrimidine_normalize(ref, alt, left, right)
        if normalised is None:
            result.n_excluded += 1
            continue
        ref_pyr, alt_pyr, left_pyr, right_pyr = normalised
        try:
            idx = channel_index(ref_pyr, alt_pyr, left_pyr, right_pyr)
        except ValueError:
            # Should be unreachable after ``pyrimidine_normalize``
            # accepted the inputs, but defensive: any mis-classification
            # lands in the exclusion bucket rather than poisoning a
            # specific channel.
            result.n_excluded += 1
            continue
        result.counts[idx] += 1
    return result


def to_jsonable(spectrum: MutationSpectrum) -> dict:
    """Render a ``MutationSpectrum`` as a JSON-serialisable dict.

    Output schema (matches what ``validate_batch.py`` writes to
    ``validation/mutation_spectrum.json``)::

        {
            "n_total": int,            # biallelic SNVs read
            "n_excluded": int,         # excluded (N flanks etc.)
            "n_binned": int,           # n_total - n_excluded
            "channels": [
                {"label": "A[C>A]A", "count": int, "fraction": float},
                …  # 96 entries in canonical channel order
            ],
        }
    """
    fractions = spectrum.fractions()
    labels = all_channel_labels()
    n_binned = sum(spectrum.counts)
    return {
        "n_total": spectrum.n_total,
        "n_excluded": spectrum.n_excluded,
        "n_binned": n_binned,
        "channels": [
            {
                "label": labels[i],
                "count": spectrum.counts[i],
                "fraction": fractions[i],
            }
            for i in range(N_CHANNELS)
        ],
    }
