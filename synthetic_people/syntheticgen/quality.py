"""Simulated sequencing quality metrics: DP, AD, GQ.

The model is deliberately simple but distributionally realistic and
multi-allelic aware:

* **DP** ~ Poisson(λ) per site; λ is a per-sample baseline with modest
  jitter so different samples in a cohort show different coverage
  profiles.
* **AD** is `Number=R` — one depth per allele (REF + every ALT). `0|0`
  lands all reads on REF; `k|k` lands all on alt k; a het `a|b` splits
  reads between allele a and b via a binomial. Hets involving REF carry
  a small reference bias (P(ref read) ≈ 0.525) matching empirical
  Illumina + BWA-MEM behaviour; `1|2`-style hets split 50/50.
  `sum(AD) == DP` always holds, which `bcftools stats` cross-checks.
* **GQ** is a Phred-like score from (AD, called GT): high when AD
  strongly supports the call, low when AD contradicts or is sparse.
  Floored at 0, capped at 99 per VCF convention.

Everything takes an explicit `random.Random` instance so seeded runs
stay deterministic.
"""

from __future__ import annotations

import math
import random


# Defaults — tunable later via CLI if needed.
DEFAULT_DP_MEAN = 30.0
# Per-sample λ jitter: spread across samples so a 500-person cohort shows
# ~25–35 mean-depth variation, not identical λ for everyone.
DEFAULT_DP_SAMPLE_JITTER_SD = 3.0
# Reference bias on heterozygote sites — ~4–6% ref-leaning is typical for
# Illumina + BWA-MEM short reads. This is the expected ALT-read fraction at
# a het (slightly below 0.5 → the REF allele is over-represented).
HET_ALT_FRAC = 0.475


def sample_lambda(base: float, jitter_sd: float,
                  rng: random.Random) -> float:
    """Per-sample mean depth: base depth + Gaussian jitter, clamped low.

    Clamped at 5 so a pathological jitter draw doesn't produce DP=0 runs.
    """
    lam = base + rng.gauss(0, jitter_sd)
    return max(5.0, lam)


def poisson(lam: float, rng: random.Random) -> int:
    """Knuth's algorithm — fine for the modest λ (5–50) we use here.

    Floors at 0. For larger λ a normal approximation would be faster,
    but we stay well inside the regime where Knuth is exact and fast.
    """
    L = math.exp(-lam)
    k = 0
    p = 1.0
    while p > L:
        k += 1
        p *= rng.random()
    return k - 1


def _binomial(n: int, p: float, rng: random.Random) -> int:
    """Binomial sampling.

    Uses Python 3.12+ `random.binomialvariate` when available, otherwise
    falls back to n coin flips (fine at our DP range).
    """
    binomfn = getattr(rng, "binomialvariate", None)
    if binomfn is not None:
        return binomfn(n, p)
    return sum(1 for _ in range(n) if rng.random() < p)


def _parse_gt_alleles(gt: str) -> tuple[int, int]:
    """Parse a GT into an ``(a, b)`` pair of allele indices.

    Handles three input shapes:

    - **Diploid phased** (``"0|2"``): returns the two allele indices.
    - **Haploid** (``"1"``) — M13.3 chrX non-PAR in males, chrY non-
      PAR in males, MT in everyone. Returned as ``(k, k)`` so the
      AD-draw + GQ-recompute treat the single haplotype as
      homozygous-equivalent at the quality-model level. This keeps
      AD consistent with the emitted GT (haploid alt → AD looks
      hom-alt, not mistakenly hom-ref).
    - **Non-numeric / malformed** (``"."`` etc.): degrades to
      ``(0, 0)``. The quality model only needs dosage info so the
      hom-ref fallback is safe.
    """
    parts = gt.split("|")
    if len(parts) == 1:
        # Haploid: interpret single-allele "k" as (k, k) for quality
        # purposes — same dosage shape as a homozygous diploid.
        try:
            k = int(parts[0])
            return k, k
        except ValueError:
            return 0, 0
    if len(parts) != 2:
        return 0, 0
    try:
        return int(parts[0]), int(parts[1])
    except ValueError:
        return 0, 0


def ad_from_gt(gt: str, n_alleles: int, dp: int,
               rng: random.Random) -> tuple[int, ...]:
    """Split DP into per-allele depths consistent with the genotype.

    Returns a tuple of length `n_alleles` (REF at index 0, ALTs after).
    `sum(...) == dp` always holds.

    Reference bias: when the GT is a het that includes REF (e.g. `0|1`,
    `0|2`), the ref slot gets a small read-share boost over a perfect
    50/50, matching empirical short-read aligner bias. Hets between two
    non-ref alleles (e.g. `1|2`) split 50/50 with no bias.
    """
    if dp == 0 or n_alleles <= 0:
        return (0,) * max(n_alleles, 0)
    a, b = _parse_gt_alleles(gt)
    # Clamp bad indices so we never index out of bounds.
    a = 0 if a < 0 or a >= n_alleles else a
    b = 0 if b < 0 or b >= n_alleles else b

    if a == b:
        counts = [0] * n_alleles
        counts[a] = dp
        return tuple(counts)

    # Heterozygote: one binomial between the two alleles present on the
    # haplotypes. If REF is one of them, apply ref-bias; otherwise 50/50.
    if a == 0:
        p_a = 1.0 - HET_ALT_FRAC  # ref (at `a`) reads a bit more
    elif b == 0:
        p_a = HET_ALT_FRAC  # ref is at `b`, so allele `a` gets the smaller share
    else:
        p_a = 0.5
    a_reads = _binomial(dp, p_a, rng)
    counts = [0] * n_alleles
    counts[a] = a_reads
    counts[b] = dp - a_reads
    return tuple(counts)


def gq_from_ad(gt: str, ad: tuple[int, ...]) -> int:
    """Phred-like GQ from (AD per allele) + called GT.

    Support = fraction of reads consistent with the called genotype:
      * `k|k`      → `AD[k] / DP`
      * `a|b` het  → how close the split between AD[a] and AD[b] is to
                     the expected 50/50, weighted by (AD[a] + AD[b]) / DP
                     so reads supporting *other* alleles count against it.

    GQ = `-10 log10(1 - support)`, depth-capped, clipped to [0, 99].
    """
    dp = sum(ad)
    if dp == 0:
        return 0
    a, b = _parse_gt_alleles(gt)
    n = len(ad)
    if a < 0 or a >= n or b < 0 or b >= n:
        return 0

    if a == b:
        support = ad[a] / dp
    else:
        used = ad[a] + ad[b]
        if used == 0:
            support = 0.0
        else:
            a_frac = ad[a] / used
            support = (1.0 - 2.0 * abs(a_frac - 0.5)) * (used / dp)

    gq = -10.0 * math.log10(max(1.0 - support, 1e-10))
    gq = min(gq, 10.0 * math.log10(max(dp, 1)) * 6.0)
    return max(0, min(99, round(gq)))


def draw_site_quality(gt: str, n_alleles: int, lam: float,
                      rng: random.Random) -> tuple[int, tuple[int, ...], int]:
    """Return `(DP, AD, GQ)` for one record with `n_alleles` total alleles.

    AD is always length `n_alleles` (REF + every ALT) per VCF `Number=R`.
    """
    dp = poisson(lam, rng)
    if dp == 0:
        return 0, (0,) * n_alleles, 0
    ad = ad_from_gt(gt, n_alleles, dp, rng)
    gq = gq_from_ad(gt, ad)
    return dp, ad, gq
