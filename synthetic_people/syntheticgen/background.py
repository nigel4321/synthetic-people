"""1000G background coordinate pool + genotype helpers.

Reservoir-samples common variants from local 1000 Genomes VCFs to build
a coordinate pool (chrom/pos/ref/alts) the cohort generator can pick
from. The source AFs are kept for back-compat but the M4+ cohort path
redraws fresh SFS-based frequencies per site and ignores them — only
coordinates and allele strings carry through.

The `phased_gt_from_af*` helpers draw phased diploid genotypes under
Hardy-Weinberg and are retained for the highlighted-variant path and
unit tests.
"""

from __future__ import annotations

import glob
import os
import random
import subprocess
import sys


def load_background_pool(globs: list[str], af_min: float,
                         per_source_limit: int,
                         rng: random.Random) -> list[dict]:
    """Reservoir-sample common variants (AF >= af_min) from each source VCF."""
    sources: list[str] = []
    for g in globs:
        sources.extend(sorted(glob.glob(g)))
    # Deduplicate in case of overlapping globs.
    seen = set()
    sources = [s for s in sources if not (s in seen or seen.add(s))]
    if not sources:
        return []

    pool: list[dict] = []
    for src in sources:
        print(f"  sampling background from {os.path.basename(src)}",
              file=sys.stderr)
        # Pull AF with `-i MAX(INFO/AF)>=af_min` so multi-allelic sites
        # qualify whenever any alt is common enough.
        cmd = [
            "bcftools", "query",
            "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n",
            "-i", f"MAX(INFO/AF)>={af_min}",
            src,
        ]
        reservoir: list[dict] = []
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True,
                              stderr=subprocess.DEVNULL) as proc:
            assert proc.stdout is not None
            i = 0
            for line in proc.stdout:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                chrom, pos, ref, alt, af_str = parts[:5]
                # Reject symbolic ALTs and length-capped records. Per-allele
                # length cap so a multi-allelic with one long alt doesn't
                # blow up record size.
                alts = alt.split(",")
                if any(a.startswith("<") or len(a) > 50 for a in alts):
                    continue
                if len(ref) > 50:
                    continue
                try:
                    afs = [float(x) for x in af_str.split(",")]
                except ValueError:
                    continue
                if len(afs) != len(alts):
                    continue
                entry = {"chrom": chrom, "pos": int(pos), "id": ".",
                         "ref": ref, "alts": alts, "afs": afs}
                if i < per_source_limit:
                    reservoir.append(entry)
                else:
                    j = rng.randint(0, i)
                    if j < per_source_limit:
                        reservoir[j] = entry
                i += 1
        pool.extend(reservoir)
    return pool


def phased_gt_from_af(af: float, rng: random.Random) -> str:
    """Phased biallelic diploid genotype under Hardy-Weinberg."""
    a = "1" if rng.random() < af else "0"
    b = "1" if rng.random() < af else "0"
    return f"{a}|{b}"


def phased_gt_from_afs(afs: list, rng: random.Random) -> str:
    """Phased diploid genotype for a multi-allelic site.

    `afs` is the list of alt-allele frequencies at the site. Each
    haplotype draws from the categorical distribution {REF, alt_1, ...}
    with P(REF) = 1 - sum(afs). Returns `"a|b"` where `a` and `b` are
    integer allele indices (0 = REF).
    """
    total_alt = sum(afs)
    if total_alt <= 0:
        return "0|0"
    # Cumulative thresholds across alleles 1..k (REF fills the remainder).
    cum = []
    running = 0.0
    for af in afs:
        running += af
        cum.append(running)

    def _draw_allele() -> int:
        r = rng.random()
        for idx, c in enumerate(cum, start=1):
            if r < c:
                return idx
        return 0  # REF

    a = _draw_allele()
    b = _draw_allele()
    return f"{a}|{b}"


def alt_dosages(gt: str, n_alts: int) -> list:
    """Per-alt dosage list: entry k is count of allele (k+1) on the two haplotypes."""
    dosages = [0] * n_alts
    for token in gt.split("|"):
        try:
            idx = int(token)
        except ValueError:
            continue
        if 1 <= idx <= n_alts:
            dosages[idx - 1] += 1
    return dosages


def alt_dosage(gt: str) -> int:
    """Total alt dosage across all alleles — for legacy biallelic call sites."""
    total = 0
    for token in gt.split("|"):
        try:
            idx = int(token)
        except ValueError:
            continue
        if idx >= 1:
            total += 1
    return total


def random_sample_id(rng: random.Random) -> str:
    """HG/NA-prefixed 5-digit ID, mirroring 1000G naming conventions.

    .. warning::
       This single-draw function is **only safe for very small
       cohorts** (n ≲ 100) — the prefix × 5-digit number gives just
       180 000 unique IDs, so by birthday-paradox math collisions
       become routine past n ≈ 600. Cohort writers must use unique
       sample names (bcftools rejects multi-sample BCFs with
       duplicate columns), so anything constructing a list of N IDs
       should call :func:`draw_sample_ids` instead.
    """
    prefix = rng.choice(("HG", "NA"))
    return f"{prefix}{rng.randint(10000, 99999)}"


# ID-pool dimensions for :func:`draw_sample_ids`. The number range
# was widened from the legacy 5-digit ``random_sample_id`` to a
# 6-digit + 7-digit superset so the cohort write succeeds at any
# realistic cohort size:
#
#   pool size = 2 prefixes × (1e7 - 1e5) = 19.8M unique IDs
#
# Plenty for n ≤ 1M (the Phase 5 stretch target). Above that, the
# function raises rather than silently risking a duplicate.
_SAMPLE_ID_PREFIXES = ("HG", "NA")
_SAMPLE_ID_NUMBER_LO = 100_000
_SAMPLE_ID_NUMBER_HI = 10_000_000
_SAMPLE_ID_POOL_SIZE = (
    len(_SAMPLE_ID_PREFIXES)
    * (_SAMPLE_ID_NUMBER_HI - _SAMPLE_ID_NUMBER_LO)
)


def draw_sample_ids(n_people: int, rng: random.Random) -> list:
    """Draw ``n_people`` distinct ``HG``/``NA``-prefixed sample IDs.

    Uniqueness is guaranteed by construction: a single
    ``rng.sample(range(pool_size), n_people)`` shuffle picks integer
    keys without replacement; each key decodes deterministically to a
    ``(prefix, number)`` pair. So ``draw_sample_ids(3000, rng)`` cannot
    produce a duplicate even after thousands of draws — the way the
    old per-call ``random_sample_id`` could and did. bcftools rejects
    multi-sample BCFs with duplicate sample columns
    (``[E::bcf_hdr_add_sample_len] Duplicated sample name``), so the
    cohort writer needs the uniqueness contract this function offers.

    Determinism: the rng consumption is one ``rng.sample`` call, so a
    re-run with the same ``--seed`` and the same ``n_people`` produces
    the same list of IDs in the same order. Output at a given seed
    differs from runs against the previous per-call
    ``random_sample_id`` because the consumption pattern changed.

    Raises ``ValueError`` if ``n_people`` exceeds the available pool
    (1.98 × 10⁷). At that scale we'd need a wider pool or a different
    naming convention — neither is needed today.
    """
    if n_people <= 0:
        return []
    if n_people > _SAMPLE_ID_POOL_SIZE:
        raise ValueError(
            f"n_people={n_people} exceeds the unique sample-ID pool "
            f"({_SAMPLE_ID_POOL_SIZE}); widen the pool in "
            f"background.py or switch to a longer naming convention"
        )
    span = _SAMPLE_ID_NUMBER_HI - _SAMPLE_ID_NUMBER_LO
    keys = rng.sample(range(_SAMPLE_ID_POOL_SIZE), n_people)
    out = []
    for k in keys:
        prefix_idx, num_in_prefix = divmod(k, span)
        prefix = _SAMPLE_ID_PREFIXES[prefix_idx]
        number = _SAMPLE_ID_NUMBER_LO + num_in_prefix
        out.append(f"{prefix}{number}")
    return out
