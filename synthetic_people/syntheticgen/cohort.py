"""Cohort-level site generation — shared coordinates, per-person genotypes.

From M4 onwards the generator simulates an N-sample cohort as a single
event: one pass selects the variable sites and their allele frequencies,
and each site is then "populated" by assigning alt alleles into specific
haplotype slots. Diploid GTs per person fall out of pairing consecutive
slot assignments. This replaces the M1–M3 per-person independent sampler.

Exact slot assignment (sampling without replacement from 2n haplotypes)
preserves the drawn minor allele count exactly, so the realised cohort
SFS matches the SFS we drew from — no smoothing from independent HWE
resampling at every site. The pairing step inherits a mild
hypergeometric-vs-HWE correction, which is if anything more realistic
for a finite cohort than strict HWE draws.

Phase 5c: per-site genotypes are stored as sparse ``carriers`` —
``[(haplotype_idx, allele_idx), ...]`` for non-zero entries only —
instead of the dense ``gts: list[str]`` shape earlier phases used.
RAM scales with alt observations rather than ``n × n_sites``, so
cohort sizes that previously OOM-killed (≥ 3 000 with full 22-chrom
× 70 Mb input) finish on a workstation. See ``cohort_sites.py`` for
the helper functions that round-trip between dense GT strings and
sparse carriers — used by tests, the BCF writer, and the in-memory
per-person fan-out.
"""

from __future__ import annotations

import random

from .sfs import DEFAULT_SFS_ALPHA, draw_allele_counts


def assign_haplotypes(n_haplotypes: int, allele_counts: list,
                      rng: random.Random) -> list:
    """Assign alt alleles to haplotype slots.

    Returns a length-n_haplotypes list where each entry is an allele
    index (0 = REF, 1..k = alt_1..alt_k). The assignment is uniformly
    random across slots without replacement, so the realised count of
    each allele equals the input `allele_counts` exactly.
    """
    total_alt = sum(allele_counts)
    if total_alt > n_haplotypes:
        raise ValueError(
            f"total alt count {total_alt} exceeds haplotype slots "
            f"{n_haplotypes}"
        )
    slots = [0] * n_haplotypes
    positions = list(range(n_haplotypes))
    rng.shuffle(positions)
    cursor = 0
    for alt_index, count in enumerate(allele_counts, start=1):
        for _ in range(count):
            slots[positions[cursor]] = alt_index
            cursor += 1
    return slots


def _carriers_from_slots(slots: list) -> list:
    """Convert a dense slot array to sparse carriers (Phase 5c).

    ``slots[i]`` is the allele index assigned to haplotype slot ``i``;
    we keep only the non-zero entries. The slot array can be released
    by the caller once carriers are extracted, so RAM is bounded by
    alt observations rather than ``n_haplotypes``.
    """
    return [
        (idx, allele) for idx, allele in enumerate(slots) if allele > 0
    ]


def draw_cohort_background(pool: list, n_people: int, n_sites: int,
                           alpha: float,
                           rng: random.Random) -> list:
    """Generate the cohort's shared background sites.

    `pool` provides candidate coordinates (chrom, pos, ref, alts, ...).
    For each chosen site, allele counts are redrawn from the power-law
    SFS (the 1000G source AFs are ignored — only its genomic coordinates
    and allele strings carry through). Alt alleles are then placed into
    specific haplotype slots so that each person's genotype at the site
    is consistent with every other person in the cohort.
    """
    if not pool or n_people < 1 or n_sites < 1:
        return []
    n_haplotypes = 2 * n_people
    sample_n = min(n_sites, len(pool))
    chosen = rng.sample(pool, sample_n)
    sites = []
    for entry in chosen:
        n_alts = len(entry["alts"])
        counts = draw_allele_counts(n_haplotypes, n_alts, alpha, rng)
        afs = [c / n_haplotypes for c in counts]
        slots = assign_haplotypes(n_haplotypes, counts, rng)
        sites.append({
            "chrom": entry["chrom"],
            "pos": entry["pos"],
            "id": entry.get("id", "."),
            "ref": entry["ref"],
            "alts": list(entry["alts"]),
            "afs": afs,
            "acs": counts,
            "n_haplotypes": n_haplotypes,
            "carriers": _carriers_from_slots(slots),
        })
    return sites


def person_records_from_cohort(sites: list, person_index: int) -> list:
    """Project cohort sites down to one person's non-hom-ref records.

    Carries the M7 annotation fields (clnsig, clndn, cosmic_id,
    cosmic_gene) through to the per-person record dict so the writer can
    emit the corresponding INFO tags for that sample.

    Phase 5c: scans the per-site sparse carriers list to find this
    person's GT rather than indexing into a dense ``gts`` list.
    Per-site cost is ``O(carriers)`` — a small constant under
    SFS-realistic cohorts where most sites carry a handful of alt
    observations.
    """
    target_lo = 2 * person_index
    target_hi = target_lo + 1
    records = []
    for site in sites:
        a1 = a2 = 0
        for hap_idx, allele_idx in site.get("carriers", ()):
            if hap_idx == target_lo:
                a1 = allele_idx
            elif hap_idx == target_hi:
                a2 = allele_idx
        if a1 == 0 and a2 == 0:
            continue
        gt = f"{a1}|{a2}"
        rec = {
            "chrom": site["chrom"],
            "pos": site["pos"],
            "id": site["id"],
            "ref": site["ref"],
            "alts": site["alts"],
            "afs": site["afs"],
            "gt": gt,
        }
        for key in ("clnsig", "clndn", "cosmic_id", "cosmic_gene"):
            if site.get(key):
                rec[key] = site[key]
        records.append(rec)
    return records
