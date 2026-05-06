"""Sparse genotype storage for cohort sites (Phase 5c).

Cohort sites carry their per-haplotype allele assignments as a list
of ``(haplotype_index, allele_index)`` tuples — the *carriers* — for
non-zero entries only. The dense list-of-GT-strings representation
that earlier phases used is recovered on demand at write time, but
never stored on the cohort_sites list itself.

Why: at ``n=100 000 × n_sites=100 000`` the dense
``gts: list[str]`` representation costs roughly 1 TB per chromosome
(``n × n_sites × ~100 B``). Sparse storage costs ``O(alt
observations)``, which the SFS makes ``O(n_sites × log(n))`` —
flat enough that n=100k fits in a few hundred MB and n=1M in a few
GB. See PERFORMANCE_PLAN §"Phase 5c" for the full memory model.

Cohort site dict shape after 5c:

::

    {
        "chrom": "22",
        "pos": 12345,
        "id": ".",
        "ref": "A",
        "alts": ["G"],
        "afs":  [0.25],
        "acs":  [10],
        "n_haplotypes": 40,            # 2 × n_people; needed for expansion
        "carriers": [(2, 1), (5, 1)],  # sparse — only non-zero entries
    }

The helper module's job is to keep tests, the writer, and the
in-memory legacy/admixture per-person fan-out comfortable with this
shape. Tests can keep constructing fixtures from dense GT lists via
:func:`carriers_from_dense_gts`; the writer expands carriers to
dense at write time via :func:`dense_gts_from_carriers`; the
in-memory per-person path uses :func:`gt_for_person` to scan one
person's GT out of a site's carriers.
"""

from __future__ import annotations


def carriers_from_dense_gts(gts: list) -> list:
    """Convert dense GT strings ``["0|0", "0|1", "1|1"]`` to sparse
    carriers ``[(hap_idx, allele_idx), ...]`` for non-zero entries.

    Phased GT strings ("0|1") and unphased ("0/1") both work — the
    separator is detected. The first allele in each string lands at
    haplotype index ``2 * person_idx``, the second at
    ``2 * person_idx + 1``.

    Used by tests that build site fixtures from human-readable dense
    GTs and by callers migrating from the old ``site["gts"]`` shape.
    """
    carriers: list = []
    for person_idx, gt in enumerate(gts):
        if "|" in gt:
            a1, _, a2 = gt.partition("|")
        elif "/" in gt:
            a1, _, a2 = gt.partition("/")
        else:
            # Single-allele or unrecognised; treat as hom-ref.
            continue
        try:
            a1_int = int(a1)
            a2_int = int(a2)
        except ValueError:
            # "./." or other missing — record as hom-ref (consistent
            # with person_records_from_cohort's drop-only-all-zero
            # semantics, which keeps "./." as a non-event).
            continue
        if a1_int != 0:
            carriers.append((2 * person_idx, a1_int))
        if a2_int != 0:
            carriers.append((2 * person_idx + 1, a2_int))
    return carriers


def dense_gts_from_carriers(carriers, n_people: int) -> list:
    """Expand sparse carriers to a length-``n_people`` list of
    ``"a|b"`` GT strings.

    Used at BCF-write time when the multi-sample VCF row needs the
    dense per-sample GT block. ``O(n_people + carriers)``: walk the
    carriers once to populate the non-zero slots, then format every
    person's pair to a string.

    Multi-allelic alleles (``allele_idx > 1``) format correctly via
    f-string — supports the legacy 1000G path's multi-alt records.
    """
    pairs = [[0, 0] for _ in range(n_people)]
    for hap_idx, allele_idx in carriers:
        person_idx, pair_idx = divmod(hap_idx, 2)
        pairs[person_idx][pair_idx] = allele_idx
    return [f"{p[0]}|{p[1]}" for p in pairs]


def gt_for_person(carriers, person_idx: int) -> str:
    """Return the GT string for one person without expanding the rest.

    Used by the in-memory legacy/admixture per-person fan-out
    (``cohort.person_records_from_cohort``). Per-site cost is
    ``O(carriers in this site)`` — typically a small constant since
    most sites carry only a handful of alt observations.
    """
    a1 = a2 = 0
    target_lo = 2 * person_idx
    target_hi = target_lo + 1
    for hap_idx, allele_idx in carriers:
        if hap_idx == target_lo:
            a1 = allele_idx
        elif hap_idx == target_hi:
            a2 = allele_idx
    return f"{a1}|{a2}"
