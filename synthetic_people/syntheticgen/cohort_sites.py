"""Sparse genotype storage for cohort sites (Phase 5c).

Cohort sites carry their per-haplotype allele assignments as a
**packed** ``np.ndarray`` of shape ``(n_carriers, 2)`` and dtype
``np.int32`` — the *carriers* — with one row per non-zero entry.
The dense list-of-GT-strings representation that earlier phases
used is recovered on demand at write time, but never stored on the
cohort_sites list itself.

Why packed: at n=1M a common-AF site (~30%) carries ~600 K rows;
as a Python ``list[tuple[int, int]]`` that was ~80 B/row ≈ 48 MB
per site, dominating parent RSS in the streaming pipeline (see
``PERFORMANCE_BUDGETS.md`` § "Known scaling ceiling" for the
2026-05-12 n=1M incident). As a 2D ``np.int32`` array it's 8 B/row
≈ 4.8 MB per site — ~10x reduction across every carriers-holding
code path (streaming heap, Arrow writer batches, materialised
sites list, per-person fan-out).

Why 2D / keep the allele column: the bcf_writer and Arrow paths
both already support multi-allelic carriers (e.g.
``[(0,1),(2,1),(5,2),(7,2)]`` for a tri-allelic site), and the
M14 mutation-model work is expected to start producing them. The
packed shape preserves that contract — ``for hap_idx, allele_idx
in carriers`` unpacks 2D rows naturally, so most consumers needed
no change when the representation flipped.

Cohort site dict shape:

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
        "carriers": np.array([[2, 1], [5, 1]], dtype=np.int32),
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

import numpy as np


# Empty-carriers sentinel returned by producers when a site has no
# alt observations (theoretically possible for the legacy 1000G
# path; monomorphic sites are filtered upstream of the coalescent
# producers). Shape (0, 2) preserves the 2D contract so consumers
# can iterate without an explicit ``if len(carriers) == 0`` guard.
_EMPTY_CARRIERS = np.zeros((0, 2), dtype=np.int32)


def carriers_from_dense_gts(gts: list) -> np.ndarray:
    """Convert dense GT strings ``["0|0", "0|1", "1|1"]`` to packed
    carriers — ``np.ndarray`` of shape ``(n_carriers, 2)``, dtype
    ``np.int32``. Each row is ``[hap_idx, allele_idx]`` for one
    non-zero entry.

    Phased GT strings ("0|1") and unphased ("0/1") both work — the
    separator is detected. The first allele in each string lands at
    haplotype index ``2 * person_idx``, the second at
    ``2 * person_idx + 1``.

    Used by tests that build site fixtures from human-readable dense
    GTs and by callers migrating from the old ``site["gts"]`` shape.
    """
    rows: list = []
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
            rows.append((2 * person_idx, a1_int))
        if a2_int != 0:
            rows.append((2 * person_idx + 1, a2_int))
    if not rows:
        return _EMPTY_CARRIERS
    return np.array(rows, dtype=np.int32)


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


def dense_gts_from_carriers_slice(carriers,
                                  slice_lo: int,
                                  slice_hi: int) -> list:
    """Expand sparse carriers to dense GT strings for a sample slice
    only — persons in ``[slice_lo, slice_hi)``.

    Phase 5e Phase A: workers in the parallel cohort BCF write each
    handle a contiguous sample slice. Per-site they need to format
    only their slice's GT block, not the whole cohort. This is the
    slice-aware sibling of :func:`dense_gts_from_carriers`.

    Returns a list of length ``slice_hi - slice_lo``. Carriers with
    haplotype indices outside the slice are skipped — single pass,
    ``O(carriers + slice_size)``.

    Slice bounds are bounded-half-open in *person* index space, not
    haplotype index space (mirrors how ``sample_ids[lo:hi]`` slices
    in the caller). The corresponding haplotype range is
    ``[2 * slice_lo, 2 * slice_hi)``.
    """
    if slice_hi <= slice_lo:
        return []
    n = slice_hi - slice_lo
    pairs = [[0, 0] for _ in range(n)]
    hap_lo = 2 * slice_lo
    hap_hi = 2 * slice_hi
    for hap_idx, allele_idx in carriers:
        if hap_lo <= hap_idx < hap_hi:
            local = hap_idx - hap_lo
            person_idx, pair_idx = divmod(local, 2)
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
