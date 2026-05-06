"""Tests for cohort-level site generation."""

from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.cohort import (  # noqa: E402
    assign_haplotypes,
    draw_cohort_background,
    person_records_from_cohort,
)
from syntheticgen.cohort_sites import (  # noqa: E402
    carriers_from_dense_gts,
    dense_gts_from_carriers,
)


def _site_with_dense_gts(gts, n_people=None, **fields):
    """Test helper — build a Phase-5c-shaped site dict from a
    human-readable list of GT strings.

    Tests stay readable (gts as ``["0|1", "1|1", ...]``) while the
    storage on the site itself matches what production code emits
    (sparse carriers).
    """
    if n_people is None:
        n_people = len(gts)
    site = {
        "n_haplotypes": 2 * n_people,
        "carriers": carriers_from_dense_gts(gts),
    }
    site.update(fields)
    return site


def _pool(n: int = 50, seed: int = 99) -> list:
    rng = random.Random(seed)
    chroms = ("1", "2", "22")
    pool = []
    for i in range(n):
        pool.append({
            "chrom": rng.choice(chroms),
            "pos": 100_000 + i * 37,
            "id": ".",
            "ref": rng.choice(("A", "C", "G", "T")),
            "alts": [rng.choice(("A", "C", "G", "T"))],
            "afs": [rng.random() * 0.5],
        })
    return pool


class TestAssignHaplotypes(unittest.TestCase):
    def test_exact_counts_preserved(self):
        rng = random.Random(0)
        slots = assign_haplotypes(20, [3, 1], rng)
        self.assertEqual(len(slots), 20)
        self.assertEqual(slots.count(1), 3)
        self.assertEqual(slots.count(2), 1)
        self.assertEqual(slots.count(0), 16)

    def test_biallelic(self):
        rng = random.Random(1)
        slots = assign_haplotypes(10, [4], rng)
        self.assertEqual(slots.count(1), 4)
        self.assertEqual(slots.count(0), 6)

    def test_rejects_overflow(self):
        rng = random.Random(0)
        with self.assertRaises(ValueError):
            assign_haplotypes(4, [3, 2], rng)

    def test_shuffles_slots(self):
        """Slot assignment should be randomised, not packed at the front."""
        rng = random.Random(2)
        observed_positions = set()
        for _ in range(50):
            slots = assign_haplotypes(30, [1], rng)
            observed_positions.add(slots.index(1))
        # With 50 independent draws across 30 slots, we expect many
        # distinct positions to show up.
        self.assertGreater(len(observed_positions), 10)


class TestDrawCohortBackground(unittest.TestCase):
    def test_returns_requested_number_of_sites(self):
        rng = random.Random(7)
        pool = _pool(100)
        sites = draw_cohort_background(pool, n_people=20, n_sites=40,
                                       alpha=2.0, rng=rng)
        self.assertEqual(len(sites), 40)

    def test_shared_coordinates_across_cohort(self):
        """Every person sees the same (chrom, pos, ref, alts) at each site."""
        rng = random.Random(3)
        pool = _pool(60)
        sites = draw_cohort_background(pool, n_people=10, n_sites=25,
                                       alpha=2.0, rng=rng)
        for site in sites:
            # Phase 5c: GTs live as sparse carriers; n_haplotypes
            # tracks the cohort cardinality.
            self.assertEqual(site["n_haplotypes"], 2 * 10)
            self.assertEqual(len(site["acs"]), len(site["alts"]))
            self.assertEqual(len(site["afs"]), len(site["alts"]))

    def test_realized_ac_matches_drawn_ac(self):
        """Slot assignment is exact — summing per-person dosages hits AC."""
        rng = random.Random(5)
        pool = _pool(40)
        n_people = 15
        sites = draw_cohort_background(pool, n_people=n_people, n_sites=20,
                                       alpha=2.0, rng=rng)
        for site in sites:
            for alt_idx, expected in enumerate(site["acs"], start=1):
                # Sparse carriers: the count of (_, allele) tuples
                # whose allele matches alt_idx is the realised AC.
                realised = sum(
                    1 for _, allele in site["carriers"]
                    if allele == alt_idx
                )
                self.assertEqual(realised, expected,
                                 f"AC mismatch at {site['chrom']}:{site['pos']}")

    def test_every_site_is_variable(self):
        """No site should be fixed for REF (total alt count ≥ 1)."""
        rng = random.Random(8)
        pool = _pool(80)
        sites = draw_cohort_background(pool, n_people=12, n_sites=40,
                                       alpha=2.0, rng=rng)
        for site in sites:
            self.assertGreaterEqual(sum(site["acs"]), 1)

    def test_empty_pool_returns_nothing(self):
        rng = random.Random(0)
        self.assertEqual(draw_cohort_background([], 10, 5, 2.0, rng), [])

    def test_reproducible_under_seed(self):
        pool = _pool(60)
        a = draw_cohort_background(pool, 8, 15, 2.0, random.Random(42))
        b = draw_cohort_background(pool, 8, 15, 2.0, random.Random(42))
        self.assertEqual(len(a), len(b))
        for s1, s2 in zip(a, b):
            self.assertEqual(s1["chrom"], s2["chrom"])
            self.assertEqual(s1["pos"], s2["pos"])
            self.assertEqual(s1["acs"], s2["acs"])
            # Phase 5c: compare carriers (sparse) instead of dense GTs.
            self.assertEqual(s1["carriers"], s2["carriers"])


class TestPersonRecordsFromCohort(unittest.TestCase):
    def test_drops_hom_ref(self):
        sites = [
            _site_with_dense_gts(
                ["0|1", "0|0", "1|1"], chrom="1", pos=100, id=".",
                ref="A", alts=["G"], afs=[0.2], acs=[1]),
            _site_with_dense_gts(
                ["0|0", "0|0", "0|1"], chrom="1", pos=200, id=".",
                ref="C", alts=["T"], afs=[0.1], acs=[1]),
        ]
        # Person 0: one het, one hom-ref → 1 record.
        self.assertEqual(len(person_records_from_cohort(sites, 0)), 1)
        # Person 1: both hom-ref → 0 records.
        self.assertEqual(len(person_records_from_cohort(sites, 1)), 0)
        # Person 2: hom-alt + het → 2 records.
        self.assertEqual(len(person_records_from_cohort(sites, 2)), 2)

    def test_carries_coords_and_gt(self):
        sites = [
            _site_with_dense_gts(
                ["0|1"], chrom="22", pos=42, id="rs1",
                ref="A", alts=["G"], afs=[0.3], acs=[2]),
        ]
        recs = person_records_from_cohort(sites, 0)
        self.assertEqual(recs[0]["chrom"], "22")
        self.assertEqual(recs[0]["pos"], 42)
        self.assertEqual(recs[0]["id"], "rs1")
        self.assertEqual(recs[0]["gt"], "0|1")


if __name__ == "__main__":
    unittest.main(verbosity=2)
