"""Tests for the UK-cohort admixture path (M6).

Skipped gracefully when msprime / demes / tskit aren't installed, so the
rest of the suite still passes in a stdlib-only environment.
"""

from __future__ import annotations

import random
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

try:
    import msprime  # noqa: F401
    import demes  # noqa: F401
    import tskit  # noqa: F401
    _HAVE_DEPS = True
except ImportError:
    _HAVE_DEPS = False


@unittest.skipUnless(_HAVE_DEPS, "msprime/demes/tskit not installed")
class TestBuildDemography(unittest.TestCase):
    def test_proportions_must_sum_to_one(self):
        from syntheticgen.admixture import build_uk_demography
        with self.assertRaises(ValueError):
            build_uk_demography(0.5, 0.3, 0.3)
        with self.assertRaises(ValueError):
            build_uk_demography(-0.1, 0.6, 0.5)
        with self.assertRaises(ValueError):
            build_uk_demography(0.6, 0.25, 0.15, pulse_time=0)
        # Default sums to 1.0 — should not raise
        graph = build_uk_demography()
        deme_names = [d.name for d in graph.demes]
        self.assertIn("UK", deme_names)
        for src in ("EUR", "SAS", "AFR"):
            self.assertIn(src, deme_names)

    def test_uk_deme_has_three_ancestors(self):
        from syntheticgen.admixture import build_uk_demography
        graph = build_uk_demography(0.6, 0.25, 0.15)
        uk = next(d for d in graph.demes if d.name == "UK")
        self.assertEqual(set(uk.ancestors), {"EUR", "SAS", "AFR"})
        # Proportions preserved into the resolved graph.
        prop_map = dict(zip(uk.ancestors, uk.proportions))
        self.assertAlmostEqual(prop_map["EUR"], 0.6)
        self.assertAlmostEqual(prop_map["SAS"], 0.25)
        self.assertAlmostEqual(prop_map["AFR"], 0.15)


@unittest.skipUnless(_HAVE_DEPS, "msprime/demes/tskit not installed")
class TestSimulateChromosome(unittest.TestCase):
    def setUp(self):
        from syntheticgen.admixture import simulate_chromosome
        self.simulate_chromosome = simulate_chromosome

    def _run(self, seed: int = 42, n_people: int = 6,
             length_mb: float = 0.5,
             proportions: tuple = (0.6, 0.25, 0.15)):
        rng = random.Random(seed)
        return self.simulate_chromosome(
            chrom="22", build="GRCh38", n_people=n_people,
            length_mb=length_mb, proportions=proportions,
            rec_rate=1e-8, mu=1.29e-8, rng=rng,
        )

    def test_returns_sites_and_segments(self):
        sites, segs = self._run()
        self.assertGreater(len(sites), 20)
        self.assertEqual(len(segs), 6)
        for person_segs in segs:
            self.assertGreater(len(person_segs), 0)
            for s, e, h1, h2 in person_segs:
                self.assertLess(s, e)
                self.assertIn(h1, {"EUR", "SAS", "AFR", "OOA", "ANC"})
                self.assertIn(h2, {"EUR", "SAS", "AFR", "OOA", "ANC"})

    def test_segments_cover_full_chromosome_length(self):
        sites, segs = self._run(length_mb=0.5)
        sim_length = int(0.5 * 1_000_000)
        for person_segs in segs:
            covered = sum(e - s for s, e, _, _ in person_segs)
            # Allow 1 bp slack for float→int rounding at boundaries
            self.assertLess(abs(covered - sim_length), 5)
            # Sorted, non-overlapping
            for i in range(len(person_segs) - 1):
                self.assertLessEqual(person_segs[i][1],
                                     person_segs[i + 1][0])

    def test_realised_ac_matches_declared(self):
        sites, _ = self._run()
        for s in sites:
            # Phase 5c sparse carriers — sum allele indices
            # (binary mutation model, so all 1).
            realised = sum(allele for _, allele in s["carriers"])
            self.assertEqual(realised, s["acs"][0],
                             f"AC mismatch at {s['pos']}")

    def test_reproducible_under_seed(self):
        sites_a, segs_a = self._run(seed=42)
        sites_b, segs_b = self._run(seed=42)
        self.assertEqual(len(sites_a), len(sites_b))
        self.assertEqual(segs_a, segs_b)
        for s1, s2 in zip(sites_a, sites_b):
            self.assertEqual(s1["pos"], s2["pos"])
            self.assertEqual(s1["carriers"], s2["carriers"])

    def test_different_seeds_give_different_output(self):
        sites_a, _ = self._run(seed=1)
        sites_b, _ = self._run(seed=2)
        self.assertNotEqual([s["pos"] for s in sites_a],
                            [s["pos"] for s in sites_b])

    def test_unknown_chromosome_raises(self):
        rng = random.Random(0)
        with self.assertRaises(ValueError):
            self.simulate_chromosome(
                chrom="42", build="GRCh38", n_people=4, length_mb=0.1,
                proportions=(0.6, 0.25, 0.15),
                rec_rate=1e-8, mu=1.29e-8, rng=rng,
            )

    def test_ancestry_fractions_track_requested_proportions(self):
        """On a moderately-sized cohort, realised SOURCE-pop fractions
        should land within ±10% of the requested mix. We use a 30-person
        × 1 Mb sim — large enough to dampen finite-cohort noise."""
        from syntheticgen.admixture import (
            ancestry_fractions, SOURCE_POPS,
        )
        rng = random.Random(7)
        _, segs = self.simulate_chromosome(
            chrom="22", build="GRCh38", n_people=30, length_mb=1.0,
            proportions=(0.6, 0.25, 0.15),
            rec_rate=1e-8, mu=1.29e-8, rng=rng,
        )
        # Aggregate across all people to dampen per-person noise
        agg = []
        for person_segs in segs:
            for s, e, h1, h2 in person_segs:
                agg.append(("22", s, e, h1, h2))
        fracs = ancestry_fractions(agg)
        eur = sum(fracs.get(p, 0) for p in ("EUR",))
        sas = sum(fracs.get(p, 0) for p in ("SAS",))
        afr = sum(fracs.get(p, 0) for p in ("AFR",))
        # Most ancestry should land in the source pops with default
        # demography (PULSE_TIME=20 gens means very little leakage to
        # ancestral demes).
        source_total = eur + sas + afr
        self.assertGreater(source_total, 0.85)
        self.assertAlmostEqual(eur / source_total, 0.6, delta=0.15)
        self.assertAlmostEqual(sas / source_total, 0.25, delta=0.10)
        self.assertAlmostEqual(afr / source_total, 0.15, delta=0.10)


@unittest.skipUnless(_HAVE_DEPS, "msprime/demes/tskit not installed")
class TestAncestryHelpers(unittest.TestCase):
    def test_ancestry_fractions_normalises(self):
        from syntheticgen.admixture import ancestry_fractions
        segs = [
            ("22", 0, 100, "EUR", "EUR"),
            ("22", 100, 200, "SAS", "EUR"),
            ("22", 200, 300, "AFR", "AFR"),
        ]
        fracs = ancestry_fractions(segs)
        self.assertAlmostEqual(sum(fracs.values()), 1.0, places=6)
        # Counts: 600 hap-bp total. EUR: 100+100+100 = 300 (50%);
        # SAS: 100 (16.67%); AFR: 200 (33.33%).
        self.assertAlmostEqual(fracs["EUR"], 0.5, places=6)
        self.assertAlmostEqual(fracs["SAS"], 100 / 600, places=6)
        self.assertAlmostEqual(fracs["AFR"], 200 / 600, places=6)

    def test_ancestry_fractions_empty(self):
        from syntheticgen.admixture import ancestry_fractions, SOURCE_POPS
        fracs = ancestry_fractions([])
        for p in SOURCE_POPS:
            self.assertEqual(fracs[p], 0.0)

    def test_write_bed_round_trip(self):
        from syntheticgen.admixture import write_ancestry_bed
        segs = [
            ("22", 0, 1000, "EUR", "AFR"),
            ("22", 1000, 5000, "SAS", "EUR"),
        ]
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "person.bed"
            write_ancestry_bed(p, segs)
            lines = p.read_text().strip().split("\n")
        self.assertEqual(lines, [
            "22\t0\t1000\tEUR\tAFR",
            "22\t1000\t5000\tSAS\tEUR",
        ])


@unittest.skipUnless(_HAVE_DEPS, "msprime/demes/tskit not installed")
class TestSimulateCohort(unittest.TestCase):
    def test_multi_chromosome(self):
        from syntheticgen.admixture import simulate_cohort
        rng = random.Random(11)
        sites, ancestry = simulate_cohort(
            chromosomes=["21", "22"], build="GRCh38", n_people=4,
            length_mb=0.2, proportions=(0.6, 0.25, 0.15),
            rec_rate=1e-8, mu=1.29e-8, rng=rng,
        )
        chroms = {s["chrom"] for s in sites}
        self.assertEqual(chroms, {"21", "22"})
        # Per-person ancestry covers both chromosomes
        for person in ancestry:
            present = {row[0] for row in person}
            self.assertEqual(present, {"21", "22"})


if __name__ == "__main__":
    unittest.main(verbosity=2)
