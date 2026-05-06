"""Tests for the coalescent backbone (M5).

Skipped gracefully when msprime / stdpopsim aren't installed, so the rest
of the suite still passes in a stdlib-only environment.
"""

from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

try:
    import msprime  # noqa: F401
    import stdpopsim  # noqa: F401
    _HAVE_DEPS = True
except ImportError:
    _HAVE_DEPS = False


@unittest.skipUnless(_HAVE_DEPS, "msprime/stdpopsim not installed")
class TestSimulateChromosome(unittest.TestCase):
    def setUp(self):
        from syntheticgen.coalescent import simulate_chromosome  # noqa
        self.simulate_chromosome = simulate_chromosome

    def _run(self, seed: int = 42, demo_model: str | None = None,
             length_mb: float = 0.3, n_people: int = 6):
        # Default to the lighter no-demo path for speed; the stdpopsim
        # path is exercised once below.
        rng = random.Random(seed)
        return self.simulate_chromosome(
            chrom="22", build="GRCh38", n_people=n_people,
            length_mb=length_mb, demo_model=demo_model,
            population="CEU", rec_rate=1e-8, mu=1.29e-8, rng=rng,
        )

    def test_returns_sites_with_correct_shape(self):
        sites = self._run()
        self.assertGreater(len(sites), 20)
        for s in sites:
            self.assertEqual(s["chrom"], "22")
            self.assertEqual(len(s["alts"]), 1)
            self.assertEqual(len(s["acs"]), 1)
            self.assertEqual(len(s["afs"]), 1)
            # Phase 5c: site shape carries n_haplotypes + sparse
            # carriers rather than a dense gts list. n_people=6 →
            # 12 haplotypes.
            self.assertEqual(s["n_haplotypes"], 12)
            self.assertIn("carriers", s)

    def test_positions_strictly_increasing(self):
        sites = self._run()
        for i in range(len(sites) - 1):
            self.assertLess(sites[i]["pos"], sites[i + 1]["pos"],
                            f"non-monotone at i={i}")

    def test_realised_ac_matches_declared(self):
        sites = self._run()
        for s in sites:
            # Sparse carriers: every entry is a non-zero allele, so
            # AC equals the carriers count for biallelic sites.
            realised = sum(
                1 for _, allele in s["carriers"] if allele >= 1
            )
            self.assertEqual(realised, s["acs"][0],
                             f"AC mismatch at {s['pos']}")

    def test_no_fixed_sites(self):
        """Every output site should be genuinely variable."""
        sites = self._run()
        for s in sites:
            ac = s["acs"][0]
            self.assertGreater(ac, 0)
            self.assertLess(ac, 2 * 6)

    def test_reproducible_under_seed(self):
        a = self._run(seed=42)
        b = self._run(seed=42)
        self.assertEqual(len(a), len(b))
        for s1, s2 in zip(a, b):
            self.assertEqual(s1["pos"], s2["pos"])
            self.assertEqual(s1["acs"], s2["acs"])
            # Phase 5c: compare carriers (sparse) instead of gts.
            self.assertEqual(s1["carriers"], s2["carriers"])

    def test_different_seeds_give_different_output(self):
        a = self._run(seed=1)
        b = self._run(seed=2)
        # At least positions or ACs should differ.
        pos_a = [s["pos"] for s in a]
        pos_b = [s["pos"] for s in b]
        self.assertNotEqual(pos_a, pos_b)

    def test_titv_in_range(self):
        """With the Ti/Tv calibrator (default target 2.1), the ratio of
        transitions to transversions should land well above 1."""
        from syntheticgen.titv import is_transition
        # Bigger sample → tighter CI; 1 Mb gives ~500-2000 SNVs.
        sites = self._run(seed=42, length_mb=1.0, n_people=10)
        ti = sum(1 for s in sites
                 if is_transition(s["ref"], s["alts"][0]))
        tv = len(sites) - ti
        self.assertGreater(tv, 0)
        ratio = ti / tv
        # Target is 2.1; accept [1.7, 2.6] to absorb small-sample noise.
        self.assertGreater(ratio, 1.7)
        self.assertLess(ratio, 2.6)

    def test_unknown_chromosome_raises(self):
        rng = random.Random(0)
        with self.assertRaises(ValueError):
            self.simulate_chromosome(
                chrom="42", build="GRCh38", n_people=4, length_mb=0.1,
                demo_model=None, population="CEU",
                rec_rate=1e-8, mu=1.29e-8, rng=rng,
            )

    def test_stdpopsim_demography_path(self):
        """One end-to-end run through the stdpopsim engine."""
        sites = self._run(seed=42, demo_model="OutOfAfrica_3G09",
                          length_mb=0.3, n_people=5)
        self.assertGreater(len(sites), 5)
        # Positions are still integer-valued; sparse carriers carry
        # only allele index 1 (multi-allelic JC69 sites are dropped
        # upstream so the cohort BCF stays biallelic-spec-clean).
        for s in sites:
            self.assertIsInstance(s["pos"], int)
            for hap_idx, allele in s["carriers"]:
                self.assertEqual(allele, 1)
                self.assertLess(hap_idx, s["n_haplotypes"])


@unittest.skipUnless(_HAVE_DEPS, "msprime/stdpopsim not installed")
class TestSimulateCohort(unittest.TestCase):
    def test_multi_chromosome(self):
        from syntheticgen.coalescent import simulate_cohort
        rng = random.Random(11)
        sites = simulate_cohort(
            chromosomes=["21", "22"], build="GRCh38", n_people=4,
            length_mb=0.2, demo_model=None, population="CEU",
            rec_rate=1e-8, mu=1.29e-8, rng=rng,
        )
        chroms = {s["chrom"] for s in sites}
        self.assertEqual(chroms, {"21", "22"})


if __name__ == "__main__":
    unittest.main(verbosity=2)
