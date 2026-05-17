"""Distribution checks for syntheticgen.quality.

Runnable either via pytest or bare `python -m unittest`. No numpy
dependency — uses statistics from the stdlib.
"""

from __future__ import annotations

import random
import statistics
import sys
import unittest
from pathlib import Path

# Make the sibling package importable when tests are run directly.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.quality import (  # noqa: E402
    HET_ALT_FRAC,
    ad_from_gt,
    draw_site_quality,
    gq_from_ad,
    poisson,
    sample_lambda,
)


class TestPoisson(unittest.TestCase):
    def test_mean_near_lambda(self):
        rng = random.Random(12345)
        samples = [poisson(30.0, rng) for _ in range(5000)]
        mean = statistics.fmean(samples)
        self.assertAlmostEqual(mean, 30.0, delta=1.0)
        # Variance of Poisson ≈ λ; allow some slack.
        var = statistics.pvariance(samples)
        self.assertAlmostEqual(var, 30.0, delta=3.0)

    def test_floors_at_zero(self):
        rng = random.Random(0)
        samples = [poisson(1.0, rng) for _ in range(2000)]
        self.assertTrue(all(s >= 0 for s in samples))


class TestBiallelicAD(unittest.TestCase):
    """AD with a single ALT — the common case."""

    def test_homref_all_ref(self):
        rng = random.Random(1)
        for _ in range(50):
            ad = ad_from_gt("0|0", 2, 30, rng)
            self.assertEqual(ad, (30, 0))

    def test_homalt_all_alt(self):
        rng = random.Random(2)
        for _ in range(50):
            ad = ad_from_gt("1|1", 2, 40, rng)
            self.assertEqual(ad, (0, 40))

    def test_het_close_to_half_with_ref_bias(self):
        rng = random.Random(3)
        n = 2000
        alt_fracs = []
        for _ in range(n):
            ref, alt = ad_from_gt("0|1", 2, 30, rng)
            self.assertEqual(ref + alt, 30)
            alt_fracs.append(alt / 30)
        mean_alt_frac = statistics.fmean(alt_fracs)
        # Expect centered around HET_ALT_FRAC with modest tolerance.
        self.assertAlmostEqual(mean_alt_frac, HET_ALT_FRAC, delta=0.03)
        # And it should be visibly below 0.5 (the "perfect" het ratio):
        # reference reads are over-represented on real hets.
        self.assertLess(mean_alt_frac, 0.5)

    def test_ad_sum_equals_dp(self):
        rng = random.Random(4)
        for gt in ("0|0", "0|1", "1|0", "1|1"):
            for dp in (0, 1, 10, 30, 100):
                ad = ad_from_gt(gt, 2, dp, rng)
                self.assertEqual(sum(ad), dp, (gt, dp))


class TestMultiallelicAD(unittest.TestCase):
    """AD with two ALTs — `0|2`, `1|2`, etc."""

    def test_zero_two_het_splits_ref_vs_alt2(self):
        rng = random.Random(10)
        dps = []
        for _ in range(1000):
            ad = ad_from_gt("0|2", 3, 40, rng)
            self.assertEqual(sum(ad), 40)
            self.assertEqual(ad[1], 0, "alt1 should have no support in a 0|2")
            dps.append(ad[2] / 40)
        mean_alt2_frac = statistics.fmean(dps)
        # Same ref-bias applies: alt2 should land near HET_ALT_FRAC.
        self.assertAlmostEqual(mean_alt2_frac, HET_ALT_FRAC, delta=0.03)

    def test_one_two_het_splits_50_50_no_ref_bias(self):
        rng = random.Random(11)
        fracs = []
        for _ in range(1000):
            ad = ad_from_gt("1|2", 3, 40, rng)
            self.assertEqual(sum(ad), 40)
            self.assertEqual(ad[0], 0, "ref should have no support in a 1|2")
            fracs.append(ad[1] / 40)
        mean_alt1_frac = statistics.fmean(fracs)
        # Two non-ref hets split 50/50 — no ref bias applies.
        self.assertAlmostEqual(mean_alt1_frac, 0.5, delta=0.03)

    def test_homalt2_all_on_alt2(self):
        rng = random.Random(12)
        for _ in range(50):
            ad = ad_from_gt("2|2", 3, 25, rng)
            self.assertEqual(ad, (0, 0, 25))

    def test_sum_always_dp_across_allele_counts(self):
        rng = random.Random(13)
        for n_alleles in (2, 3, 4):
            for gt in ("0|0", "0|1", "1|2", "2|2"):
                for dp in (0, 5, 30):
                    ad = ad_from_gt(gt, n_alleles, dp, rng)
                    self.assertEqual(len(ad), n_alleles)
                    self.assertEqual(sum(ad), dp, (gt, n_alleles, dp))


class TestGQ(unittest.TestCase):
    def test_range_clamped(self):
        # High-support cases should go high, contradictions should go low.
        self.assertGreater(gq_from_ad("1|1", (0, 40)), 80)
        self.assertGreater(gq_from_ad("0|0", (40, 0)), 80)
        # Hom-alt call with zero alt reads → disagreement → low GQ.
        self.assertLess(gq_from_ad("1|1", (40, 0)), 10)
        # Depth=0 → GQ=0.
        self.assertEqual(gq_from_ad("0|1", (0, 0)), 0)

    def test_het_peaks_near_half(self):
        # GQ for a het call should be highest when AD is roughly 50/50
        # and lower at the tails.
        center = gq_from_ad("0|1", (20, 20))
        skewed = gq_from_ad("0|1", (35, 5))
        self.assertGreater(center, skewed)

    def test_multiallelic_support(self):
        # Het 1|2 with perfect split between alt1 and alt2 and no ref
        # contamination should have high GQ.
        good = gq_from_ad("1|2", (0, 20, 20))
        # Same call but with noisy ref reads should drop GQ.
        noisy = gq_from_ad("1|2", (10, 15, 15))
        self.assertGreater(good, noisy)
        self.assertGreater(good, 80)

    def test_bounds(self):
        for ref in range(0, 30, 5):
            for alt in range(0, 30, 5):
                for gt in ("0|0", "0|1", "1|1"):
                    gq = gq_from_ad(gt, (ref, alt))
                    self.assertGreaterEqual(gq, 0)
                    self.assertLessEqual(gq, 99)


class TestHaploidGtSupport(unittest.TestCase):
    """M13.3 review (Copilot PR #107): quality helpers must treat
    haploid GTs ("0" / "1") as homozygous-equivalent for the AD draw
    and GQ recompute. Pre-fix _parse_gt_alleles fell back to (0, 0)
    for any non-2-part input, so a haploid alt call ("1") emitted by
    write_person_vcf got AD that looked hom-ref — inconsistent with
    the called GT."""

    def test_parse_gt_alleles_haploid_alt(self):
        from syntheticgen.quality import _parse_gt_alleles
        # Single token "1" must parse as (1, 1), not (0, 0).
        self.assertEqual(_parse_gt_alleles("1"), (1, 1))

    def test_parse_gt_alleles_haploid_ref(self):
        from syntheticgen.quality import _parse_gt_alleles
        self.assertEqual(_parse_gt_alleles("0"), (0, 0))

    def test_parse_gt_alleles_haploid_multiallelic(self):
        from syntheticgen.quality import _parse_gt_alleles
        # Haploid alt-2 (multi-allelic at a haploid locus).
        self.assertEqual(_parse_gt_alleles("2"), (2, 2))

    def test_haploid_alt_ad_looks_hom_alt(self):
        # Pre-fix: AD for haploid "1" looked hom-ref (all reads on
        # the ref slot). Now it looks hom-alt — every read on the
        # alt slot — consistent with the called GT.
        rng = random.Random(7)
        for _ in range(50):
            ad = ad_from_gt("1", 2, 30, rng)
            self.assertEqual(ad[0], 0,
                             "haploid '1' should have no ref reads")
            self.assertEqual(ad[1], 30,
                             "haploid '1' should put all reads on alt")

    def test_haploid_ref_ad_looks_hom_ref(self):
        rng = random.Random(7)
        for _ in range(50):
            ad = ad_from_gt("0", 2, 30, rng)
            self.assertEqual(ad[0], 30)
            self.assertEqual(ad[1], 0)

    def test_haploid_alt_gq_high_when_reads_agree(self):
        # GQ for a haploid alt with all reads on alt should be high
        # (confident call). Pre-fix it would have been 0 because the
        # support formula used ad[0] (ref slot) on a (0, 0) parse.
        gq = gq_from_ad("1", (0, 40))
        self.assertGreater(gq, 80)

    def test_haploid_alt_gq_low_when_reads_disagree(self):
        # The mirror: 0 alt reads on a haploid "1" call should be
        # near 0 (high disagreement).
        gq = gq_from_ad("1", (40, 0))
        self.assertLess(gq, 10)


class TestDrawSiteQuality(unittest.TestCase):
    def test_tuple_consistency(self):
        rng = random.Random(7)
        for _ in range(200):
            for gt in ("0|0", "0|1", "1|1"):
                dp, ad, gq = draw_site_quality(gt, 2, 30.0, rng)
                self.assertEqual(len(ad), 2)
                self.assertEqual(dp, sum(ad))
                self.assertGreaterEqual(gq, 0)
                self.assertLessEqual(gq, 99)

    def test_multiallelic_consistency(self):
        rng = random.Random(8)
        for _ in range(200):
            for gt in ("0|1", "0|2", "1|2", "2|2"):
                dp, ad, gq = draw_site_quality(gt, 3, 30.0, rng)
                self.assertEqual(len(ad), 3)
                self.assertEqual(dp, sum(ad))
                self.assertGreaterEqual(gq, 0)
                self.assertLessEqual(gq, 99)

    def test_cohort_dp_mean(self):
        rng = random.Random(11)
        dps = []
        for _ in range(500):
            lam = sample_lambda(30.0, 3.0, rng)
            for _ in range(10):
                dp, _, _ = draw_site_quality("0|1", 2, lam, rng)
                dps.append(dp)
        mean_dp = statistics.fmean(dps)
        self.assertAlmostEqual(mean_dp, 30.0, delta=1.5)

    def test_sample_lambda_clamp(self):
        rng = random.Random(13)
        for _ in range(100):
            lam = sample_lambda(30.0, 3.0, rng)
            self.assertGreaterEqual(lam, 5.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
