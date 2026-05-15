"""Tests for the M13.1 ploidy + PAR helpers in ``syntheticgen.builds``.

These are pure-Python lookups against the BUILDS table — no fixtures
needed, fast to run, and they pin the load-bearing per-chromosome
ploidy contract that M13.3+ will read.

Coverage:

- ``ploidy_for`` returns the correct ploidy for every (chrom, sex,
  pos) combination, including PAR / non-PAR distinction on chrX/chrY.
- ``is_in_par`` recognises the documented PAR ranges and rejects
  positions outside them.
- The function rejects invalid sex inputs.
- Unknown chromosomes default to autosomal (2) — defensive.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.builds import (  # noqa: E402
    BUILDS,
    GRCH37_PAR_REGIONS,
    GRCH38_PAR_REGIONS,
    is_in_par,
    ploidy_for,
)


class IsInParTest(unittest.TestCase):
    """``is_in_par`` is consulted at every chrX/chrY variant in M13.3+
    so its boundary behaviour is load-bearing. Pin both endpoints of
    each PAR range — off-by-ones here would silently flip ploidy on
    edge positions."""

    def test_par1_chrx_grch38_endpoints_inclusive(self):
        # GRCh38 PAR1 on chrX: 10_001 to 2_781_479 (inclusive).
        self.assertTrue(is_in_par("X", 10_001, "GRCh38"))
        self.assertTrue(is_in_par("X", 2_781_479, "GRCh38"))
        # Just outside both ends:
        self.assertFalse(is_in_par("X", 10_000, "GRCh38"))
        self.assertFalse(is_in_par("X", 2_781_480, "GRCh38"))

    def test_par2_chrx_grch38(self):
        # GRCh38 PAR2 on chrX: 155_701_383 to 156_030_895.
        self.assertTrue(is_in_par("X", 155_701_383, "GRCh38"))
        self.assertTrue(is_in_par("X", 156_030_895, "GRCh38"))
        self.assertFalse(is_in_par("X", 155_701_382, "GRCh38"))

    def test_chry_par_endpoints_grch38(self):
        # On GRCh38, chrY PAR1 shares bp coordinates with chrX PAR1
        # (both 10_001-2_781_479) — the PAR1 sequence is identical
        # between the two chromosomes in the assembly. PAR2 sits at
        # the end of each chromosome and therefore lands at
        # different bp on the (shorter) chrY: 56_887_903-57_217_415,
        # versus chrX's ~155.7-156.0 Mb. Lock both ends of both PARs.
        self.assertTrue(is_in_par("Y", 10_001, "GRCh38"))
        self.assertTrue(is_in_par("Y", 2_781_479, "GRCh38"))
        self.assertTrue(is_in_par("Y", 56_887_903, "GRCh38"))
        self.assertTrue(is_in_par("Y", 57_217_415, "GRCh38"))

    def test_non_par_chrx_position_returns_false(self):
        # Mid-X is non-PAR.
        self.assertFalse(is_in_par("X", 80_000_000, "GRCh38"))

    def test_non_xy_chrom_returns_false(self):
        # PAR is an XY concept; autosomes / MT shouldn't trigger.
        self.assertFalse(is_in_par("22", 1_000_000, "GRCh38"))
        self.assertFalse(is_in_par("MT", 1, "GRCh38"))
        self.assertFalse(is_in_par("1", 1, "GRCh38"))

    def test_grch37_par1_chrx_endpoints(self):
        self.assertTrue(is_in_par("X", 60_001, "GRCh37"))
        self.assertTrue(is_in_par("X", 2_699_520, "GRCh37"))
        self.assertFalse(is_in_par("X", 60_000, "GRCh37"))

    def test_par_regions_published_in_builds_dict(self):
        # PAR data must be reachable via the BUILDS dict (downstream
        # M13.3 code threads it through there). Defensive check that
        # the wiring survives any future refactor.
        self.assertIs(BUILDS["GRCh37"]["par_regions"], GRCH37_PAR_REGIONS)
        self.assertIs(BUILDS["GRCh38"]["par_regions"], GRCH38_PAR_REGIONS)


class PloidyForTest(unittest.TestCase):
    """``ploidy_for(chrom, sex, build, pos)`` is the lookup M13.3+
    will use at every variant in the simulator's hot loop. The
    decision tree:

      - Autosomes 1..22 → 2 (both sexes).
      - chrX in female → 2.
      - chrX in male → 1 (non-PAR) or 2 (PAR position).
      - chrY in female → 0 (chromosome absent).
      - chrY in male → 1 (non-PAR) or 2 (PAR position).
      - MT → 1 (both sexes).
      - Unknown chrom → 2 (autosomal default).
    """

    def test_autosomes_always_diploid(self):
        for chrom in ("1", "10", "22"):
            for sex in ("m", "f"):
                self.assertEqual(ploidy_for(chrom, sex), 2,
                                 f"{chrom}/{sex} should be diploid")

    def test_chrx_female_always_diploid(self):
        # No PAR distinction in females — X is always 2.
        self.assertEqual(ploidy_for("X", "f"), 2)
        self.assertEqual(ploidy_for("X", "f", pos=1_000_000), 2)
        self.assertEqual(ploidy_for("X", "f", pos=80_000_000), 2)

    def test_chrx_male_non_par_haploid(self):
        # Non-PAR pos on X: haploid. Includes the no-pos default
        # path (which conservatively returns the non-PAR answer).
        self.assertEqual(ploidy_for("X", "m"), 1)
        self.assertEqual(ploidy_for("X", "m", pos=80_000_000), 1)

    def test_chrx_male_par_diploid(self):
        # PAR pos on X in males: diploid (recombines with Y PAR).
        self.assertEqual(ploidy_for("X", "m", pos=1_000_000), 2)
        # PAR2:
        self.assertEqual(
            ploidy_for("X", "m", pos=155_800_000, build="GRCh38"), 2,
        )

    def test_chry_female_absent(self):
        # Y is absent in females — every position returns 0.
        self.assertEqual(ploidy_for("Y", "f"), 0)
        self.assertEqual(ploidy_for("Y", "f", pos=1_000_000), 0)
        self.assertEqual(ploidy_for("Y", "f", pos=10_001), 0)

    def test_chry_male_non_par_haploid(self):
        self.assertEqual(ploidy_for("Y", "m"), 1)
        self.assertEqual(
            ploidy_for("Y", "m", pos=20_000_000, build="GRCh38"), 1,
        )

    def test_chry_male_par_diploid(self):
        self.assertEqual(
            ploidy_for("Y", "m", pos=1_000_000, build="GRCh38"), 2,
        )

    def test_mt_always_haploid(self):
        # MT is haploid, maternally inherited — sex independent.
        self.assertEqual(ploidy_for("MT", "m"), 1)
        self.assertEqual(ploidy_for("MT", "f"), 1)

    def test_unknown_chrom_defaults_to_diploid(self):
        # Custom contigs shouldn't crash the simulator. Default to 2.
        self.assertEqual(ploidy_for("chr_custom", "m"), 2)

    def test_invalid_sex_raises(self):
        with self.assertRaises(ValueError):
            ploidy_for("1", "male")
        with self.assertRaises(ValueError):
            ploidy_for("X", "")


if __name__ == "__main__":
    unittest.main()
