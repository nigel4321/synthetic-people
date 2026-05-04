"""Unit tests for the --chromosomes spec parser."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.cli import parse_chromosomes


class TestParseChromosomes(unittest.TestCase):
    def test_single(self):
        self.assertEqual(parse_chromosomes("22", "GRCh38"), ["22"])

    def test_comma_list(self):
        self.assertEqual(
            parse_chromosomes("19,20,21,22", "GRCh38"),
            ["19", "20", "21", "22"],
        )

    def test_range(self):
        self.assertEqual(
            parse_chromosomes("1-5", "GRCh38"),
            ["1", "2", "3", "4", "5"],
        )

    def test_full_autosome_range(self):
        self.assertEqual(
            parse_chromosomes("1-22", "GRCh38"),
            [str(i) for i in range(1, 23)],
        )

    def test_mixed_range_and_singletons(self):
        self.assertEqual(
            parse_chromosomes("1-3,5,19-22,X", "GRCh38"),
            ["1", "2", "3", "5", "19", "20", "21", "22", "X"],
        )

    def test_dedupes_preserving_first_occurrence(self):
        self.assertEqual(
            parse_chromosomes("22,1-3,2,22", "GRCh38"),
            ["22", "1", "2", "3"],
        )

    def test_strips_whitespace(self):
        self.assertEqual(
            parse_chromosomes(" 19 , 1 - 3 , X ", "GRCh38"),
            ["19", "1", "2", "3", "X"],
        )

    def test_single_member_range(self):
        self.assertEqual(parse_chromosomes("7-7", "GRCh38"), ["7"])

    def test_non_numeric_range_rejected(self):
        with self.assertRaisesRegex(ValueError, "must be numeric"):
            parse_chromosomes("X-Y", "GRCh38")

    def test_inverted_range_rejected(self):
        with self.assertRaisesRegex(ValueError, "empty"):
            parse_chromosomes("10-5", "GRCh38")

    def test_unknown_chromosome_rejected(self):
        with self.assertRaisesRegex(ValueError, "unknown chromosome"):
            parse_chromosomes("23", "GRCh38")

    def test_empty_spec_rejected(self):
        with self.assertRaisesRegex(ValueError, "empty"):
            parse_chromosomes("", "GRCh38")

    def test_grch37_also_supported(self):
        self.assertEqual(
            parse_chromosomes("1-3", "GRCh37"),
            ["1", "2", "3"],
        )


if __name__ == "__main__":
    unittest.main()
