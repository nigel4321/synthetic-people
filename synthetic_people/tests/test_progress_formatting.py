"""Unit tests for the progress-log duration formatter.

The progress logs in cli.py print elapsed/eta values that can run to
multi-hour totals on big cohort runs (e.g. ``elapsed 9913s, eta
27970s``); the helper renders them as ``Hh Mm Ss`` so the user can
read them at a glance instead of doing seconds-to-hours arithmetic
in their head.
"""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.cli import _format_duration


class FormatDurationTest(unittest.TestCase):

    def test_seconds_only_for_under_one_minute(self):
        self.assertEqual(_format_duration(0), "0s")
        self.assertEqual(_format_duration(1), "1s")
        self.assertEqual(_format_duration(45), "45s")
        self.assertEqual(_format_duration(59), "59s")

    def test_minutes_and_seconds_for_under_one_hour(self):
        self.assertEqual(_format_duration(60), "1m 0s")
        self.assertEqual(_format_duration(90), "1m 30s")
        self.assertEqual(_format_duration(3599), "59m 59s")

    def test_hours_minutes_seconds(self):
        self.assertEqual(_format_duration(3600), "1h 0m 0s")
        # Real value pulled straight from a recent fanout run so the
        # test pins the human-friendly transformation we wanted.
        self.assertEqual(_format_duration(9913), "2h 45m 13s")
        self.assertEqual(_format_duration(27970), "7h 46m 10s")
        # Long runs (e.g. n=100k stretch target): make sure we keep
        # rendering hours rather than rolling into days.
        self.assertEqual(_format_duration(100_000), "27h 46m 40s")

    def test_rounds_to_nearest_second(self):
        # We track elapsed/eta as floats; the formatter rounds rather
        # than truncates so ``29.9s`` doesn't display as ``29s``.
        self.assertEqual(_format_duration(29.9), "30s")
        self.assertEqual(_format_duration(0.4), "0s")
        self.assertEqual(_format_duration(0.6), "1s")

    def test_infinity_renders_as_question_mark(self):
        # The eta computation returns ``inf`` when rate=0 (no progress
        # yet); callers used to special-case this themselves with a
        # guard branch — the helper now owns the convention.
        self.assertEqual(_format_duration(float("inf")), "?")

    def test_nan_renders_as_question_mark(self):
        self.assertEqual(_format_duration(float("nan")), "?")

    def test_negative_clamps_to_zero(self):
        # Defensive — a clock-skew or off-by-one elapsed shouldn't
        # surface as ``-1s`` in user-facing output.
        self.assertEqual(_format_duration(-5), "0s")


if __name__ == "__main__":
    unittest.main()
