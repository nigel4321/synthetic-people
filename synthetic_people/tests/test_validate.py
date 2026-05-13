"""Tests for the M10 validation suite.

Math primitives (per-record classification, dosage parsing, AF
binning, het/hom counts, r²) run with numpy installed (a hard dep
since M2). PCA / plot smoke-tests are gated on `sklearn` and
`matplotlib` so the file still loads in stripped environments.
"""

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import sklearn  # noqa: F401
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

try:
    import matplotlib  # noqa: F401
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

from syntheticgen.validate import (
    Record,
    SampleStats,
    _classify_record,
    _gt_dosage,
    _is_dropout,
    _parse_info,
    aggregate_indel_lengths,
    aggregate_sv_summary,
    check_ref_against_fasta,
    cohort_chrom_stats,
    cohort_overlay_density,
    het_hom_ratio,
    titv_from_stats,
)


class TestParseInfo(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(_parse_info(""), {})
        self.assertEqual(_parse_info("."), {})

    def test_single_kv(self):
        self.assertEqual(_parse_info("AC=1"), {"AC": "1"})

    def test_multiple_kv(self):
        d = _parse_info("AC=1;AN=2;AF=0.5")
        self.assertEqual(d, {"AC": "1", "AN": "2", "AF": "0.5"})

    def test_flag(self):
        d = _parse_info("HIGHLIGHT;AC=1")
        self.assertTrue(d["HIGHLIGHT"])
        self.assertEqual(d["AC"], "1")


class TestClassifyRecord(unittest.TestCase):
    def _rec(self, ref, alt):
        return Record(chrom="22", pos=1, ref=ref, alt=alt, gt="0|1",
                      dp=10, gq=99, ad_ref=5, ad_alt=5, info={})

    def test_snv(self):
        self.assertEqual(_classify_record(self._rec("A", "G")), "snv")

    def test_indel_insertion(self):
        self.assertEqual(_classify_record(self._rec("A", "AGT")), "indel")

    def test_indel_deletion(self):
        self.assertEqual(_classify_record(self._rec("ACGT", "A")), "indel")

    def test_sv_symbolic_alt(self):
        for alt in ("<DEL>", "<DUP>", "<INV>", "<INS>"):
            self.assertEqual(_classify_record(self._rec("A", alt)), "sv")


class TestGtDosage(unittest.TestCase):
    def test_hom_ref(self):
        self.assertEqual(_gt_dosage("0|0"), 0)
        self.assertEqual(_gt_dosage("0/0"), 0)

    def test_het(self):
        self.assertEqual(_gt_dosage("0|1"), 1)
        self.assertEqual(_gt_dosage("1|0"), 1)

    def test_hom_alt(self):
        self.assertEqual(_gt_dosage("1|1"), 2)

    def test_multiallelic_het_treated_as_two_alts(self):
        # 1|2 has both alt alleles set — per-ALT splitting upstream
        # means each Record sees one ALT, but the dosage helper still
        # returns "non-zero alleles", so 1|2 reports 2.
        self.assertEqual(_gt_dosage("1|2"), 2)

    def test_missing_returns_negative(self):
        self.assertEqual(_gt_dosage("./."), -1)
        self.assertEqual(_gt_dosage(".|."), -1)
        self.assertEqual(_gt_dosage("./1"), -1)


class TestIsDropout(unittest.TestCase):
    def test_clean_calls(self):
        for gt in ("0|0", "0|1", "1|1", "1|2"):
            self.assertFalse(_is_dropout(gt))

    def test_dropouts(self):
        for gt in ("./.", ".|.", "./0"):
            self.assertTrue(_is_dropout(gt))


class TestSampleStatsAggregation(unittest.TestCase):
    def _stats(self, ti, tv, het, hom_alt):
        s = SampleStats(name="x")
        s.n_ti, s.n_tv, s.n_het, s.n_hom_alt = ti, tv, het, hom_alt
        return s

    def test_titv_from_stats(self):
        a = self._stats(ti=21, tv=10, het=0, hom_alt=0)
        b = self._stats(ti=15, tv=10, het=0, hom_alt=0)
        self.assertAlmostEqual(titv_from_stats([a, b]), (21 + 15) / 20)

    def test_titv_zero_tv(self):
        a = self._stats(ti=5, tv=0, het=0, hom_alt=0)
        self.assertEqual(titv_from_stats([a]), float("inf"))

    def test_titv_empty(self):
        self.assertEqual(titv_from_stats([]), 0.0)

    def test_het_hom_ratio(self):
        s = self._stats(ti=0, tv=0, het=20, hom_alt=10)
        self.assertAlmostEqual(het_hom_ratio(s), 2.0)

    def test_het_hom_zero_hom(self):
        s = self._stats(ti=0, tv=0, het=5, hom_alt=0)
        self.assertEqual(het_hom_ratio(s), float("inf"))

    def test_het_hom_zero_both(self):
        s = self._stats(ti=0, tv=0, het=0, hom_alt=0)
        self.assertEqual(het_hom_ratio(s), 0.0)

    def test_aggregate_indel_lengths(self):
        a = SampleStats(name="a")
        a.indel_lengths = [-3, -1, 2]
        b = SampleStats(name="b")
        b.indel_lengths = [2, 5]
        h = aggregate_indel_lengths([a, b])
        self.assertEqual(h, {-3: 1, -1: 1, 2: 2, 5: 1})

    def test_aggregate_sv_summary(self):
        from collections import defaultdict
        a = SampleStats(name="a")
        a.sv_by_type = defaultdict(int, {"DEL": 3, "DUP": 1})
        b = SampleStats(name="b")
        b.sv_by_type = defaultdict(int, {"DEL": 2, "INV": 1})
        out = aggregate_sv_summary([a, b])
        self.assertEqual(out, {"DEL": 5, "DUP": 1, "INV": 1})


class TestCohortChromStats(unittest.TestCase):
    """Tier 1: per-chromosome breakouts surface chrom-specific
    regressions that cohort-wide aggregates hide."""

    def _stats_with_chroms(self, chrom_counts: dict) -> SampleStats:
        """``chrom_counts`` is ``{chrom: {n_records, n_ti, n_tv, ...}}``
        — populates a SampleStats with the requested per-chrom
        buckets, skipping the iter_records / bcftools path."""
        s = SampleStats(name="x")
        for chrom, fields in chrom_counts.items():
            for k, v in fields.items():
                s.by_chrom[chrom][k] = v
        return s

    def test_per_chrom_counts_aggregate_across_samples(self):
        a = self._stats_with_chroms(
            {"22": {"n_records": 10, "n_snv": 8, "n_ti": 5, "n_tv": 3}}
        )
        b = self._stats_with_chroms(
            {"22": {"n_records": 4, "n_snv": 3, "n_ti": 2, "n_tv": 1},
             "X":  {"n_records": 7, "n_snv": 7, "n_ti": 5, "n_tv": 2}}
        )
        out = cohort_chrom_stats([a, b])
        self.assertEqual(out["22"]["n_records"], 14)
        self.assertEqual(out["22"]["n_snv"], 11)
        self.assertEqual(out["22"]["n_ti"], 7)
        self.assertEqual(out["22"]["n_tv"], 4)
        self.assertEqual(out["X"]["n_records"], 7)

    def test_per_chrom_titv_attached_to_each_row(self):
        a = self._stats_with_chroms(
            {"22": {"n_ti": 21, "n_tv": 10}}
        )
        out = cohort_chrom_stats([a])
        self.assertAlmostEqual(out["22"]["titv"], 2.1)

    def test_per_chrom_titv_handles_zero_tv(self):
        a = self._stats_with_chroms(
            {"Y": {"n_ti": 5, "n_tv": 0}}
        )
        out = cohort_chrom_stats([a])
        self.assertEqual(out["Y"]["titv"], float("inf"))

    def test_chrom_order_canonical(self):
        # Mix integer chroms with X/Y/MT; verify the order is
        # 1-22, X, Y, MT — the standard VCF convention.
        a = self._stats_with_chroms({
            "MT": {"n_records": 1},
            "Y": {"n_records": 1},
            "X": {"n_records": 1},
            "2": {"n_records": 1},
            "22": {"n_records": 1},
            "10": {"n_records": 1},
        })
        out = cohort_chrom_stats([a])
        self.assertEqual(
            list(out.keys()), ["2", "10", "22", "X", "Y", "MT"],
        )

    def test_chr_prefix_sorts_with_unprefixed(self):
        # "chr22" should sort like "22"; both common in real VCFs.
        a = self._stats_with_chroms({
            "chr22": {"n_records": 1},
            "chr2": {"n_records": 1},
        })
        out = cohort_chrom_stats([a])
        self.assertEqual(list(out.keys()), ["chr2", "chr22"])

    def test_mixed_prefix_has_deterministic_tertiary_tiebreak(self):
        # Regression for PR #74 review #3: when ``"22"`` and
        # ``"chr22"`` are both present, the primary sort key is
        # identical (both → ``(0, 22, ...)``). The tertiary
        # tie-breaker (raw chrom string) must give a deterministic
        # order — lexicographic, so ``"22"`` < ``"chr22"``.
        a = self._stats_with_chroms({
            "chr22": {"n_records": 1},
            "22": {"n_records": 1},
        })
        out = cohort_chrom_stats([a])
        self.assertEqual(list(out.keys()), ["22", "chr22"])
        # And the reverse insertion order must produce the same
        # final order — the tie-breaker dominates over insertion.
        b = self._stats_with_chroms({
            "22": {"n_records": 1},
            "chr22": {"n_records": 1},
        })
        out2 = cohort_chrom_stats([b])
        self.assertEqual(list(out2.keys()), ["22", "chr22"])

    def test_cohort_overlay_density_handles_generator_input(self):
        # Regression for PR #74 review #1: ``cohort_overlay_density``
        # previously made four separate ``sum()`` calls, exhausting
        # a generator after the first. Now it's a single-pass loop;
        # passing a one-shot generator must yield correct counts.
        from syntheticgen.validate import cohort_overlay_density
        s1 = SampleStats(name="a")
        s1.n_records = 10
        s1.n_with_rs = 3
        s1.n_with_clnsig = 1
        s1.n_with_cosmic_id = 0
        s2 = SampleStats(name="b")
        s2.n_records = 20
        s2.n_with_rs = 5
        s2.n_with_clnsig = 2
        s2.n_with_cosmic_id = 1
        # ``iter(...)`` makes it a true one-shot generator —
        # exhausted after a single pass.
        out = cohort_overlay_density(iter([s1, s2]))
        self.assertEqual(out["n_records"], 30)
        self.assertEqual(out["rsid"]["n"], 8)
        self.assertEqual(out["clinvar"]["n"], 3)
        self.assertEqual(out["cosmic"]["n"], 1)


class TestJsonable(unittest.TestCase):
    """Regression for PR #74 review #5. ``_jsonable`` in
    ``validate_batch.py`` must replace non-finite floats with
    ``None`` so the emitted ``summary.json`` is RFC-8259 valid
    (``Infinity`` / ``NaN`` are forbidden by strict JSON parsers,
    and per-chrom Ti/Tv goes to ``inf`` when ``tv == 0``)."""

    def setUp(self):
        # ``validate_batch.py`` lives alongside the ``syntheticgen``
        # package in ``synthetic_people/``. The top-of-file
        # ``sys.path.insert(0, str(Path(__file__).resolve().parent.parent))``
        # adds that directory to ``sys.path``, so it's importable
        # by its bare module name.
        import validate_batch  # noqa: F401
        self._jsonable = validate_batch._jsonable

    def test_inf_replaced_with_none(self):
        self.assertIsNone(self._jsonable(float("inf")))
        self.assertIsNone(self._jsonable(float("-inf")))

    def test_nan_replaced_with_none(self):
        self.assertIsNone(self._jsonable(float("nan")))

    def test_finite_floats_pass_through(self):
        for v in (0.0, -1.5, 3.14, 1e308):
            self.assertEqual(self._jsonable(v), v)

    def test_non_float_types_unchanged(self):
        # ints, strings, bools, None — none of these are
        # non-finite floats, so they pass through unchanged.
        self.assertEqual(self._jsonable(42), 42)
        self.assertEqual(self._jsonable("hello"), "hello")
        self.assertIs(self._jsonable(True), True)
        self.assertIsNone(self._jsonable(None))

    def test_walks_nested_dicts(self):
        # Mimics the shape of the cohort_chrom_stats output.
        inp = {
            "22": {"n_records": 100, "titv": 2.1},
            "Y":  {"n_records": 5, "titv": float("inf")},
            "X":  {"n_records": 50, "titv": float("nan")},
        }
        out = self._jsonable(inp)
        self.assertEqual(out["22"]["titv"], 2.1)
        self.assertIsNone(out["Y"]["titv"])
        self.assertIsNone(out["X"]["titv"])

    def test_walks_nested_lists(self):
        # Mimics the ld_decay list-of-dicts shape.
        inp = [{"mean_r2": 0.5}, {"mean_r2": float("nan")}]
        out = self._jsonable(inp)
        self.assertEqual(out[0]["mean_r2"], 0.5)
        self.assertIsNone(out[1]["mean_r2"])

    def test_round_trips_through_strict_json(self):
        # End-to-end: a sanitised payload with an inf-typed value
        # must serialise cleanly under strict JSON (``allow_nan=False``).
        import json
        inp = {"chrom_stats": {"Y": {"titv": float("inf")}}}
        out = self._jsonable(inp)
        # allow_nan=False raises on inf/nan; success proves
        # _jsonable replaced everything non-finite.
        json.dumps(out, allow_nan=False)


class TestCohortOverlayDensity(unittest.TestCase):
    """Tier 1: realised overlay-density counts let the validator
    compare what landed in the VCFs against what the manifest
    requested — catches drift in the overlay pipeline."""

    def _stats(self, name: str, n: int, rs: int = 0,
               cln: int = 0, cos: int = 0) -> SampleStats:
        s = SampleStats(name=name)
        s.n_records = n
        s.n_with_rs = rs
        s.n_with_clnsig = cln
        s.n_with_cosmic_id = cos
        return s

    def test_aggregates_counts_and_fractions_across_samples(self):
        a = self._stats("a", n=100, rs=20, cln=1, cos=0)
        b = self._stats("b", n=100, rs=15, cln=2, cos=0)
        out = cohort_overlay_density([a, b])
        self.assertEqual(out["n_records"], 200)
        self.assertEqual(out["rsid"]["n"], 35)
        self.assertAlmostEqual(out["rsid"]["fraction"], 0.175)
        self.assertEqual(out["clinvar"]["n"], 3)
        self.assertAlmostEqual(out["clinvar"]["fraction"], 0.015)
        self.assertEqual(out["cosmic"]["n"], 0)
        self.assertEqual(out["cosmic"]["fraction"], 0.0)

    def test_zero_records_returns_zero_fractions(self):
        # Defensive: no division-by-zero on empty cohorts.
        out = cohort_overlay_density([self._stats("e", n=0)])
        self.assertEqual(out["n_records"], 0)
        self.assertEqual(out["rsid"]["fraction"], 0.0)
        self.assertEqual(out["clinvar"]["fraction"], 0.0)
        self.assertEqual(out["cosmic"]["fraction"], 0.0)


class TestCheckRefAgainstFasta(unittest.TestCase):
    """Tier 1: REF-matches-FASTA gate. Today's synthetic output
    uses fabricated REF so the gate fails on every record — the
    test below verifies the helper returns a structured failure
    rather than crashing, so downstream consumers can display it."""

    def _mock_popen(self, stderr_bytes: bytes, returncode: int = 0):
        """Build a MagicMock that simulates ``subprocess.Popen``'s
        contract for the streaming-stderr ref-check flow.

        ``stderr`` is iterated line-by-line by the helper, so the
        mock exposes a byte-stream iterable rather than a single
        captured blob. ``wait`` returns the configured returncode.
        """
        from unittest.mock import MagicMock
        import io
        fake = MagicMock()
        fake.stderr = io.BytesIO(stderr_bytes)
        fake.wait = MagicMock(return_value=returncode)
        return fake

    def test_missing_bcftools_returns_errored(self):
        # When bcftools isn't on PATH, the helper must surface a
        # structured ``errored=True`` result rather than crash. Mock
        # Popen to raise FileNotFoundError as if the executable
        # couldn't be located.
        from unittest.mock import patch
        with patch(
            "syntheticgen.validate.subprocess.Popen",
            side_effect=FileNotFoundError("bcftools"),
        ):
            r = check_ref_against_fasta(
                Path("/nonexistent.vcf.gz"), Path("/nonexistent.fa"),
            )
        self.assertTrue(r["errored"])
        self.assertFalse(r["passed"])
        self.assertEqual(r["mismatches"], 0)
        self.assertIn("bcftools", r["stderr_tail"])

    def test_mismatches_counted_from_stderr(self):
        # When bcftools writes REF_MISMATCH warnings, the helper
        # streams stderr line-by-line and counts them. Mock Popen so
        # ``proc.stderr`` is an iterable byte-stream of the canonical
        # REF_MISMATCH-line format.
        from unittest.mock import patch
        stderr = (
            b"REF_MISMATCH\t20\t12345\tA\tG\n"
            b"REF_MISMATCH\t20\t67890\tT\tC\n"
            b"Lines reformatted: 2\n"
        )
        with patch(
            "syntheticgen.validate.subprocess.Popen",
            return_value=self._mock_popen(stderr, returncode=0),
        ):
            r = check_ref_against_fasta(
                Path("/nonexistent.vcf.gz"), Path("/nonexistent.fa"),
            )
        self.assertEqual(r["mismatches"], 2)
        # passed=False because mismatches > 0, even though
        # returncode==0 (bcftools warned, didn't error).
        self.assertFalse(r["passed"])
        self.assertFalse(r["errored"])

    def test_clean_pass(self):
        # No REF_MISMATCH lines + zero exit = passed.
        from unittest.mock import patch
        with patch(
            "syntheticgen.validate.subprocess.Popen",
            return_value=self._mock_popen(
                b"Lines total: 1000\n", returncode=0,
            ),
        ):
            r = check_ref_against_fasta(
                Path("/nonexistent.vcf.gz"), Path("/nonexistent.fa"),
            )
        self.assertTrue(r["passed"])
        self.assertEqual(r["mismatches"], 0)
        self.assertFalse(r["errored"])

    def test_stderr_tail_bounded_under_high_mismatch_volume(self):
        # Regression for PR #75 review #1: ``--check-ref w`` emits one
        # REF_MISMATCH line per mismatched record; on today's
        # fabricated-REF synthetic output that's millions of lines per
        # chrom. The helper must stream stderr rather than capturing
        # it in full — feeding 100K mismatch lines must yield the full
        # count without retaining the full stderr in memory.
        from unittest.mock import patch
        n_lines = 100_000
        stderr = b"".join(
            f"REF_MISMATCH\t1\t{i}\tA\tG\n".encode()
            for i in range(n_lines)
        )
        with patch(
            "syntheticgen.validate.subprocess.Popen",
            return_value=self._mock_popen(stderr, returncode=0),
        ):
            r = check_ref_against_fasta(
                Path("/nonexistent.vcf.gz"), Path("/nonexistent.fa"),
            )
        # All mismatches counted.
        self.assertEqual(r["mismatches"], n_lines)
        # Diagnostic tail is bounded — not the full stderr blob.
        # The helper truncates to ~1500 bytes for human-readable
        # output; verify we're not retaining megabytes.
        self.assertLess(len(r["stderr_tail"]), 2000)
        # Tail contains the LAST mismatch lines, not the first
        # (deque-with-maxlen semantics).
        self.assertIn(f"\t{n_lines - 1}\t", r["stderr_tail"])


class TestSummariseVcfOverlayCounters(unittest.TestCase):
    """Tier 1: confirm ``summarise_vcf`` increments the overlay-
    density counters on ``SampleStats`` for records that carry
    INFO/RS, INFO/CLNSIG, or INFO/COSMIC_ID — and ignores empty,
    dotted, or absent values.

    Exercises the *production* counter logic by patching
    ``iter_records`` to yield ``Record`` fixtures rather than
    re-implementing the predicate inline (PR #74 review caught
    that the original tests passed even if ``summarise_vcf`` was
    broken). End-to-end through ``SampleStats``.
    """

    def _record(self, info_str: str) -> Record:
        return Record(
            chrom="22", pos=100, ref="A", alt="C", gt="0|1",
            dp=30, gq=40, ad_ref=15, ad_alt=15,
            info=_parse_info(info_str),
        )

    def _summarise(self, records: list):
        # Patch iter_records to yield the fixture records,
        # bypassing the bcftools subprocess in summarise_vcf.
        from syntheticgen.validate import summarise_vcf
        with patch(
            "syntheticgen.validate.iter_records",
            return_value=iter(records),
        ):
            return summarise_vcf(Path("/nonexistent.vcf.gz"), name="t")

    def test_info_rs_increments_counter(self):
        stats = self._summarise([self._record("RS=12345")])
        self.assertEqual(stats.n_with_rs, 1)
        self.assertEqual(stats.n_with_clnsig, 0)
        self.assertEqual(stats.n_with_cosmic_id, 0)

    def test_info_rs_empty_not_counted(self):
        # "." and "" both mean "no value" in VCF INFO conventions
        # — neither should increment the rsID counter.
        for empty in (".", ""):
            stats = self._summarise([self._record(f"RS={empty}")])
            self.assertEqual(stats.n_with_rs, 0)

    def test_info_rs_absent_not_counted(self):
        stats = self._summarise([self._record("AC=1;AN=2")])
        self.assertEqual(stats.n_with_rs, 0)

    def test_info_flag_form_counted(self):
        # ``INFO=FLAG`` (no "=") is parsed as ``True`` by
        # ``_parse_info``; the counter must accept it as set.
        stats = self._summarise([self._record("CLNSIG;AC=1")])
        self.assertEqual(stats.n_with_clnsig, 1)

    def test_clinvar_and_cosmic_counters(self):
        # Three records: rsID-only, ClinVar-only, COSMIC-only —
        # each should increment exactly one counter.
        stats = self._summarise([
            self._record("RS=12345"),
            self._record("CLNSIG=Pathogenic"),
            self._record("COSMIC_ID=COSV12345"),
        ])
        self.assertEqual(stats.n_with_rs, 1)
        self.assertEqual(stats.n_with_clnsig, 1)
        self.assertEqual(stats.n_with_cosmic_id, 1)
        self.assertEqual(stats.n_records, 3)

    def test_multi_marker_record_counts_each_channel(self):
        # A record carrying both rs and CLNSIG (e.g. a ClinVar
        # entry with a known rsID) must increment both counters.
        stats = self._summarise([
            self._record("RS=12345;CLNSIG=Pathogenic"),
        ])
        self.assertEqual(stats.n_with_rs, 1)
        self.assertEqual(stats.n_with_clnsig, 1)


@unittest.skipUnless(HAS_NUMPY, "numpy required")
class TestAfHistogram(unittest.TestCase):
    def test_distribution_lands_in_correct_bins(self):
        from syntheticgen.validate import af_histogram

        s = SampleStats(name="x")
        s.af_values = [0.0, 0.05, 0.5, 0.95, 1.0]
        edges, counts = af_histogram([s], n_bins=10)
        # 10 bins of width 0.1 over [0, 1]
        self.assertEqual(len(counts), 10)
        # 0.0 in [0, 0.1); 0.05 in [0, 0.1); 0.5 in [0.5, 0.6);
        # 0.95 in [0.9, 1.0); 1.0 in last bin (numpy includes upper edge)
        self.assertEqual(counts[0], 2)
        self.assertEqual(counts[5], 1)
        self.assertEqual(counts[9], 2)
        self.assertEqual(sum(counts), 5)

    def test_empty_returns_zeros(self):
        from syntheticgen.validate import af_histogram
        edges, counts = af_histogram([], n_bins=5)
        self.assertEqual(counts, [0, 0, 0, 0, 0])
        self.assertEqual(len(edges), 6)


@unittest.skipUnless(HAS_NUMPY, "numpy required")
class TestR2Pair(unittest.TestCase):
    def test_perfect_correlation(self):
        from syntheticgen.validate import _r2_pair
        a = [0, 1, 2, 0, 1, 2, 0, 1]
        b = [0, 1, 2, 0, 1, 2, 0, 1]
        self.assertAlmostEqual(_r2_pair(a, b), 1.0)

    def test_perfect_anticorrelation(self):
        from syntheticgen.validate import _r2_pair
        a = [0, 1, 2, 0, 1, 2, 0, 1]
        b = [2, 1, 0, 2, 1, 0, 2, 1]
        # r = -1 → r² = 1
        self.assertAlmostEqual(_r2_pair(a, b), 1.0)

    def test_uncorrelated(self):
        from syntheticgen.validate import _r2_pair
        # Build a clear independent pair
        a = [0, 0, 1, 1, 2, 2, 0, 1]
        b = [1, 2, 0, 2, 1, 0, 2, 1]
        # Just check the value is finite and within [0, 1]
        r2 = _r2_pair(a, b)
        self.assertFalse(math.isnan(r2))
        self.assertGreaterEqual(r2, 0.0)
        self.assertLessEqual(r2, 1.0)

    def test_constant_vector_returns_nan(self):
        from syntheticgen.validate import _r2_pair
        a = [1] * 8
        b = [0, 1, 2, 0, 1, 2, 0, 1]
        self.assertTrue(math.isnan(_r2_pair(a, b)))

    def test_too_few_samples_returns_nan(self):
        from syntheticgen.validate import _r2_pair
        a = [0, 1, 2]
        b = [0, 1, 2]
        self.assertTrue(math.isnan(_r2_pair(a, b)))

    def test_handles_missing(self):
        from syntheticgen.validate import _r2_pair
        a = [0, 1, 2, -1, 1, 2, 0, 1]
        b = [0, 1, 2, 0, 1, 2, 0, 1]
        # One pair masked out; remaining pairs are perfectly correlated.
        self.assertAlmostEqual(_r2_pair(a, b), 1.0)


@unittest.skipUnless(HAS_NUMPY, "numpy required")
class TestLdDecay(unittest.TestCase):
    def test_returns_one_dict_per_bin(self):
        from syntheticgen.validate import ld_decay

        # 4 samples, 6 SNPs at increasing positions on chr22 — set up
        # so the decay direction is clear: short-range pairs are
        # perfectly correlated, long-range pairs anticorrelated.
        # Three perfectly-correlated rows (cols 0..2 share a column),
        # then col 3 is anticorrelated.
        mat = np.array([
            [0, 0, 0, 2, 0, 1],
            [1, 1, 1, 1, 2, 0],
            [2, 2, 2, 0, 1, 2],
            [0, 0, 0, 2, 2, 0],
        ], dtype=np.int8)
        positions = [100, 600, 1100, 5_000_000, 6_000_000, 6_500_000]
        chroms = ["22"] * 6
        bins = ld_decay(mat, positions, chroms,
                        distance_bins_kb=((0.1, 1.0), (1.0, 10.0)),
                        pairs_per_bin=100)
        self.assertEqual(len(bins), 2)
        self.assertEqual({"low_kb", "high_kb", "n_pairs", "mean_r2"},
                         set(bins[0]))

    def test_short_range_higher_than_long(self):
        from syntheticgen.validate import ld_decay

        # Three "blocks" of perfectly-correlated SNPs at short range,
        # then far apart from each other — short-range r² should be
        # ~1, long-range much lower.
        rng_data = np.array([
            [0, 0, 0,  2, 2, 2,  0, 0, 0],
            [1, 1, 1,  0, 0, 0,  2, 2, 2],
            [2, 2, 2,  1, 1, 1,  1, 1, 1],
            [0, 0, 0,  2, 2, 2,  0, 0, 0],
            [1, 1, 1,  1, 1, 1,  2, 2, 2],
        ], dtype=np.int8)
        positions = [100, 200, 300,
                     1_000_000, 1_000_100, 1_000_200,
                     5_000_000, 5_000_100, 5_000_200]
        chroms = ["22"] * 9
        bins = ld_decay(rng_data, positions, chroms,
                        distance_bins_kb=((0.05, 1.0), (500.0, 10_000.0)))
        short = bins[0]["mean_r2"]
        long_ = bins[1]["mean_r2"]
        self.assertGreater(short, 0.95)  # within-block ≈ 1
        # Cross-block r² should be lower (often much lower)
        self.assertLess(long_, short)


@unittest.skipUnless(HAS_NUMPY and HAS_SKLEARN, "numpy + sklearn required")
class TestCohortPca(unittest.TestCase):
    def test_basic_pca_runs(self):
        from syntheticgen.validate import cohort_pca

        # Two clear clusters in dosage space: rows 0..3 hom-ref on
        # cols 0..4 and hom-alt on cols 5..9; rows 4..7 vice versa.
        mat = np.zeros((8, 10), dtype=np.int8)
        mat[:4, :5] = 0
        mat[:4, 5:] = 2
        mat[4:, :5] = 2
        mat[4:, 5:] = 0
        transformed, evr, kept = cohort_pca(mat, n_components=2)
        self.assertEqual(transformed.shape, (8, 2))
        self.assertEqual(len(evr), 2)
        # PC1 should pick up almost all the variance
        self.assertGreater(evr[0], 0.95)

    def test_too_few_columns_returns_none(self):
        from syntheticgen.validate import cohort_pca
        mat = np.zeros((4, 1), dtype=np.int8)
        mat[:, 0] = 1  # constant → pruned → 0 columns left
        out = cohort_pca(mat, n_components=2)
        self.assertEqual(out, (None, None, None))


@unittest.skipUnless(HAS_MPL, "matplotlib required")
class TestPlots(unittest.TestCase):
    """Plot helpers should write a PNG without raising. We don't
    inspect the image bytes — just that the file exists, is non-empty,
    and starts with the PNG magic header."""

    def _check_png(self, path):
        self.assertTrue(path.exists())
        self.assertGreater(path.stat().st_size, 100)
        with open(path, "rb") as fh:
            self.assertEqual(fh.read(8),
                             b"\x89PNG\r\n\x1a\n")

    def test_plot_ld_decay(self):
        from syntheticgen import plots
        import tempfile
        bins = [
            {"low_kb": 0.1, "high_kb": 1.0, "n_pairs": 10, "mean_r2": 0.6},
            {"low_kb": 1.0, "high_kb": 10.0, "n_pairs": 5, "mean_r2": 0.2},
        ]
        with tempfile.TemporaryDirectory() as d:
            out = plots.plot_ld_decay(bins, Path(d) / "ld.png")
            self._check_png(out)

    def test_plot_af_histogram(self):
        from syntheticgen import plots
        import tempfile
        edges = [0.0, 0.5, 1.0]
        counts = [10, 4]
        with tempfile.TemporaryDirectory() as d:
            out = plots.plot_af_histogram(
                edges, counts, Path(d) / "af.png")
            self._check_png(out)

    def test_plot_indel_lengths(self):
        from syntheticgen import plots
        import tempfile
        h = {-2: 5, -1: 2, 1: 3, 5: 1}
        with tempfile.TemporaryDirectory() as d:
            out = plots.plot_indel_lengths(h, Path(d) / "indel.png")
            self._check_png(out)

    def test_plot_pca_handles_none(self):
        from syntheticgen import plots
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            out = plots.plot_pca(None, [], Path(d) / "pca.png")
            self._check_png(out)


class TestCohortPerRegionDensity(unittest.TestCase):
    """Tier 2 #6: per-Mb variant density bins.

    Flat under today's uniform-μ simulation; should show gene-density
    structure post-M14. The validator surfaces both the per-chrom
    table and a coefficient-of-variation diagnostic.
    """

    def _stats(self, name: str, density: dict) -> SampleStats:
        s = SampleStats(name=name)
        for chrom, bins in density.items():
            for bin_idx, count in bins.items():
                s.density_bin_counts[chrom][bin_idx] = count
        return s

    def test_aggregates_counts_per_bin_across_samples(self):
        from syntheticgen.validate import cohort_per_region_density
        a = self._stats("a", {"22": {0: 5, 1: 3}})
        b = self._stats("b", {"22": {0: 7, 2: 2}, "1": {0: 10}})
        out = cohort_per_region_density([a, b])
        # Sorted by canonical chromosome order (1 before 22).
        self.assertEqual(list(out.keys()), ["1", "22"])
        # Bins sorted within each chrom; counts summed across samples.
        self.assertEqual(out["22"], [
            {"start_mb": 0, "end_mb": 1, "count": 12},
            {"start_mb": 1, "end_mb": 2, "count": 3},
            {"start_mb": 2, "end_mb": 3, "count": 2},
        ])
        self.assertEqual(out["1"], [
            {"start_mb": 0, "end_mb": 1, "count": 10},
        ])

    def test_empty_density_returns_empty_dict(self):
        from syntheticgen.validate import cohort_per_region_density
        self.assertEqual(cohort_per_region_density([]), {})

    def test_canonical_chrom_order(self):
        from syntheticgen.validate import cohort_per_region_density
        s = self._stats("a", {
            "Y": {0: 1}, "X": {0: 1}, "2": {0: 1}, "22": {0: 1},
        })
        out = cohort_per_region_density([s])
        self.assertEqual(list(out.keys()), ["2", "22", "X", "Y"])


class TestCohortQualityMetrics(unittest.TestCase):
    """Tier 2 #7: DP / GQ / AD-ref-fraction distribution sanity.

    Surfaces drift in the empirical targets baked into quality.py.
    Targets are imported from ``quality.py`` (``DEFAULT_DP_MEAN``,
    ``1 - HET_ALT_FRAC``) so a future re-calibration of those
    constants automatically updates the validator's expectations.
    """

    def _stats(self, name: str, dp: list, gq: list, ad: list) -> SampleStats:
        s = SampleStats(name=name)
        s.dp_samples = list(dp)
        s.gq_samples = list(gq)
        s.ad_het_ref_frac_samples = list(ad)
        return s

    def test_summary_stats_attached_per_metric(self):
        from syntheticgen.validate import cohort_quality_metrics
        from syntheticgen.quality import (
            DEFAULT_DP_MEAN, HET_ALT_FRAC,
        )
        a = self._stats("a", dp=[28, 30, 32], gq=[99, 99, 95],
                        ad=[0.48, 0.50, 0.45])
        out = cohort_quality_metrics([a])
        # Targets come directly from quality.py — no magic numbers.
        self.assertEqual(out["dp"]["target"], DEFAULT_DP_MEAN)
        self.assertAlmostEqual(out["dp"]["mean"], 30.0)
        self.assertEqual(out["dp"]["n"], 3)
        # AD ref-fraction target is ``1 - HET_ALT_FRAC`` because
        # ``HET_ALT_FRAC`` is the ALT share at hets (under reference
        # bias) and we record the REF share. Off-by-direction in
        # the original would have given 0.475 instead of 0.525.
        self.assertAlmostEqual(
            out["ad_het_ref_fraction"]["target"],
            1.0 - HET_ALT_FRAC,
        )
        self.assertAlmostEqual(
            out["ad_het_ref_fraction"]["mean"], 0.4767, places=4,
        )

    def test_empty_samples_returns_none_stats(self):
        from syntheticgen.validate import cohort_quality_metrics
        out = cohort_quality_metrics([])
        self.assertEqual(out["dp"]["n"], 0)
        self.assertIsNone(out["dp"]["mean"])

    def test_percentiles_are_true_quantiles_not_index_picks(self):
        # Regression for PR #80 review #1: the original
        # ``_summarise`` computed ``p90 = values_sorted[int(0.90 *
        # n)]`` which for n=10 returned the max (index 9). True p90
        # via linear interpolation on [1..10] is 9.1. Locking in
        # that the new implementation uses numpy percentile rather
        # than an index pick.
        from syntheticgen.validate import cohort_quality_metrics
        # 10 values: 1, 2, ..., 10. Old code would have given
        # p90 = 10 (the max); true p90 ≈ 9.1.
        a = self._stats(
            "a", dp=list(range(1, 11)), gq=[], ad=[],
        )
        out = cohort_quality_metrics([a])
        self.assertEqual(out["dp"]["n"], 10)
        # numpy's default linear interpolation gives 9.1 for p90
        # of [1..10]. Old index pick would have returned 10.
        self.assertAlmostEqual(out["dp"]["p90"], 9.1, places=3)
        # p10 of [1..10] ≈ 1.9; old code gave 2.
        self.assertAlmostEqual(out["dp"]["p10"], 1.9, places=3)
        # Median of [1..10] is 5.5; old ``values_sorted[n // 2]``
        # picked index 5 → 6 (the upper-middle).
        self.assertAlmostEqual(out["dp"]["median"], 5.5)

    def test_aggregates_across_samples(self):
        from syntheticgen.validate import cohort_quality_metrics
        a = self._stats("a", dp=[30], gq=[99], ad=[0.5])
        b = self._stats("b", dp=[30, 30], gq=[99, 99], ad=[0.4, 0.5])
        out = cohort_quality_metrics([a, b])
        self.assertEqual(out["dp"]["n"], 3)
        self.assertEqual(out["ad_het_ref_fraction"]["n"], 3)


class TestCohortFStatistic(unittest.TestCase):
    """Tier 2 #8: per-sample inbreeding coefficient F."""

    def test_balanced_hwe_matrix_gives_per_sample_f_zero(self):
        # 4 samples × 4 variants. Each column sums to 4 so cohort
        # AF = 0.5 at every variant (expected_het per variant =
        # 2·0.5·0.5 = 0.5; expected_het per sample = 4 × 0.5 = 2).
        # Each row has exactly 2 hets (dosage=1), so observed_het
        # per sample = 2. F = 1 − 2/2 = 0 for every sample.
        # This is the property the name asserts.
        from syntheticgen.validate import cohort_f_statistic
        import numpy as np
        matrix = np.array([
            [1, 1, 0, 2],
            [1, 0, 2, 1],
            [0, 2, 1, 1],
            [2, 1, 1, 0],
        ])
        out = cohort_f_statistic(matrix)
        self.assertEqual(len(out["per_sample"]), 4)
        for entry in out["per_sample"]:
            self.assertAlmostEqual(entry["f"], 0.0, places=10)
            self.assertEqual(entry["observed_het"], 2)
            self.assertAlmostEqual(entry["expected_het"], 2.0)
        self.assertAlmostEqual(out["cohort_mean"], 0.0)
        self.assertAlmostEqual(out["cohort_median"], 0.0)

    def test_homozygous_extreme_gives_f_near_one(self):
        # All hom-ref + all hom-alt samples at variants where the
        # cohort has intermediate AF. Each sample's observed_het is
        # 0; expected_het > 0; F = 1.0. Documents the extreme-end
        # of the F-statistic interpretation table.
        from syntheticgen.validate import cohort_f_statistic
        import numpy as np
        matrix = np.array([
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [2, 2, 2, 2],
            [2, 2, 2, 2],
        ])
        out = cohort_f_statistic(matrix)
        for entry in out["per_sample"]:
            self.assertEqual(entry["observed_het"], 0)
            self.assertAlmostEqual(entry["f"], 1.0)
        self.assertAlmostEqual(out["cohort_mean"], 1.0)

    def test_empty_matrix_returns_empty_result(self):
        from syntheticgen.validate import cohort_f_statistic
        import numpy as np
        out = cohort_f_statistic(np.zeros((0, 0), dtype=int))
        self.assertEqual(out["per_sample"], [])
        self.assertIsNone(out["cohort_mean"])

    def test_per_sample_struct_shape(self):
        from syntheticgen.validate import cohort_f_statistic
        import numpy as np
        # Simple 2x2 with one het each side: AF = 0.5 at each variant.
        matrix = np.array([[1, 0], [0, 1]])
        out = cohort_f_statistic(matrix)
        for entry in out["per_sample"]:
            self.assertIn("sample_idx", entry)
            self.assertIn("observed_het", entry)
            self.assertIn("expected_het", entry)
            self.assertIn("f", entry)


class TestCohortAncestryTracts(unittest.TestCase):
    """Tier 2 #9: tract-length distribution from per-person
    ancestry BEDs (admixture mode)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def _write_bed(self, name: str, rows: list) -> Path:
        path = self.dir / name
        path.write_text(
            "".join(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\n"
                    for r in rows),
        )
        return path

    def test_single_population_single_tract(self):
        from syntheticgen.validate import cohort_ancestry_tracts
        # One person, one row: hap1 and hap2 both EUR.
        p = self._write_bed("p.bed", [
            ("22", 0, 1_000_000, "EUR", "EUR"),
        ])
        out = cohort_ancestry_tracts([p])
        # Two tracts (one per haplotype), each 1 Mb of EUR.
        self.assertEqual(out["by_population"]["EUR"]["n"], 2)
        self.assertEqual(
            out["by_population"]["EUR"]["mean_bp"], 1_000_000,
        )

    def test_consecutive_same_pop_rows_merge_into_one_tract(self):
        from syntheticgen.validate import cohort_ancestry_tracts
        # Adjacent rows with the same hap1_pop but different hap2_pop
        # — hap1's tract spans both rows, hap2 has two tracts.
        p = self._write_bed("p.bed", [
            ("22", 0, 1_000_000, "EUR", "EUR"),
            ("22", 1_000_000, 2_000_000, "EUR", "SAS"),
        ])
        out = cohort_ancestry_tracts([p])
        # hap1: one EUR tract of 2 Mb (rows merged).
        # hap2: one EUR tract of 1 Mb + one SAS tract of 1 Mb.
        eur = out["by_population"]["EUR"]
        sas = out["by_population"]["SAS"]
        # Two EUR tracts total: hap1's 2Mb + hap2's 1Mb.
        self.assertEqual(eur["n"], 2)
        self.assertEqual(sas["n"], 1)
        self.assertEqual(sas["mean_bp"], 1_000_000)

    def test_chrom_boundary_breaks_tract(self):
        from syntheticgen.validate import cohort_ancestry_tracts
        # Same hap1_pop across two chroms — must be reported as two
        # tracts, not one (a chr1 and chr22 stretch are not
        # contiguous even if the ancestry happens to match).
        p = self._write_bed("p.bed", [
            ("1", 0, 1_000_000, "EUR", "EUR"),
            ("22", 0, 1_000_000, "EUR", "EUR"),
        ])
        out = cohort_ancestry_tracts([p])
        # 2 EUR tracts × 2 haplotypes = 4 tracts.
        self.assertEqual(out["by_population"]["EUR"]["n"], 4)

    def test_missing_file_skipped_silently(self):
        # An OSError on the BED read shouldn't tank the whole run.
        from syntheticgen.validate import cohort_ancestry_tracts
        out = cohort_ancestry_tracts([Path("/nonexistent.bed")])
        self.assertEqual(out["by_population"], {})


class TestSummariseVcfTier2(unittest.TestCase):
    """End-to-end Tier 2 wiring through ``summarise_vcf`` —
    confirms the new SampleStats fields are populated correctly
    when records flow through the real (mocked-iter_records)
    summarisation path."""

    def _record(self, chrom: str, pos: int, gt: str,
                dp: int = 30, gq: int = 99,
                ad_ref: int = 15, ad_alt: int = 15) -> Record:
        return Record(
            chrom=chrom, pos=pos, ref="A", alt="C", gt=gt,
            dp=dp, gq=gq, ad_ref=ad_ref, ad_alt=ad_alt,
            info={},
        )

    def _summarise(self, records: list):
        from syntheticgen.validate import summarise_vcf
        with patch(
            "syntheticgen.validate.iter_records",
            return_value=iter(records),
        ):
            return summarise_vcf(Path("/nonexistent.vcf.gz"), name="t")

    def test_density_bins_increment_per_pos(self):
        # Three records on chr22 at positions 100, 1_500_000,
        # 1_999_999. Bin 0 gets 1, bin 1 gets 2.
        stats = self._summarise([
            self._record("22", 100, "0|1"),
            self._record("22", 1_500_000, "0|1"),
            self._record("22", 1_999_999, "0|1"),
        ])
        self.assertEqual(stats.density_bin_counts["22"][0], 1)
        self.assertEqual(stats.density_bin_counts["22"][1], 2)

    def test_dp_gq_samples_captured(self):
        stats = self._summarise([
            self._record("22", 100, "0|1", dp=28, gq=99),
            self._record("22", 200, "0|1", dp=32, gq=99),
        ])
        self.assertEqual(stats.dp_samples, [28, 32])
        self.assertEqual(stats.gq_samples, [99, 99])

    def test_ad_ref_fraction_recorded_only_at_hets(self):
        # Two hets (dosage=1) and one hom-alt (dosage=2) — the
        # AD ref-fraction is only recorded at hets.
        stats = self._summarise([
            self._record("22", 100, "0|1", ad_ref=14, ad_alt=16),
            self._record("22", 200, "1|1", ad_ref=0, ad_alt=30),
            self._record("22", 300, "1|0", ad_ref=15, ad_alt=15),
        ])
        self.assertEqual(len(stats.ad_het_ref_frac_samples), 2)
        # First het: 14/30 ≈ 0.467.
        self.assertAlmostEqual(
            stats.ad_het_ref_frac_samples[0], 14 / 30, places=4,
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
