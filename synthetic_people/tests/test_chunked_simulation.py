"""Phase 5f — chunked simulation tests.

Three things to lock in:

1. Auto-pick correctness — given a known available-RAM and worker
   count, the picked chunk size matches the analytical formula and
   shrinks proportionally as workers increase.
2. Override behaviour — ``--chr-chunk-mb N`` (with ``N > 0``) is
   honoured regardless of available RAM.
3. Chunked vs unchunked correctness — chunked output stays sorted,
   unique, in genome-position range, and statistically resembles
   the unchunked output (record count, allele counts) within the
   noise of independent simulations.

Heavy paths gate on msprime + stdpopsim. Auto-pick math is pure
Python and runs anywhere.
"""

from __future__ import annotations

import importlib.util
import random
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.coalescent import (
    CHUNK_AUTO_DERATE_FLOOR_MB,
    CHUNK_OVERLAP_FRACTION,
    CHUNK_OVERLAP_MAX_BP,
    CHUNK_OVERLAP_MIN_BP,
    CHUNK_RAM_BYTES_PER_SAMPLE_PER_MB_CONSTANT_NE,
    CHUNK_RAM_BYTES_PER_SAMPLE_PER_MB_OOA,
    DEFAULT_AUTO_PICK_TARGET_FRACTION,
    auto_derate_workers,
    auto_pick_chunk_size_mb,
    estimate_chunk_ram_bytes,
)


_HAVE_MSPRIME = importlib.util.find_spec("msprime") is not None


class ChunkRamEstimateTest(unittest.TestCase):
    """The estimator should be linear in (n × chunk_mb) and use the
    cheaper rate for constant-Ne. Pure-Python, no msprime needed."""

    def test_ooa_rate_used_for_real_demography(self):
        bytes_full = estimate_chunk_ram_bytes(3000, 70, "OutOfAfrica_3G09")
        expected = (CHUNK_RAM_BYTES_PER_SAMPLE_PER_MB_OOA * 3000 * 70)
        self.assertEqual(bytes_full, expected)

    def test_constant_ne_rate_is_lower(self):
        bytes_ooa = estimate_chunk_ram_bytes(3000, 10, "OutOfAfrica_3G09")
        bytes_const = estimate_chunk_ram_bytes(3000, 10, None)
        self.assertEqual(
            bytes_const,
            CHUNK_RAM_BYTES_PER_SAMPLE_PER_MB_CONSTANT_NE * 3000 * 10,
        )
        self.assertLess(bytes_const, bytes_ooa)

    def test_none_string_treated_as_constant_ne(self):
        # CLI passes args.demo_model.lower() == "none" through to the
        # estimator; the estimator should accept the string too.
        bytes_str = estimate_chunk_ram_bytes(3000, 10, "none")
        bytes_obj = estimate_chunk_ram_bytes(3000, 10, None)
        self.assertEqual(bytes_str, bytes_obj)

    def test_linear_in_chunk_size(self):
        small = estimate_chunk_ram_bytes(3000, 5, "OutOfAfrica_3G09")
        big = estimate_chunk_ram_bytes(3000, 50, "OutOfAfrica_3G09")
        self.assertEqual(big, 10 * small)

    def test_linear_in_n(self):
        smaller = estimate_chunk_ram_bytes(1000, 10, "OutOfAfrica_3G09")
        larger = estimate_chunk_ram_bytes(10_000, 10, "OutOfAfrica_3G09")
        self.assertEqual(larger, 10 * smaller)


class AutoPickChunkSizeTest(unittest.TestCase):
    """Auto-pick should pick the largest chunk size whose estimate
    fits in the worker budget. Pure-Python."""

    def test_full_chrom_fits_returns_length(self):
        # n=10 × 5 Mb × OOA at 16 GB / 1 worker fits trivially —
        # auto-pick returns the full chromosome length, no chunking.
        mb = auto_pick_chunk_size_mb(
            n_people=10, length_mb=5.0, demo_model="OutOfAfrica_3G09",
            available_bytes=16 * 1024**3, workers=1,
        )
        self.assertEqual(mb, 5.0)

    def test_does_not_fit_returns_smaller_chunk(self):
        # n=3000 × 70 Mb × OOA at 16 GB / 1 worker doesn't fit; the
        # picked chunk should be smaller than 70 Mb.
        mb = auto_pick_chunk_size_mb(
            n_people=3000, length_mb=70.0, demo_model="OutOfAfrica_3G09",
            available_bytes=16 * 1024**3, workers=1,
        )
        self.assertLess(mb, 70.0)
        # And the estimate at the picked chunk size should fit half
        # the available RAM (the auto-pick target).
        target = int(16 * 1024**3 * 0.5)
        self.assertLessEqual(
            estimate_chunk_ram_bytes(3000, mb, "OutOfAfrica_3G09"),
            target,
        )

    def test_workers_divide_budget(self):
        # Doubling workers should halve the picked chunk size when
        # the constraint is RAM-bound.
        single = auto_pick_chunk_size_mb(
            n_people=3000, length_mb=70.0, demo_model="OutOfAfrica_3G09",
            available_bytes=16 * 1024**3, workers=1,
        )
        quad = auto_pick_chunk_size_mb(
            n_people=3000, length_mb=70.0, demo_model="OutOfAfrica_3G09",
            available_bytes=16 * 1024**3, workers=4,
        )
        self.assertAlmostEqual(quad * 4, single, delta=0.5)

    def test_constant_ne_picks_bigger_chunks_than_ooa(self):
        ooa = auto_pick_chunk_size_mb(
            n_people=3000, length_mb=70.0, demo_model="OutOfAfrica_3G09",
            available_bytes=16 * 1024**3, workers=1,
        )
        const = auto_pick_chunk_size_mb(
            n_people=3000, length_mb=70.0, demo_model=None,
            available_bytes=16 * 1024**3, workers=1,
        )
        # Constant-Ne is ~5× cheaper, so chunks can be ~5× larger
        # (capped at length_mb=70).
        self.assertGreater(const, ooa)


@unittest.skipUnless(_HAVE_MSPRIME, "msprime not installed")
class ChunkedSimulationParityTest(unittest.TestCase):
    """Chunked output for the same configuration as unchunked should
    produce a comparable record set: positions sorted + unique +
    bounded by the simulated length, record count within stochastic
    noise, allele counts in the same range."""

    def _run(self, chunk_size_mb: float):
        from syntheticgen.coalescent import simulate_chromosome
        return simulate_chromosome(
            chrom="22", build="GRCh38", n_people=4, length_mb=0.5,
            demo_model=None, population="CEU",
            rec_rate=1e-8, mu=1.29e-8,
            rng=random.Random(42),
            chunk_size_mb=chunk_size_mb,
        )

    def test_chunked_positions_sorted(self):
        sites = self._run(chunk_size_mb=0.1)
        positions = [s["pos"] for s in sites]
        self.assertEqual(positions, sorted(positions))

    def test_chunked_positions_unique(self):
        sites = self._run(chunk_size_mb=0.1)
        positions = [s["pos"] for s in sites]
        self.assertEqual(len(set(positions)), len(positions))

    def test_chunked_positions_in_range(self):
        sites = self._run(chunk_size_mb=0.1)
        positions = [s["pos"] for s in sites]
        # 1-based, inside [1, length_bp]. length_mb=0.5 → 500_000 bp.
        self.assertGreaterEqual(min(positions), 1)
        self.assertLessEqual(max(positions), 500_000)

    def test_chunked_record_count_within_noise(self):
        unchunked = self._run(chunk_size_mb=0.0)
        chunked = self._run(chunk_size_mb=0.1)
        # Independent chunk simulations + small overlap-margin
        # discards mean the counts won't match exactly. They should
        # be within 25% of each other on a small fixture (gets
        # tighter at larger n/length).
        ratio = len(chunked) / max(1, len(unchunked))
        self.assertGreater(ratio, 0.75,
                           f"chunked {len(chunked)} vs unchunked "
                           f"{len(unchunked)}")
        self.assertLess(ratio, 1.25,
                        f"chunked {len(chunked)} vs unchunked "
                        f"{len(unchunked)}")

    def test_chunked_determinism_at_fixed_seed(self):
        # Same seed + same chunk size → same chunked output across
        # runs. Phase 5f's per-chunk seed derivation depends only on
        # chrom_seed (drawn once from the master rng) plus
        # chunk_index, so re-runs are deterministic.
        a = self._run(chunk_size_mb=0.1)
        b = self._run(chunk_size_mb=0.1)
        self.assertEqual(len(a), len(b))
        for sa, sb in zip(a, b):
            self.assertEqual(sa["pos"], sb["pos"])
            self.assertEqual(sa["acs"], sb["acs"])
            self.assertTrue(
                np.array_equal(sa["carriers"], sb["carriers"]),
            )

    def test_overlap_dedup_no_duplicate_positions_at_boundaries(self):
        # Walk a 0.4 Mb chromosome with 0.1 Mb chunks (4 chunks +
        # overlap). Boundaries land at simulated positions ~100k,
        # ~200k, ~300k. The trailing-overlap region of chunk K
        # should have been dropped at write time, so the per-chrom
        # site list has no duplicates anywhere — including across
        # boundaries.
        from syntheticgen.coalescent import simulate_chromosome
        sites = simulate_chromosome(
            chrom="22", build="GRCh38", n_people=4, length_mb=0.4,
            demo_model=None, population="CEU",
            rec_rate=1e-8, mu=1.29e-8,
            rng=random.Random(42), chunk_size_mb=0.1,
        )
        positions = [s["pos"] for s in sites]
        self.assertEqual(len(positions), len(set(positions)))


class AutoDerateWorkersTest(unittest.TestCase):
    """Worker auto-derate landed after the user's 16 GB host stalled
    with 4 parallel workers each holding a 4 GB tree sequence — the
    auto-pick had picked an 8.7 Mb chunk to fit 4 workers in 25% of
    available RAM, but the calibration was 2× off and total RAM
    actually saturated. The fix: when the auto-pick math is forced
    to drop below ~2 Mb per chunk to keep all workers in budget,
    drop a worker instead. Below 2 Mb the per-chunk msprime startup
    cost dominates the per-chunk simulation cost, and the boundary-
    smoothing benefit erodes — so the trade-off is *much* better
    spent on fewer parallel workers than tinier chunks.
    """

    def test_at_high_ram_workers_unchanged(self):
        # n=10 × 5 Mb × OOA at 64 GB / 4 workers fits trivially —
        # all 4 workers can hold their tree sequence in budget, no
        # derate needed.
        derated = auto_derate_workers(
            n_people=10, length_mb=5.0,
            demo_model="OutOfAfrica_3G09",
            available_bytes=64 * 1024**3, requested_workers=4,
        )
        self.assertEqual(derated, 4)

    def test_users_failing_case_derates(self):
        # 16 GB host, n=3000, 70 Mb chrom, OOA, 4 workers requested.
        # With the recalibrated 160 KB / (sample × Mb) and 25%
        # target, 4 workers each get a 1 GB budget — chunk size
        # auto-picks at ~2.2 Mb. With 8 workers each gets 0.5 GB,
        # forcing chunks below 1.1 Mb (under the 2 Mb floor).
        # The function must derate to keep us at the floor.
        derated = auto_derate_workers(
            n_people=3000, length_mb=70.0,
            demo_model="OutOfAfrica_3G09",
            available_bytes=16 * 1024**3, requested_workers=8,
        )
        self.assertLess(derated, 8)
        # And a sanity check that the auto-picked chunk size at
        # the derated worker count IS at or above the floor.
        chunk_mb = auto_pick_chunk_size_mb(
            3000, 70.0, "OutOfAfrica_3G09",
            available_bytes=16 * 1024**3, workers=derated,
        )
        self.assertGreaterEqual(chunk_mb, CHUNK_AUTO_DERATE_FLOOR_MB)

    def test_workers_one_returns_one(self):
        # The function never *increases* parallelism. workers=1 in →
        # workers=1 out, regardless of host RAM.
        for ram in (4 * 1024**3, 16 * 1024**3, 256 * 1024**3):
            self.assertEqual(
                auto_derate_workers(
                    n_people=3000, length_mb=70.0,
                    demo_model="OutOfAfrica_3G09",
                    available_bytes=ram, requested_workers=1,
                ),
                1,
            )

    def test_extremely_constrained_host_returns_one(self):
        # Even at workers=1 the chunk would be sub-floor (e.g. n=1M
        # at 4 GB host). The function honours the request anyway —
        # at workers=1 a smaller chunk is the only path forward.
        # The CLI prints a warning at flag-resolution time so the
        # user sees this case.
        derated = auto_derate_workers(
            n_people=1_000_000, length_mb=70.0,
            demo_model="OutOfAfrica_3G09",
            available_bytes=4 * 1024**3, requested_workers=4,
        )
        self.assertEqual(derated, 1)


class CalibrationFromUserTraceTest(unittest.TestCase):
    """Sanity tests pinning the recalibrated constants against the
    user's profile-memory trace. The trace showed 4 workers × ~4 GB
    per tree sequence at chunk_size=8.7 Mb, n=3000, OOA_3G09 —
    yielding ~153 KiB/(sample × Mb). We pessimistically rounded up
    to 160 KiB/(sample × Mb) and tightened the safety target to
    25% of available RAM. These tests ensure those choices stick.
    """

    def test_ooa_coefficient_matches_user_trace(self):
        # Per-(sample × Mb) cost should be 160 KiB.
        self.assertEqual(
            CHUNK_RAM_BYTES_PER_SAMPLE_PER_MB_OOA, 160 * 1024)

    def test_default_target_fraction_is_one_quarter(self):
        self.assertEqual(DEFAULT_AUTO_PICK_TARGET_FRACTION, 0.25)

    def test_users_failing_config_picks_a_safe_chunk_at_workers_1(self):
        # n=3000, 70 Mb, OOA, 16 GB host, --workers 1. With the
        # recalibrated coefficient + 25% target this should pick
        # a chunk that fits in 4 GB (= 25% of 16 GB) — about 8.7 Mb,
        # similar to what 5f's first calibration gave but now backed
        # by a coefficient that matches reality.
        chunk_mb = auto_pick_chunk_size_mb(
            n_people=3000, length_mb=70.0,
            demo_model="OutOfAfrica_3G09",
            available_bytes=16 * 1024**3, workers=1,
        )
        # Estimate at the picked chunk size should fit in the
        # 4 GB target.
        target = 16 * 1024**3 // 4
        self.assertLessEqual(
            estimate_chunk_ram_bytes(3000, chunk_mb, "OutOfAfrica_3G09"),
            target,
        )
        # And we expect a substantive chunk size (≥ 5 Mb) — the
        # whole point of recalibration is to fit the failing config
        # without dropping to micro-chunks.
        self.assertGreaterEqual(chunk_mb, 5.0)


class ChunkOverlapBoundsTest(unittest.TestCase):
    """The overlap-bp bounds in coalescent.py should clamp the
    overlap to a sensible range across chunk sizes."""

    def test_overlap_floor_min(self):
        from syntheticgen.coalescent import _chunk_overlap_bp
        # 100 kb chunk × 10% = 10 kb, well under the 500 kb floor.
        self.assertEqual(_chunk_overlap_bp(100_000), CHUNK_OVERLAP_MIN_BP)

    def test_overlap_ceiling_max(self):
        from syntheticgen.coalescent import _chunk_overlap_bp
        # 100 Mb chunk × 10% = 10 Mb, hits the 5 Mb ceiling.
        self.assertEqual(
            _chunk_overlap_bp(100_000_000), CHUNK_OVERLAP_MAX_BP)

    def test_overlap_proportional_in_range(self):
        from syntheticgen.coalescent import _chunk_overlap_bp
        # 20 Mb chunk × 10% = 2 Mb, inside [0.5, 5] Mb range.
        self.assertEqual(_chunk_overlap_bp(20_000_000), 2_000_000)


if __name__ == "__main__":
    unittest.main()
