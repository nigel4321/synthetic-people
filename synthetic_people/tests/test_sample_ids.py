"""Tests for sample-ID generation in ``syntheticgen.background``.

Regression-driven: the user hit
``[E::bcf_hdr_add_sample_len] Duplicated sample name 'NA44856'`` at
``--n 3000`` because the legacy single-call ``random_sample_id``
drew from a 180k-name pool (HG/NA × 5-digit number), and at n≈600+
the birthday paradox routinely produced duplicates. ``bcftools``
rejects multi-sample BCFs with duplicate sample columns, so the
cohort-write path needs guaranteed uniqueness.

The new :func:`draw_sample_ids` batch-draws via
``rng.sample(range(pool_size), n)``, which picks integer keys
without replacement and decodes each to a unique
``(prefix, number)`` pair.
"""

from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.background import (
    _SAMPLE_ID_POOL_SIZE,
    draw_sample_ids,
    random_sample_id,
)


class DrawSampleIdsUniquenessTest(unittest.TestCase):
    """The headline contract: no duplicates across any n we'd
    realistically run."""

    def test_3000_samples_no_duplicates(self):
        # The user-reported failure point.
        ids = draw_sample_ids(3000, random.Random(42))
        self.assertEqual(len(ids), 3000)
        self.assertEqual(len(set(ids)), 3000)

    def test_100k_samples_no_duplicates(self):
        # Phase 5 stretch target. rng.sample on a ~20M-size range
        # is a constant-time-per-element shuffle; this should run in
        # a fraction of a second even though n is large.
        ids = draw_sample_ids(100_000, random.Random(0))
        self.assertEqual(len(ids), 100_000)
        self.assertEqual(len(set(ids)), 100_000)

    def test_at_n_600_legacy_function_would_have_collided(self):
        # Sanity check that the OLD per-call ``random_sample_id``
        # would have produced duplicates at this scale — i.e. the
        # bug was real, not theoretical. Birthday paradox math says
        # at n=600 we expect ~1 collision; using a fresh deterministic
        # seed, we should be able to reproduce that.
        rng = random.Random(2023)
        seen: set = set()
        collisions = 0
        for _ in range(600):
            sid = random_sample_id(rng)
            if sid in seen:
                collisions += 1
            seen.add(sid)
        self.assertGreater(collisions, 0,
                           "expected at least one duplicate from the "
                           "legacy random_sample_id at n=600 — if this "
                           "passes the bug isn't reproducing, check the "
                           "pool size constants")


class DrawSampleIdsShapeTest(unittest.TestCase):
    """Each emitted ID is HG/NA-prefixed and has the expected
    digit width."""

    def test_each_id_has_valid_prefix_and_number(self):
        ids = draw_sample_ids(50, random.Random(1))
        for sid in ids:
            self.assertIn(sid[:2], ("HG", "NA"))
            self.assertTrue(sid[2:].isdigit())
            number = int(sid[2:])
            self.assertGreaterEqual(number, 100_000)
            self.assertLess(number, 10_000_000)

    def test_zero_n_returns_empty_list(self):
        self.assertEqual(draw_sample_ids(0, random.Random(0)), [])

    def test_negative_n_returns_empty_list(self):
        # Defensive — caller bug shouldn't blow up here.
        self.assertEqual(draw_sample_ids(-5, random.Random(0)), [])


class DrawSampleIdsDeterminismTest(unittest.TestCase):
    """Same seed + same n → same list. The cohort runs depend on this
    for resume + reproducibility."""

    def test_same_seed_same_output(self):
        a = draw_sample_ids(100, random.Random(42))
        b = draw_sample_ids(100, random.Random(42))
        self.assertEqual(a, b)

    def test_different_seeds_different_output(self):
        a = draw_sample_ids(100, random.Random(1))
        b = draw_sample_ids(100, random.Random(2))
        self.assertNotEqual(a, b)


class DrawSampleIdsLimitsTest(unittest.TestCase):
    """Above the pool size the function refuses rather than
    silently risking a duplicate."""

    def test_pool_overflow_raises(self):
        # Asking for more IDs than the pool can offer should raise,
        # not loop forever or silently return a list with duplicates.
        with self.assertRaisesRegex(ValueError, "pool"):
            draw_sample_ids(_SAMPLE_ID_POOL_SIZE + 1, random.Random(0))

    def test_at_pool_size_works(self):
        # Pulling exactly the pool size is rng.sample-equivalent; it
        # should still succeed and return that many distinct IDs.
        # Skipped if the pool is larger than 2M to keep the test
        # fast — the contract is the same.
        if _SAMPLE_ID_POOL_SIZE > 2_000_000:
            self.skipTest("pool too large to materialise in a unit test")
        ids = draw_sample_ids(_SAMPLE_ID_POOL_SIZE, random.Random(0))
        self.assertEqual(len(set(ids)), _SAMPLE_ID_POOL_SIZE)


if __name__ == "__main__":
    unittest.main()
