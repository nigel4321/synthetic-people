"""Tests for the carriers spill-to-disk sidecar (Fix B.1).

The sidecar's contract is straightforward — write packed carriers,
read them back at known offsets, cleanup on close — but it's the
load-bearing primitive for keeping parent RSS bounded at WGS scale
under realistic overlay densities. Locks in:

- Round-trip identity for the packed (n_rows, 2) int32 shape.
- Empty-carriers handling (the no-non-zero-haps case).
- Out-of-order reads (writes are sequential, reads are random).
- Cleanup: file unlinked on close(), even after an exception.
- Idempotent close() — safe to call twice.
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.carriers_sidecar import CarriersSidecar


class CarriersSidecarRoundTripTest(unittest.TestCase):
    """Single-block round trips for several typical shapes."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "rt.spill"

    def tearDown(self):
        self.tmp.cleanup()

    def test_typical_biallelic_carriers(self):
        with CarriersSidecar(self.path) as sc:
            c = np.array([[0, 1], [3, 1], [5, 1]], dtype=np.int32)
            offset, n_bytes = sc.write(c)
            self.assertEqual(offset, 0)
            self.assertEqual(n_bytes, 3 * 8)  # 3 rows × 8 B/row
            back = sc.read(offset, n_bytes)
            self.assertTrue(np.array_equal(back, c))

    def test_multi_allelic_preserves_allele_column(self):
        # The (n_rows, 2) shape exists precisely to keep multi-
        # allelic carriers (allele_idx > 1) intact — Fix B.1 must
        # not lose that.
        with CarriersSidecar(self.path) as sc:
            c = np.array(
                [[0, 1], [3, 2], [5, 1], [7, 2]], dtype=np.int32,
            )
            ref = sc.write(c)
            back = sc.read(*ref)
            self.assertTrue(np.array_equal(back, c))

    def test_empty_carriers_round_trip(self):
        # Sites with no non-zero haplotypes serialise as a 0-byte
        # write; read returns a fresh (0, 2) array (no allocation
        # surprises for callers iterating empties).
        with CarriersSidecar(self.path) as sc:
            empty = np.zeros((0, 2), dtype=np.int32)
            offset, n_bytes = sc.write(empty)
            self.assertEqual(n_bytes, 0)
            back = sc.read(offset, n_bytes)
            self.assertEqual(back.shape, (0, 2))
            self.assertEqual(back.dtype, np.int32)

    def test_large_block_round_trip(self):
        # n=1M analogue: ~600K carriers per common-AF site at WGS.
        # 100K rows is plenty to exercise the read/write paths
        # without bloating test runtime.
        with CarriersSidecar(self.path) as sc:
            n_rows = 100_000
            c = np.column_stack((
                np.arange(n_rows, dtype=np.int32),
                np.ones(n_rows, dtype=np.int32),
            )).astype(np.int32)
            ref = sc.write(c)
            self.assertEqual(ref[1], n_rows * 8)
            back = sc.read(*ref)
            self.assertTrue(np.array_equal(back, c))


class CarriersSidecarMultiBlockTest(unittest.TestCase):
    """Sequential writes + random-order reads. Mirrors the streaming
    pipeline: every site's carriers are appended; reads happen on
    heap-pop in position-sorted order, which is unrelated to
    insertion (tree-walk) order."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "multi.spill"

    def tearDown(self):
        self.tmp.cleanup()

    def test_sequential_write_random_read(self):
        with CarriersSidecar(self.path) as sc:
            blocks = [
                np.array([[i, 1]], dtype=np.int32)
                for i in range(10)
            ]
            refs = [sc.write(b) for b in blocks]
            # Read back in reverse order — exercises os.pread's
            # positional-read semantics independent of write head.
            for ref, expected in reversed(list(zip(refs, blocks))):
                back = sc.read(*ref)
                self.assertTrue(np.array_equal(back, expected))

    def test_interleaved_write_read(self):
        # Mid-pass reads during writes — the pattern hot-path
        # streaming uses (some sites are popped before all sites
        # are pushed). Reads must see writes that occurred earlier
        # in the same pass.
        with CarriersSidecar(self.path) as sc:
            ref_a = sc.write(np.array([[1, 1]], dtype=np.int32))
            ref_b = sc.write(np.array([[2, 2]], dtype=np.int32))
            back_a = sc.read(*ref_a)
            self.assertEqual(back_a[0, 0], 1)
            ref_c = sc.write(np.array([[3, 1]], dtype=np.int32))
            # All three reads still correct after the third write.
            self.assertEqual(sc.read(*ref_a)[0, 0], 1)
            self.assertEqual(sc.read(*ref_b)[0, 1], 2)
            self.assertEqual(sc.read(*ref_c)[0, 0], 3)


class CarriersSidecarLifecycleTest(unittest.TestCase):
    """The sidecar is transient scratch: every code path must leave
    the filesystem clean."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "lifecycle.spill"

    def tearDown(self):
        self.tmp.cleanup()

    def test_context_manager_unlinks_on_normal_exit(self):
        with CarriersSidecar(self.path) as sc:
            sc.write(np.array([[1, 1]], dtype=np.int32))
            self.assertTrue(self.path.exists())
        self.assertFalse(self.path.exists())

    def test_context_manager_unlinks_on_exception(self):
        with self.assertRaises(RuntimeError):
            with CarriersSidecar(self.path) as sc:
                sc.write(np.array([[1, 1]], dtype=np.int32))
                self.assertTrue(self.path.exists())
                raise RuntimeError("simulated stream failure")
        self.assertFalse(self.path.exists())

    def test_close_is_idempotent(self):
        sc = CarriersSidecar(self.path)
        sc.write(np.array([[1, 1]], dtype=np.int32))
        sc.close()
        # Second close must not raise even though fds + file are
        # gone — important because the cli wraps close() in a
        # try/finally that fires regardless of prior state.
        sc.close()
        self.assertFalse(self.path.exists())

    def test_n_bytes_written_tracks_appends(self):
        with CarriersSidecar(self.path) as sc:
            self.assertEqual(sc.n_bytes_written, 0)
            sc.write(np.array([[1, 1], [2, 1]], dtype=np.int32))
            self.assertEqual(sc.n_bytes_written, 16)
            sc.write(np.zeros((0, 2), dtype=np.int32))
            self.assertEqual(sc.n_bytes_written, 16)
            sc.write(np.array([[3, 1]], dtype=np.int32))
            self.assertEqual(sc.n_bytes_written, 24)


class CarriersSidecarErrorPathsTest(unittest.TestCase):
    """Pathological inputs surface clear errors rather than corrupt
    silently."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "errors.spill"

    def tearDown(self):
        self.tmp.cleanup()

    def test_wrong_dtype_rejected(self):
        with CarriersSidecar(self.path) as sc:
            bad = np.array([[1, 1]], dtype=np.int64)
            with self.assertRaises(AssertionError) as ctx:
                sc.write(bad)
            self.assertIn("int32", str(ctx.exception))

    def test_wrong_shape_rejected(self):
        # 1-D carriers (the pre-Fix-A list-of-tuples shape) must be
        # caught — they'd silently corrupt the sidecar's byte
        # layout.
        with CarriersSidecar(self.path) as sc:
            bad = np.array([0, 1, 2], dtype=np.int32)
            with self.assertRaises(AssertionError):
                sc.write(bad)

    def test_short_read_raises(self):
        # Reading more bytes than written must raise OSError so
        # programmer bugs in offset/length tracking surface
        # immediately rather than producing garbage carriers.
        with CarriersSidecar(self.path) as sc:
            sc.write(np.array([[1, 1]], dtype=np.int32))
            with self.assertRaises(OSError) as ctx:
                sc.read(0, 999)
            self.assertIn("short read", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
