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

    def test_init_failure_cleans_up_partial_state(self):
        # PR #77 review: ``__init__`` opens two fds sequentially —
        # if the second open fails after the first succeeded, the
        # first must be released and the file unlinked. Patch
        # ``os.open`` to raise so the second fd acquisition fails;
        # ``open()`` (the first call, for the write fd) is left
        # alone. After the exception, no file should remain.
        from unittest.mock import patch
        with patch(
            "syntheticgen.carriers_sidecar.os.open",
            side_effect=OSError("simulated read-fd open failure"),
        ):
            with self.assertRaises(OSError):
                CarriersSidecar(self.path)
        # File created by the write-fd open() must have been
        # unlinked when the read-fd open() raised. Otherwise the
        # cli's per-chrom retry path would see a phantom sidecar
        # from a prior crash.
        self.assertFalse(self.path.exists())


class CarriersSidecarShortWriteTest(unittest.TestCase):
    """PR #77 review #3: the unbuffered write fd's ``write()`` can
    return fewer bytes than requested under EINTR / disk-full /
    signal interruption. The implementation loops until all bytes
    are written, but if no progress is made (true 0-byte write)
    it must raise rather than spin forever."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "short_write.spill"

    def tearDown(self):
        self.tmp.cleanup()

    def test_partial_writes_are_retried(self):
        # Simulate a series of partial writes (e.g. one byte at a
        # time) — ``write()`` must loop until all bytes land, and
        # the read-back must produce the full payload.
        from unittest.mock import patch
        with CarriersSidecar(self.path) as sc:
            real_write = sc._write_fh.write
            call_count = [0]
            def short_write(data: memoryview) -> int:
                # First three calls advance one byte each; later
                # calls write everything remaining.
                call_count[0] += 1
                if call_count[0] <= 3:
                    return real_write(bytes(data[:1]))
                return real_write(bytes(data))
            with patch.object(sc._write_fh, "write",
                              side_effect=short_write):
                c = np.array([[1, 1], [2, 1], [3, 1]], dtype=np.int32)
                # 24 bytes total; first 3 calls write 1 byte each
                # (3 bytes), call 4 finishes the remaining 21.
                ref = sc.write(c)
                self.assertEqual(ref, (0, 24))
                self.assertGreaterEqual(call_count[0], 4)
            back = sc.read(*ref)
            self.assertTrue(np.array_equal(back, c))

    def test_zero_byte_write_raises(self):
        # A write that makes literally no progress would loop
        # forever without the guard. Stub returns 0 to simulate
        # this pathological case.
        from unittest.mock import patch
        with CarriersSidecar(self.path) as sc:
            with patch.object(sc._write_fh, "write", return_value=0):
                with self.assertRaises(OSError) as ctx:
                    sc.write(np.array([[1, 1]], dtype=np.int32))
                self.assertIn("no progress", str(ctx.exception))


class CarriersSidecarErrorPathsTest(unittest.TestCase):
    """Pathological inputs surface clear errors rather than corrupt
    silently."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "errors.spill"

    def tearDown(self):
        self.tmp.cleanup()

    def test_wrong_dtype_rejected(self):
        # PR #77 review: dtype/shape validation must raise real
        # exceptions (not ``assert``) so they fire under
        # ``python -O`` too. Silent corruption-mode is not OK
        # because the sidecar's byte layout depends on the input
        # shape.
        with CarriersSidecar(self.path) as sc:
            bad = np.array([[1, 1]], dtype=np.int64)
            with self.assertRaises(TypeError) as ctx:
                sc.write(bad)
            self.assertIn("int32", str(ctx.exception))

    def test_wrong_shape_rejected(self):
        # 1-D carriers (the pre-Fix-A list-of-tuples shape) must be
        # caught — they'd silently corrupt the sidecar's byte
        # layout. ``ValueError`` not ``AssertionError`` so the
        # check stands under ``python -O``.
        with CarriersSidecar(self.path) as sc:
            bad = np.array([0, 1, 2], dtype=np.int32)
            with self.assertRaises(ValueError):
                sc.write(bad)

    def test_short_read_raises(self):
        # Reading more bytes than written must raise OSError so
        # programmer bugs in offset/length tracking surface
        # immediately rather than producing garbage carriers.
        with CarriersSidecar(self.path) as sc:
            sc.write(np.array([[1, 1]], dtype=np.int32))
            with self.assertRaises(OSError) as ctx:
                sc.read(0, 24)  # 3 rows requested, only 1 written
            self.assertIn("short read", str(ctx.exception))

    def test_non_multiple_read_length_raises(self):
        # ``read(offset, n_bytes)`` requires ``n_bytes`` be a
        # multiple of 8 (the per-row size: 2 × int32). A non-
        # multiple would previously have hit numpy's reshape with
        # a confusing "cannot reshape array of size N into shape
        # (M, 2)" error; now it raises ValueError with a clear
        # message.
        with CarriersSidecar(self.path) as sc:
            sc.write(np.array([[1, 1], [2, 1]], dtype=np.int32))
            with self.assertRaises(ValueError) as ctx:
                sc.read(0, 12)  # not a multiple of 8
            self.assertIn("multiple of", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
