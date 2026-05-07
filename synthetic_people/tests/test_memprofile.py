"""Tests for syntheticgen/memprofile.py — opt-in RSS sampler.

The profiler is a diagnostic tool, not part of the production data
flow. The contract worth pinning down:

- The TSV header is well-formed (six columns including the
  children/total rss).
- Periodic samples land at roughly the requested cadence.
- ``mark(label)`` from the parent emits a labelled row.
- ``mark(label)`` from a forked child is a no-op (the file handle
  isn't process-safe to share; we don't want children racing on
  writes).
- File data is durable mid-run — every sample is flushed + fsynced
  so an OOM kill doesn't lose the trace.

Tests gate on ``psutil`` being importable. The whole module is
opt-in via the CLI flag, so a host without psutil simply can't use
the profiler — that's the documented behaviour.
"""

from __future__ import annotations

import importlib.util
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

_HAVE_PSUTIL = importlib.util.find_spec("psutil") is not None


def _read_tsv(path: Path) -> list:
    """Parse a memprofile TSV into a list of (header, rows) pairs."""
    text = path.read_text()
    lines = [l for l in text.splitlines() if l]
    header = lines[0].split("\t")
    rows = [l.split("\t") for l in lines[1:]]
    return header, rows


@unittest.skipUnless(_HAVE_PSUTIL, "psutil not installed")
class MemoryProfilerHeaderTest(unittest.TestCase):
    """Header columns and ordering are part of the contract — a
    downstream pandas/awk consumer relies on them."""

    def test_header_has_expected_columns(self):
        from syntheticgen.memprofile import MemoryProfiler
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "mem.tsv"
            with MemoryProfiler(out, sample_interval_s=0.1):
                time.sleep(0.05)
            header, _ = _read_tsv(out)
            self.assertEqual(
                header,
                ["elapsed_s", "rss_mb", "vms_mb",
                 "children_rss_mb", "total_rss_mb", "label"],
            )


@unittest.skipUnless(_HAVE_PSUTIL, "psutil not installed")
class MemoryProfilerSamplingTest(unittest.TestCase):
    """The background thread should fire roughly at the configured
    cadence and produce parseable rows."""

    def test_periodic_samples_within_cadence(self):
        from syntheticgen.memprofile import MemoryProfiler
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "mem.tsv"
            with MemoryProfiler(out, sample_interval_s=0.1):
                # Run for ~0.5 s; expect at least 3 periodic samples
                # plus the start + stop marks.
                time.sleep(0.5)
            _, rows = _read_tsv(out)
            self.assertGreater(len(rows), 3,
                               f"expected >3 rows, got {len(rows)}")
            # Every row parses to (float, float, float, float, float, str).
            for row in rows:
                self.assertEqual(len(row), 6)
                float(row[0])  # elapsed
                float(row[1])  # rss
                float(row[2])  # vms
                float(row[3])  # children rss
                float(row[4])  # total rss
                # row[5] is a string, may be empty


@unittest.skipUnless(_HAVE_PSUTIL, "psutil not installed")
class MemoryProfilerMarkTest(unittest.TestCase):
    """Manual mark() calls land in the TSV with the right label."""

    def test_marks_round_trip_via_module_level_install(self):
        from syntheticgen.memprofile import (
            MemoryProfiler, install, mark,
        )
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "mem.tsv"
            profiler = MemoryProfiler(out, sample_interval_s=10.0)
            profiler.start()
            try:
                install(profiler)
                mark("hello world")
                mark("phase 1 done")
            finally:
                profiler.stop()
                install(None)
            _, rows = _read_tsv(out)
            labels = [row[5] for row in rows]
            self.assertIn("start", labels)  # auto-mark from start()
            self.assertIn("hello world", labels)
            self.assertIn("phase 1 done", labels)
            self.assertIn("stop", labels)

    def test_mark_when_no_profiler_installed_is_noop(self):
        # The module-level mark() should not raise when no profiler
        # is registered — that's the production fast-path when the
        # CLI flag isn't used.
        from syntheticgen.memprofile import install, mark
        install(None)
        mark("nobody home")  # should silently return

    def test_label_with_tabs_and_newlines_sanitised(self):
        # TSV is the wire format; tab/newline in a label would break
        # downstream parsing. memprofile strips them before writing.
        from syntheticgen.memprofile import MemoryProfiler
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "mem.tsv"
            profiler = MemoryProfiler(out, sample_interval_s=10.0)
            profiler.start()
            try:
                profiler.mark("with\ttab\nand\rnewline")
            finally:
                profiler.stop()
            text = out.read_text()
            # No data line should contain a literal tab inside the
            # label column. Easiest check: every data row splits to
            # exactly 6 tab-separated fields.
            for line in text.splitlines()[1:]:
                if not line:
                    continue
                self.assertEqual(len(line.split("\t")), 6,
                                 f"row split count != 6 in {line!r}")


@unittest.skipUnless(_HAVE_PSUTIL, "psutil not installed")
class MemoryProfilerForkSafetyTest(unittest.TestCase):
    """Forked workers must not write to the parent's TSV — the
    file handle isn't process-safe and concurrent writes from
    both sides would interleave randomly.

    The module-level ``mark()`` checks ``os.getpid()`` against the
    profiler's recorded parent pid; mark calls from a fork are
    silently ignored.
    """

    def test_mark_from_fork_is_noop(self):
        import multiprocessing as mp
        from syntheticgen.memprofile import (
            MemoryProfiler, install, mark,
        )
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "mem.tsv"
            profiler = MemoryProfiler(out, sample_interval_s=10.0)
            profiler.start()
            try:
                install(profiler)

                # Mark from parent — should land.
                mark("parent before fork")

                # Fork a worker that calls mark; it should be a no-op.
                ctx = mp.get_context("fork")
                proc = ctx.Process(target=_child_mark)
                proc.start()
                proc.join()

                mark("parent after fork")
            finally:
                profiler.stop()
                install(None)

            _, rows = _read_tsv(out)
            labels = [row[5] for row in rows]
            self.assertIn("parent before fork", labels)
            self.assertIn("parent after fork", labels)
            # No "from child" label — it was suppressed in the fork.
            self.assertNotIn("from child", labels)


def _child_mark() -> None:
    """Top-level so the fork can pickle it (under spawn) or import it
    fresh (under fork)."""
    from syntheticgen.memprofile import mark
    mark("from child")


@unittest.skipUnless(_HAVE_PSUTIL, "psutil not installed")
class MemoryProfilerDurabilityTest(unittest.TestCase):
    """Every write should be flushed + fsynced so an OOM mid-run
    preserves the trace up to the latest sample."""

    def test_data_visible_to_other_readers_after_each_write(self):
        # We can't easily simulate an OOM, but we can check that
        # data is on disk *before* stop() runs — i.e. another reader
        # can read the file mid-run.
        from syntheticgen.memprofile import MemoryProfiler
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "mem.tsv"
            profiler = MemoryProfiler(out, sample_interval_s=0.05)
            profiler.start()
            try:
                # Wait for at least one periodic sample to land.
                time.sleep(0.2)
                # Read the file *now*, while the profiler is still
                # running. Should already see ≥ 2 rows (header +
                # at least one sample plus the start mark).
                with open(out) as fh:
                    contents = fh.read()
                self.assertGreaterEqual(
                    contents.count("\n"), 2,
                    f"expected at least 2 lines visible mid-run, "
                    f"got {contents!r}",
                )
            finally:
                profiler.stop()


if __name__ == "__main__":
    unittest.main()
