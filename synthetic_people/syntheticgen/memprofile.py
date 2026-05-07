"""Lightweight memory-usage profiler for diagnosing OOM in
``generate_people``.

Opt-in via ``--profile-memory PATH``. When active, a background
thread samples ``psutil.Process().memory_info()`` once per second
and appends one TSV row per sample to ``PATH``. The rest of the
code can call :meth:`MemoryProfiler.mark` from key transition
points (after a costly load, before/after each chromosome's
simulation, around per-chunk msprime calls, …) to drop labelled
checkpoint rows into the same TSV — those rows have a non-empty
``label`` column and capture the RSS at exactly that point.

The TSV is flushed + fsynced after **every** write so an OOM kill
preserves all data up to the point the kernel reaped the process.
That's the usual reason people turn this on, and missing the last
few seconds of data is exactly when you most want them.

Output schema (tab-separated, one header row):

::

    elapsed_s   rss_mb   vms_mb   children_rss_mb   total_rss_mb   label

``rss_mb`` / ``vms_mb`` are the parent process alone;
``children_rss_mb`` is the sum across all live descendants (fork
workers, bcftools subprocesses, etc.); ``total_rss_mb`` is parent +
children — usually what you want for OOM diagnosis since the kernel
budgets against the whole process tree. ``label`` is empty for the
periodic samples and carries the caller-supplied string for
explicit marks. Plot with whatever's handy:

.. code-block:: python

    import pandas as pd, matplotlib.pyplot as plt
    df = pd.read_csv("memprofile.tsv", sep="\\t")
    plt.plot(df.elapsed_s, df.rss_mb)
    for _, row in df[df.label != ""].iterrows():
        plt.axvline(row.elapsed_s, color="red", alpha=0.3)
        plt.text(row.elapsed_s, row.rss_mb, row.label, rotation=90)
    plt.xlabel("seconds"); plt.ylabel("RSS (MB)")
    plt.show()

The profiler is a no-op when ``psutil`` isn't importable — the
flag's help text points the user at ``pip install psutil`` if they
hit it.
"""

from __future__ import annotations

import os
import threading
import time
from pathlib import Path


class MemoryProfiler:
    """Background-thread RSS sampler with mark API.

    Use as a context manager or pair start/stop manually. The mark
    method is thread-safe via an in-memory queue + lock — callers
    don't need to coordinate with the sampling thread directly.
    """

    def __init__(self, out_path: Path, sample_interval_s: float = 1.0):
        self.out_path = Path(out_path)
        self.sample_interval_s = sample_interval_s
        self._fh = None
        self._proc = None
        self._t0 = 0.0
        self._stop_event = threading.Event()
        self._thread: threading.Thread | None = None
        self._write_lock = threading.Lock()
        # Lazy-import psutil so the rest of the module imports
        # without it. The CLI checks at flag-parse time.
        import psutil
        self._psutil = psutil

    def __enter__(self) -> "MemoryProfiler":
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.stop()

    def start(self) -> None:
        """Open the TSV, start the sampler thread."""
        self.out_path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(self.out_path, "w")
        # Header row. Tab-separated for grep / awk / pandas friendliness.
        self._fh.write(
            "elapsed_s\trss_mb\tvms_mb\t"
            "children_rss_mb\ttotal_rss_mb\tlabel\n")
        self._fh.flush()
        os.fsync(self._fh.fileno())
        self._proc = self._psutil.Process(os.getpid())
        self._t0 = time.monotonic()
        self.mark("start")
        self._thread = threading.Thread(
            target=self._sample_loop, name="memprofile", daemon=True,
        )
        self._thread.start()

    def stop(self) -> None:
        """Signal the sampler thread, wait for it, close the file.

        Idempotent — safe to call from a finally block whether or
        not start() succeeded.
        """
        if self._thread is None:
            return
        self.mark("stop")
        self._stop_event.set()
        self._thread.join(timeout=self.sample_interval_s * 2)
        self._thread = None
        if self._fh is not None:
            self._fh.close()
            self._fh = None

    def mark(self, label: str) -> None:
        """Drop a labelled checkpoint row at the current RSS.

        Safe to call from any thread. The label column lets a
        downstream plotter draw a vertical guide at exactly this
        point.
        """
        self._sample_and_write(label)

    def _sample_loop(self) -> None:
        """Periodic RSS sample, runs in a daemon thread."""
        while not self._stop_event.is_set():
            self._sample_and_write("")
            # Wait with stop-aware timeout so a fast stop() doesn't
            # block on a long sleep.
            self._stop_event.wait(self.sample_interval_s)

    def _sample_and_write(self, label: str) -> None:
        """Read parent + children RSS and emit one TSV row.

        Children-RSS is summed across all live descendants — fork
        workers carrying their own tree sequences plus the
        ``bcftools view -O b`` subprocesses we pipe BCFs into. For
        OOM diagnosis the total (parent + children) is what
        matters, since the kernel budgets against the whole tree.
        """
        if self._fh is None or self._proc is None:
            return
        try:
            mem = self._proc.memory_info()
        except Exception:
            # Self-process going away (post-stop). Bail.
            return
        elapsed = time.monotonic() - self._t0
        rss_mb = mem.rss / (1024 * 1024)
        vms_mb = mem.vms / (1024 * 1024)
        children_rss_mb = self._children_rss_mb()
        total_rss_mb = rss_mb + children_rss_mb
        clean = (label.replace("\t", " ").replace("\n", " ")
                      .replace("\r", " "))
        with self._write_lock:
            if self._fh is None:
                return
            self._fh.write(
                f"{elapsed:.3f}\t{rss_mb:.2f}\t{vms_mb:.2f}\t"
                f"{children_rss_mb:.2f}\t{total_rss_mb:.2f}\t"
                f"{clean}\n"
            )
            try:
                self._fh.flush()
                os.fsync(self._fh.fileno())
            except (OSError, ValueError):
                # File may have been closed by a concurrent stop().
                pass

    def _children_rss_mb(self) -> float:
        """Sum RSS of all live descendants in MB.

        Walks the full tree (``recursive=True``) so we catch
        msprime worker processes, bcftools sub-subprocesses,
        anything else spawned. A child going away mid-walk is
        ignored — its RSS at that moment is gone anyway.
        """
        try:
            children = self._proc.children(recursive=True)
        except Exception:
            return 0.0
        total = 0
        for child in children:
            try:
                total += child.memory_info().rss
            except Exception:
                continue
        return total / (1024 * 1024)


# Module-level singleton populated by main() when --profile-memory is
# set. Worker functions in fork-based ProcessPoolExecutor inherit it
# via copy-on-write but should NOT call mark from inside workers —
# the file handle isn't safe to share across processes. Mark calls
# from worker code are silently ignored.
_PROFILER: MemoryProfiler | None = None


def install(profiler: MemoryProfiler | None) -> None:
    """Register a profiler at module level so unrelated callers can
    reach it via :func:`mark` without threading the object through
    every signature."""
    global _PROFILER
    _PROFILER = profiler


def mark(label: str) -> None:
    """Drop a checkpoint via the module-level profiler if one is
    active. No-op otherwise; safe to call from anywhere in the
    parent process."""
    if _PROFILER is None or _PROFILER._proc is None:
        return
    # Defensive: ignore mark calls from forked workers. The file
    # handle isn't process-safe to share, and worker-side marks
    # would mostly capture the worker's COW-shared baseline anyway.
    if os.getpid() != _PROFILER._proc.pid:
        return
    _PROFILER.mark(label)
