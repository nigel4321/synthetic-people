#!/usr/bin/env python3
"""
Spike 1 — OS-level mmap+fork smoke test.

Validates that np.memmap + multiprocessing.fork shares physical pages
across workers via the OS page cache. This is the foundational
assumption behind Phase 5d in PERFORMANCE_PLAN.md (streaming-mmap
cohort intermediate, path to n=1M).

PASS:  total system RAM during worker reads grows by ~1 × file_size
       above baseline. Workers all read shared page-cache pages.
FAIL:  total system RAM grows by ~n_workers × file_size. Workers each
       hold private copies; 5d's mmap-share assumption does not hold
       on this host.

See README.md alongside this script for full context.
"""

import argparse
import multiprocessing as mp
import sys
import threading
import time
from pathlib import Path

import numpy as np
import psutil

CHUNK_BYTES = 64 * 1024 * 1024


def make_test_file(path: Path, size_gb: float) -> int:
    """Write `size_gb` GB of random bytes to `path`. Reuse if size matches."""
    size_bytes = int(size_gb * 1024**3)
    if path.exists() and path.stat().st_size == size_bytes:
        print(f"reusing existing {path} ({size_bytes / 1024**3:.2f} GB)")
        return size_bytes
    print(f"writing {size_bytes / 1024**3:.2f} GB of random bytes to {path}...")
    rng = np.random.default_rng(42)
    written = 0
    with path.open("wb") as f:
        while written < size_bytes:
            remaining = size_bytes - written
            n = min(remaining, CHUNK_BYTES)
            chunk = rng.integers(0, 256, size=n, dtype=np.uint8)
            f.write(chunk.tobytes())
            written += n
    print(f"wrote {written / 1024**3:.2f} GB")
    return size_bytes


def worker(path: str, slice_lo: int, slice_hi: int, queue: mp.Queue, idx: int) -> None:
    """Worker: mmap, read its byte slice, return checksum + own RSS."""
    arr = np.memmap(path, dtype=np.uint8, mode="r")
    s = arr[slice_lo:slice_hi]
    checksum = int(s.sum(dtype=np.uint64))
    rss_mb = psutil.Process().memory_info().rss / (1024 * 1024)
    queue.put((idx, checksum, rss_mb))


def sample_system_ram(stop: threading.Event, samples: list, interval_s: float = 0.5) -> None:
    """Background sampler: appends (t, used_gb, available_gb) until stop is set."""
    t0 = time.time()
    while not stop.is_set():
        vm = psutil.virtual_memory()
        samples.append((time.time() - t0, vm.used / (1024**3), vm.available / (1024**3)))
        time.sleep(interval_s)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--size-gb", type=float, default=4.0,
                        help="test file size in GB (default: 4)")
    parser.add_argument("--workers", type=int, default=8,
                        help="number of forked workers (default: 8)")
    parser.add_argument("--path", type=Path, default=Path("/tmp/spike1_mmap_fork.bin"),
                        help="path for the test file (default: /tmp/spike1_mmap_fork.bin)")
    args = parser.parse_args()

    file_size = make_test_file(args.path, args.size_gb)
    file_size_gb = file_size / (1024**3)

    vm0 = psutil.virtual_memory()
    print(f"\nBaseline system RAM: used = {vm0.used / 1024**3:.2f} GB, "
          f"available = {vm0.available / 1024**3:.2f} GB")

    arr = np.memmap(args.path, dtype=np.uint8, mode="r")
    parent_rss_gb = psutil.Process().memory_info().rss / (1024**3)
    vm_after = psutil.virtual_memory()
    print(f"After parent np.memmap: parent RSS = {parent_rss_gb:.2f} GB, "
          f"system used = {vm_after.used / 1024**3:.2f} GB")
    _ = arr  # keep mmap alive across fork

    mp.set_start_method("fork", force=True)
    queue: mp.Queue = mp.Queue()
    slice_size = file_size // args.workers
    procs = []
    for i in range(args.workers):
        lo = i * slice_size
        hi = (i + 1) * slice_size if i < args.workers - 1 else file_size
        proc = mp.Process(target=worker, args=(str(args.path), lo, hi, queue, i))
        procs.append(proc)

    samples: list = []
    stop = threading.Event()
    sampler = threading.Thread(target=sample_system_ram, args=(stop, samples), daemon=True)
    sampler.start()

    print(f"\nForking {args.workers} workers...")
    t0 = time.time()
    for proc in procs:
        proc.start()

    results = []
    for _ in range(args.workers):
        results.append(queue.get())

    for proc in procs:
        proc.join()
    elapsed = time.time() - t0

    stop.set()
    sampler.join(timeout=2)

    print(f"All workers done in {elapsed:.1f} s\n")

    print("Per-worker reports:")
    apparent_total_gb = 0.0
    for idx, checksum, rss_mb in sorted(results):
        rss_gb = rss_mb / 1024
        apparent_total_gb += rss_gb
        print(f"  worker {idx}: checksum={checksum}, reported RSS = {rss_gb:.2f} GB")
    print(f"\nSum of per-worker reported RSS (apparent): {apparent_total_gb:.2f} GB")

    if not samples:
        print("ERROR: no system RAM samples captured; workers ran too fast?")
        return 3

    peak_used_gb = max(s[1] for s in samples)
    min_available_gb = min(s[2] for s in samples)
    baseline_gb = vm0.used / 1024**3
    delta_gb = peak_used_gb - baseline_gb

    print("\nSystem RAM during worker phase:")
    print(f"  baseline used:           {baseline_gb:.2f} GB")
    print(f"  peak used:               {peak_used_gb:.2f} GB")
    print(f"  min available:           {min_available_gb:.2f} GB")
    print(f"  delta vs baseline:       +{delta_gb:.2f} GB")
    print(f"  file size:               {file_size_gb:.2f} GB")
    print(f"  n_workers × file size:   {args.workers * file_size_gb:.2f} GB "
          f"(what we'd see if no sharing)")
    print(f"  apparent (sum of RSS):   {apparent_total_gb:.2f} GB")

    print("\n=== VERDICT ===")
    pass_threshold = file_size_gb * 1.5
    fail_threshold = args.workers * file_size_gb * 0.5

    if delta_gb < pass_threshold:
        print(f"PASS: system RAM grew by {delta_gb:.2f} GB during worker phase,")
        print(f"      below the {pass_threshold:.2f} GB threshold "
              f"(1.5 × file_size).")
        print(f"      Apparent total ({apparent_total_gb:.2f} GB) >> physical "
              f"delta ({delta_gb:.2f} GB)")
        print(f"      confirms fork-COW sharing via the OS page cache.")
        print(f"      Phase 5d's OS-level assumption holds on this host.")
        return 0
    elif delta_gb < fail_threshold:
        print(f"AMBIGUOUS: system RAM grew by {delta_gb:.2f} GB.")
        print(f"      Above the pass threshold ({pass_threshold:.2f} GB) but below")
        print(f"      the clean-fail threshold ({fail_threshold:.2f} GB).")
        print(f"      Some sharing, some divergence. Investigate kernel settings:")
        print(f"      swappiness, transparent hugepages, KSM. See README.md.")
        return 1
    else:
        print(f"FAIL: system RAM grew by {delta_gb:.2f} GB, close to")
        print(f"      n_workers × file_size = {args.workers * file_size_gb:.2f} GB.")
        print(f"      Workers appear to hold private copies. Phase 5d's")
        print(f"      mmap-share assumption does NOT hold on this host.")
        print(f"      Halt 5d planning until the cause is understood.")
        return 2


if __name__ == "__main__":
    sys.exit(main())
