#!/usr/bin/env python3
"""
Spike 2b — Phase 5d sub-spike.

Characterise parent peak RSS during Arrow IPC streaming write as a
function of three knobs:

  * batch_size in {256, 1024, 4096}
  * ARROW_DEFAULT_MEMORY_POOL in {jemalloc, system}
  * release_unused-after-each-write in {False, True}

Goal: identify a known-good (batch_size, pool, release_unused) default
for Phase 5d.1 build, before that default is wired into production
code.

Each matrix cell runs in a fresh subprocess so the pyarrow memory pool
starts clean — the pool is a process-global singleton and would
otherwise carry high-water mark across cells.

Default scale: 50k samples x 200k sites = ~10 GB Arrow file. Big
enough to amortise allocator effects above noise; small enough that
the full 12-cell matrix finishes in roughly 15 minutes on a
workstation. Pass --samples / --sites to scale up.

Usage:
    python spike2b_pool_batch_matrix.py
    python spike2b_pool_batch_matrix.py --samples 100000 --sites 500000
    python spike2b_pool_batch_matrix.py --path-dir /mnt/nvme

Exits 3 if pyarrow is missing.
"""

import argparse
import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np

DEFAULT_SAMPLES = 50_000
DEFAULT_SITES = 200_000
DEFAULT_BATCH_SIZES = [256, 1024, 4096]
DEFAULT_POOLS = ["jemalloc", "system"]
DEFAULT_PATH_DIR = "/tmp"
CELL_TIMEOUT_S = 1200


def cell_main(args):
    try:
        import pyarrow as pa
        import pyarrow.ipc as ipc
    except ImportError:
        sys.stderr.write("pyarrow not installed\n")
        sys.exit(3)

    import psutil

    samples = args.samples
    sites = args.sites
    batch_size = args.batch_size
    release_unused = args.release_unused
    path = Path(args.path)

    proc = psutil.Process()
    baseline_mb = proc.memory_info().rss / 1e6

    peak_holder = [baseline_mb]
    stop = threading.Event()

    def sampler():
        while not stop.is_set():
            r = proc.memory_info().rss / 1e6
            if r > peak_holder[0]:
                peak_holder[0] = r
            time.sleep(0.2)

    sampler_thread = threading.Thread(target=sampler, daemon=True)
    sampler_thread.start()

    schema = pa.schema(
        [
            ("pos", pa.int64()),
            ("genotypes", pa.list_(pa.int8(), samples)),
        ]
    )

    pool_backend = pa.default_memory_pool().backend_name

    rng = np.random.default_rng(seed=42)
    n_batches = (sites + batch_size - 1) // batch_size

    t0 = time.time()
    with ipc.new_file(str(path), schema) as writer:
        for b in range(n_batches):
            this_batch = min(batch_size, sites - b * batch_size)
            pos_np = np.arange(
                b * batch_size, b * batch_size + this_batch, dtype=np.int64
            )
            gt_2d = (rng.random((this_batch, samples)) < 0.05).astype(np.int8)
            gt_flat = gt_2d.reshape(-1)

            pos_arr = pa.array(pos_np)
            gt_values = pa.array(gt_flat, type=pa.int8())
            gt_arr = pa.FixedSizeListArray.from_arrays(gt_values, samples)

            batch = pa.RecordBatch.from_arrays(
                [pos_arr, gt_arr], schema=schema
            )
            writer.write_batch(batch)

            del batch, gt_2d, gt_flat, pos_np, pos_arr, gt_values, gt_arr
            if release_unused:
                pa.default_memory_pool().release_unused()

    elapsed = time.time() - t0
    stop.set()
    sampler_thread.join()

    file_size_gb = path.stat().st_size / 1e9
    expected_bytes = samples * sites
    write_thr = expected_bytes / elapsed / 1e6 if elapsed > 0 else 0.0

    print(
        json.dumps(
            {
                "samples": samples,
                "sites": sites,
                "batch_size": batch_size,
                "pool_requested": os.environ.get(
                    "ARROW_DEFAULT_MEMORY_POOL", "default"
                ),
                "pool_backend": pool_backend,
                "release_unused": release_unused,
                "n_batches": n_batches,
                "baseline_mb": round(baseline_mb, 1),
                "peak_mb": round(peak_holder[0], 1),
                "growth_mb": round(peak_holder[0] - baseline_mb, 1),
                "elapsed_s": round(elapsed, 2),
                "file_size_gb": round(file_size_gb, 2),
                "write_throughput_mb_s": round(write_thr, 1),
            }
        )
    )


def run_matrix(args):
    matrix = []
    for pool in args.pools:
        for batch_size in args.batch_sizes:
            release_axis = (False, True) if args.with_release_unused else (False,)
            for release_unused in release_axis:
                matrix.append((pool, batch_size, release_unused))

    raw_gb = args.samples * args.sites / 1e9
    print(f"Spike 2b — pool x batch_size matrix")
    print(f"  cells: {len(matrix)}")
    print(f"  scale: {args.samples} samples x {args.sites} sites = ~{raw_gb:.1f} GB raw GT bytes")
    print(f"  pools: {args.pools}")
    print(f"  batch sizes: {args.batch_sizes}")
    print(f"  release_unused axis: {'on' if args.with_release_unused else 'off'}")
    print(f"  per-cell timeout: {CELL_TIMEOUT_S} s")
    print(f"  cell file: {args.path_dir}/spike2b_cell.arrow (deleted between cells)")
    print()

    path = Path(args.path_dir) / "spike2b_cell.arrow"
    results = []

    for i, (pool, batch_size, release_unused) in enumerate(matrix):
        print(
            f"--- cell {i+1}/{len(matrix)}: pool={pool} batch_size={batch_size} "
            f"release_unused={release_unused} ---"
        )

        env = os.environ.copy()
        env["ARROW_DEFAULT_MEMORY_POOL"] = pool

        cmd = [
            sys.executable,
            __file__,
            "--cell",
            "--samples", str(args.samples),
            "--sites", str(args.sites),
            "--batch-size", str(batch_size),
            "--path", str(path),
        ]
        if release_unused:
            cmd.append("--release-unused")

        cell_t0 = time.time()
        try:
            r = subprocess.run(
                cmd,
                env=env,
                capture_output=True,
                text=True,
                timeout=CELL_TIMEOUT_S,
            )
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT after {CELL_TIMEOUT_S} s")
            results.append(
                {
                    "pool_requested": pool,
                    "batch_size": batch_size,
                    "release_unused": release_unused,
                    "error": "timeout",
                }
            )
            if path.exists():
                path.unlink()
            print()
            continue
        cell_elapsed = time.time() - cell_t0

        if r.returncode != 0:
            print(f"  FAILED rc={r.returncode}")
            if r.stderr:
                print(f"  stderr: {r.stderr.strip()[:500]}")
            if r.returncode == 3:
                print("Aborting matrix — pyarrow missing in this env.")
                sys.exit(3)
            results.append(
                {
                    "pool_requested": pool,
                    "batch_size": batch_size,
                    "release_unused": release_unused,
                    "error": f"rc={r.returncode}",
                }
            )
            if path.exists():
                path.unlink()
            print()
            continue

        try:
            cell_result = json.loads(r.stdout.strip().splitlines()[-1])
        except (json.JSONDecodeError, IndexError) as e:
            print(f"  PARSE ERROR: {e}")
            print(f"  stdout (first 500 chars): {r.stdout[:500]}")
            results.append(
                {
                    "pool_requested": pool,
                    "batch_size": batch_size,
                    "release_unused": release_unused,
                    "error": "parse",
                }
            )
            if path.exists():
                path.unlink()
            print()
            continue

        results.append(cell_result)
        print(
            f"  pool_backend={cell_result['pool_backend']}  "
            f"peak={cell_result['peak_mb']} MB  "
            f"growth=+{cell_result['growth_mb']} MB  "
            f"thr={cell_result['write_throughput_mb_s']} MB/s  "
            f"elapsed={cell_result['elapsed_s']} s  "
            f"(wallclock {cell_elapsed:.1f} s)"
        )

        if path.exists():
            path.unlink()
        print()

    print_summary_table(results, args)


def print_summary_table(results, args):
    print("=" * 78)
    print("=== Summary ===")
    raw_gb = args.samples * args.sites / 1e9
    print(f"  scale: {args.samples} samples x {args.sites} sites = ~{raw_gb:.1f} GB raw")
    print()

    header = (
        f"  {'pool':<10} {'backend':<10} {'batch':>6} {'rel_unused':>10} "
        f"{'peak_MB':>10} {'growth_MB':>10} {'thr_MB/s':>10} {'elapsed_s':>10}"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))
    for r in results:
        if "error" in r:
            print(
                f"  {r['pool_requested']:<10} {'-':<10} {r['batch_size']:>6} "
                f"{str(r['release_unused']):>10}  ERROR ({r['error']})"
            )
        else:
            print(
                f"  {r['pool_requested']:<10} {r['pool_backend']:<10} "
                f"{r['batch_size']:>6} {str(r['release_unused']):>10} "
                f"{r['peak_mb']:>10} {r['growth_mb']:>10} "
                f"{r['write_throughput_mb_s']:>10} {r['elapsed_s']:>10}"
            )

    valid = [r for r in results if "error" not in r]
    if valid:
        lowest_growth = min(valid, key=lambda r: r["growth_mb"])
        fastest = max(valid, key=lambda r: r["write_throughput_mb_s"])
        print()
        print(
            f"  Lowest peak growth: pool={lowest_growth['pool_requested']} "
            f"batch_size={lowest_growth['batch_size']} "
            f"release_unused={lowest_growth['release_unused']} "
            f"-> +{lowest_growth['growth_mb']} MB at "
            f"{lowest_growth['write_throughput_mb_s']} MB/s"
        )
        print(
            f"  Fastest write:      pool={fastest['pool_requested']} "
            f"batch_size={fastest['batch_size']} "
            f"release_unused={fastest['release_unused']} "
            f"-> {fastest['write_throughput_mb_s']} MB/s at "
            f"+{fastest['growth_mb']} MB peak"
        )
    print()
    print("Note: cells where pool_backend != pool_requested mean pyarrow")
    print("fell back to a different pool because the requested one is not")
    print("compiled in. release_unused=True on the system pool is a no-op.")


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--cell", action="store_true", help=argparse.SUPPRESS
    )
    parser.add_argument("--samples", type=int, default=DEFAULT_SAMPLES)
    parser.add_argument("--sites", type=int, default=DEFAULT_SITES)
    parser.add_argument(
        "--batch-size", type=int, default=1024, help="cell mode only"
    )
    parser.add_argument(
        "--release-unused", action="store_true", help="cell mode only"
    )
    parser.add_argument("--path", type=str, help="cell mode only")
    parser.add_argument(
        "--path-dir",
        type=str,
        default=DEFAULT_PATH_DIR,
        help="driver: directory for the per-cell Arrow file (deleted between cells)",
    )
    parser.add_argument(
        "--batch-sizes",
        type=int,
        nargs="+",
        default=DEFAULT_BATCH_SIZES,
        help="driver: batch sizes to sweep",
    )
    parser.add_argument(
        "--pools",
        type=str,
        nargs="+",
        default=DEFAULT_POOLS,
        help="driver: memory pools to sweep (jemalloc, system, mimalloc)",
    )
    parser.add_argument(
        "--with-release-unused",
        dest="with_release_unused",
        action="store_true",
        default=True,
        help="driver: include release_unused True/False axis (default on; 12 cells)",
    )
    parser.add_argument(
        "--no-release-unused-axis",
        dest="with_release_unused",
        action="store_false",
        help="driver: omit release_unused axis (6 cells, all release_unused=False)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.cell:
        if not args.path:
            sys.stderr.write("--path required in cell mode\n")
            sys.exit(2)
        cell_main(args)
    else:
        run_matrix(args)


if __name__ == "__main__":
    main()
