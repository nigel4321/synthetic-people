#!/usr/bin/env python3
"""Phase 3 measure-first baseline: peak RSS and per-site gts overhead.

Runs ``simulate_cohort`` for one ``--n`` value, then reports:

* the number of cohort sites msprime emitted;
* the parent process peak RSS (``getrusage(RUSAGE_SELF).ru_maxrss``);
* the bytes occupied by each site's ``gts: list[str]`` field —
  deduplicating string objects by ``id()`` so CPython's small-string
  interning is reflected in the realised cost;
* the gts share of peak RSS as a percentage.

The decision gate in ``PERFORMANCE_PLAN.md`` is: if the gts share is
below 20%, the numpy-uint8 refactor is not worth its tradeoffs and we
land ``sys.intern`` instead. Run for ``n ∈ {200, 500, 1000}`` and
fill in the Baseline table.

Each run is its own process so peak RSS is per-run (process water-mark,
not session-lifetime). Invoke once per ``n``::

    .venv/bin/python synthetic_people/scripts/profile_memory.py --n 200
    .venv/bin/python synthetic_people/scripts/profile_memory.py --n 500
    .venv/bin/python synthetic_people/scripts/profile_memory.py --n 1000

Output is one ``key=value`` line per metric so the wrapper script can
grep / parse it without a full JSON dependency.
"""

from __future__ import annotations

import argparse
import gc
import json
import platform
import random
import resource
import sys
import time
from pathlib import Path

# Allow `from syntheticgen.* import ...` when invoked as a plain script.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from syntheticgen.coalescent import simulate_cohort  # noqa: E402


def gts_overhead_bytes(cohort_sites: list) -> tuple[int, int]:
    """Return ``(total_bytes, unique_strings)`` for the gts field.

    Deduplicates string objects by ``id()`` so an interned ``"0|1"``
    that appears across thousands of sites only costs us once. This is
    the realised RAM cost — the worst-case (no interning) figure can
    be derived as ``sum(sys.getsizeof(s) for site in sites for s in
    site["gts"])`` if needed.
    """
    seen_strings: set = set()
    total = 0
    for site in cohort_sites:
        gts = site.get("gts")
        if gts is None:
            continue
        total += sys.getsizeof(gts)
        for s in gts:
            sid = id(s)
            if sid in seen_strings:
                continue
            seen_strings.add(sid)
            total += sys.getsizeof(s)
    return total, len(seen_strings)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, required=True,
                    help="Cohort size")
    ap.add_argument("--chr-length-mb", type=float, default=5.0,
                    help="Simulated prefix per chromosome (Mb). 0 = full")
    ap.add_argument("--chromosomes", default="22",
                    help="Comma-separated chromosomes (default: 22)")
    ap.add_argument("--demo-model", default="OutOfAfrica_3G09",
                    help="stdpopsim model id; 'none' = constant-Ne msprime")
    ap.add_argument("--population", default="CEU")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    chromosomes = [c.strip() for c in args.chromosomes.split(",")
                   if c.strip()]
    demo_model = None if args.demo_model.lower() == "none" \
        else args.demo_model

    print(f"# host: {platform.platform()}", file=sys.stderr)
    print(f"# python: {sys.version.split()[0]}", file=sys.stderr)
    print(f"# n={args.n} chrs={chromosomes} length={args.chr_length_mb}Mb "
          f"model={demo_model} pop={args.population}",
          file=sys.stderr)

    rng = random.Random(args.seed)

    t0 = time.monotonic()
    cohort_sites = simulate_cohort(
        chromosomes=chromosomes, build="GRCh38",
        n_people=args.n, length_mb=args.chr_length_mb,
        demo_model=demo_model, population=args.population,
        rec_rate=1e-8, mu=1.29e-8, rng=rng,
        verbose=False, workers=1,
    )
    sim_secs = time.monotonic() - t0

    # Force a collection so any temporary tree-sequence garbage is
    # released before we read the high-water mark.
    gc.collect()
    peak_rss_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    n_sites = len(cohort_sites)
    gts_bytes, n_unique_gts = gts_overhead_bytes(cohort_sites)

    peak_rss_bytes = peak_rss_kb * 1024
    gts_share_pct = (gts_bytes / peak_rss_bytes * 100.0
                     if peak_rss_bytes > 0 else 0.0)

    summary = {
        "n": args.n,
        "n_sites": n_sites,
        "sim_secs": round(sim_secs, 2),
        "peak_rss_mb": round(peak_rss_kb / 1024.0, 1),
        "gts_overhead_mb": round(gts_bytes / (1024 * 1024), 2),
        "gts_share_pct": round(gts_share_pct, 1),
        "unique_gt_strings": n_unique_gts,
        "host": platform.node(),
    }
    # Human-readable lines on stderr, machine-readable JSON on stdout.
    for k, v in summary.items():
        print(f"  {k}: {v}", file=sys.stderr)
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
