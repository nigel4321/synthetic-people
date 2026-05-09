# Phase 5d viability spikes

Small standalone scripts that validate the load-bearing assumptions
behind **Phase 5d** (`PERFORMANCE_PLAN.md` — _streaming-mmap cohort
intermediate, path to n=1M_) before we commit to the ~600–800 LOC
build.

The 5d design rests on two assumptions that are easy to *believe*
and hard to *prove* without running real code:

1. **`mmap` + `fork` shares physical pages across workers** via the
   OS page cache. If true, parent + W workers can all read the same
   large file backed by ~1 × file-size of physical RAM, not
   W × file-size.
2. **`pyarrow` IPC reader actually does zero-copy mmap** for our
   schema and access pattern. If a worker's `read_all()` quietly
   loads buffers into anonymous memory, we're back to the
   refcount-COW problem that killed Phase 5e Phase A at n=3000.

If either assumption fails, 5d's architecture needs rethinking
before we write 600+ LOC of changes to the synthetic_people
pipeline. These spikes test each assumption in isolation so we
learn early.

## Spike layout

| Spike | Status | What it tests | Cost |
|---|---|---|---|
| 1 | ✅ shipped, **PASS** (2026-05-09) | OS-level: `np.memmap` + `fork` + 8 workers shares physical RAM | ~120 LOC |
| 2 | ✅ shipped | `pyarrow` IPC streaming write + zero-copy mmap read across workers | ~280 LOC |

Run Spike 1 first. If it fails (workers don't share physical RAM),
Spike 2 is moot — 5d's foundation is wrong and we rethink. If
Spike 1 passes, run Spike 2 to validate the Arrow-specific layer.

Spike 1 result on the user's host (32 GB workstation) recorded in
`spike1_results_2026-05-09.txt`: clean PASS — total system RAM
delta during the 8-worker read phase was 0 GB, apparent per-worker
RSS sum (4.21 GB) tracked file size rather than `n_workers ×
file_size`. Green-lit Spike 2.

## Spike 1 — OS-level mmap+fork smoke test

**File:** [`spike1_mmap_fork_smoke.py`](spike1_mmap_fork_smoke.py)

**Result:** [`spike1_results_2026-05-09.txt`](spike1_results_2026-05-09.txt)
— **PASS** on the user's 32 GB workstation (2026-05-09).

**What it does:**

1. Writes a 4 GB file of random bytes to `/tmp` (or the path you
   pass via `--path`).
2. Parent process opens it via `np.memmap`.
3. Forks N workers (default 8) using
   `multiprocessing.set_start_method("fork")`.
4. Each worker reads its own non-overlapping byte slice and
   computes a checksum (forces page-in of the slice's pages).
5. While workers run, samples `psutil.virtual_memory().used`
   every 500 ms in a background thread.
6. Reports: per-worker reported RSS, peak total system RAM during
   read phase, baseline-relative delta, and a pass / ambiguous /
   fail verdict.

**Pass criteria:**

- Peak total system RAM during the worker phase grows by less than
  ~1.5 × file_size above baseline. That means the kernel is
  sharing the mmap'd pages across workers via the page cache —
  exactly what 5d needs.

**Fail criteria:**

- Peak total system RAM grows by close to `n_workers × file_size`
  above baseline. Workers each materialised a private copy of the
  file. 5d's `mmap` story is broken on this host and we do not
  proceed with 5d.1.

**Ambiguous:**

- Total RAM growth between 1.5× and `0.5 × n_workers ×` file_size.
  Some sharing, some divergence. Investigate kernel version,
  swappiness, transparent hugepage settings, KSM, before deciding.

### How to run

```bash
# from repo root, with synthetic_people venv active (numpy + psutil)
python synthetic_people/scripts/spikes/spike1_mmap_fork_smoke.py

# larger / smaller test
python synthetic_people/scripts/spikes/spike1_mmap_fork_smoke.py \
    --size-gb 8 --workers 16

# custom path (e.g. on a fast NVMe vs default /tmp)
python synthetic_people/scripts/spikes/spike1_mmap_fork_smoke.py \
    --path /mnt/nvme/spike1.bin --size-gb 16
```

The script reuses the test file across runs if its size matches,
so re-running is cheap.

### Interpreting noisy results

Linux page-cache behaviour can be affected by:

- **Swappiness:** if `cat /proc/sys/vm/swappiness` > 60 and the
  host is RAM-constrained, the kernel may evict mmap'd pages
  during the test, causing repeated re-page-in that inflates
  apparent RAM. Run on an idle host or temporarily set
  `swappiness=10`.
- **Transparent hugepages (THP):** can cause coarser-grained
  page accounting. Check `cat /sys/kernel/mm/transparent_hugepage/
  enabled`. `madvise` or `never` is the cleanest test bed; `always`
  may inflate the per-worker reported RSS even though physical use
  is small.
- **KSM (Kernel Same-page Merging):** opposite direction —
  deduplicates identical pages across processes, which would
  *help* sharing but isn't enabled by default on most distros.
- **Other workloads on the host:** any process churning memory
  during the test will pollute the `psutil.virtual_memory()`
  delta. Run on a quiet host.

### What "pass" tells us about 5d

A clean pass means the kernel-level mechanism behind 5d works on
your host. It does *not* yet validate:

- That `pyarrow.ipc` actually exposes mmap'd buffers without
  copying (Spike 2's job).
- That worker-side iteration over numpy mmap views doesn't
  trigger COW divergence through some unexpected write (e.g.,
  a numpy operation that writes to a temporary buffer in shared
  pages). Spike 2 covers this with realistic per-site iteration.

So Spike 1 is necessary but not sufficient. Pass → green-light
Spike 2; fail → halt 5d planning until we understand why.

## Spike 2 — Arrow IPC streaming write + zero-copy mmap

**File:** [`spike2_arrow_streaming.py`](spike2_arrow_streaming.py)

**Result:** _not yet run_ — paste stdout into
`spike2_results_<date>.txt` alongside this README and link it here
once available.

**What it does:**

1. Synthesises data of realistic shape on the fly:
   - Default: 10,000 samples × 100,000 sites of int8 genotypes
   - Schema: `pos int64 | genotypes FixedSizeList(int8, n_samples)`
     — matches the layout tentatively locked in by the 5d plan
     (per-row fixed-size list = effectively a row-major int8
     matrix, with workers slicing the columns dimension via
     numpy's `.reshape(-1, n_samples)[:, lo:hi]` view)
   - ~1 GB raw GT bytes; Arrow IPC file slightly larger due to
     metadata
   - Allele frequency: ~5 % non-ref (sparse-ish, realistic for
     common variants)
2. **Streaming write phase:** parent generates one record batch
   at a time (`--batch-size`, default 1024 sites) and writes
   each via `pyarrow.ipc.new_file(path, schema).write_batch(...)`.
   A background thread samples parent RSS at 200 ms intervals
   throughout the write to catch any growth.
3. **Worker mmap-read phase:** fork N workers (default 8); each:
   - Opens the file via `pa.memory_map(path, "r")` →
     `pa.ipc.open_file(mm)` (explicit mmap, not silently buffered)
   - Iterates record batches; for each batch, takes the
     `genotypes` FixedSizeListArray, calls
     `.values.to_numpy(zero_copy_only=True)` (raises if a copy
     would be needed — the explicit guarantee), reshapes to
     `(batch_n_sites, n_samples)`, and slices its sample range
     (zero-copy view).
   - Iterates per-site (Python `for` loop over the matrix rows)
     and computes per-site alt counts. **The per-site loop is
     deliberate** — it mirrors the realistic 5d.1 worker pattern
     (per-site BCF write) and validates that this access pattern
     doesn't trigger PyObject allocation per element or COW
     divergence.
4. **Measures:**
   - Parent peak RSS during streaming write (target: < 500 MB)
   - Total system RSS during worker phase (target: ≈ file size,
     not W × file size)
   - Per-worker reported RSS
   - Write throughput (MB/s logical)
   - Aggregate read throughput across workers (MB/s)

**Pass criteria:**

- Parent peak RSS during streaming write < 500 MB regardless of
  how many sites are streamed.
- Total system RSS during worker reads stays close to the Arrow
  file size; does NOT scale with worker count.
- Per-worker iteration completes without OOM.
- Write throughput > 200 MB/s on NVMe (so an 80 GB write at n=1M
  finishes in < 10 minutes per chromosome).
- Aggregate read throughput > 500 MB/s (so workers complete I/O
  in reasonable time).

**Fail modes that would change the 5d design:**

1. **Parent RAM grows with sites written** during streaming write
   → `pyarrow.ipc.new_file` doesn't actually stream. Investigate
   alternative writers (chunked Arrow files, raw struct format,
   bypass Arrow entirely).
2. **Total system RAM grows per worker** during read → "zero-copy
   mmap" isn't, for our schema. Investigate: schema choice
   (per-sample columns vs `ListArray<int32>` for carriers),
   compression off, explicit `memory_map=True` flag, page-cache
   prewarm.
3. **Read throughput < 100 MB/s aggregate** → mmap contention or
   page-cache thrash under W workers. Might need explicit `pread`,
   one Arrow file per chromosome (already the design), or a
   smaller record-batch size.

### Decision tree after Spike 2

| Result | Action |
|---|---|
| All pass | Green-light **5d.1** branch. Build with confidence. |
| Streaming-write fails | Investigate Arrow writer options; if no fix, consider raw memory-mapped numpy arrays without Arrow. |
| Mmap-share fails (per-worker copies) | Re-test with explicit `memory_map=True` and compression-off; if still failing, drop Arrow and use raw numpy `np.memmap` directly with our own column-offset format. |
| Throughput too low | Profile to identify whether bottleneck is mmap read, decode, or worker-side iteration; consider one-file-per-chromosome layout or smaller batches. |

### How to run

Requires `pyarrow` to be installed in the active env (`pip install
pyarrow`). The script exits with status 3 and a clear message if
it isn't.

```bash
# default (10k samples × 100k sites, 8 workers, ~1 GB Arrow file)
python synthetic_people/scripts/spikes/spike2_arrow_streaming.py

# scale-up smoke (~10 GB Arrow file)
python synthetic_people/scripts/spikes/spike2_arrow_streaming.py \
    --samples 50000 --sites 200000

# stress test (~50 GB Arrow file — only on hosts with adequate disk)
python synthetic_people/scripts/spikes/spike2_arrow_streaming.py \
    --samples 100000 --sites 500000

# alternate path / fast NVMe
python synthetic_people/scripts/spikes/spike2_arrow_streaming.py \
    --path /mnt/nvme/spike2.arrow

# fewer workers (test slice-vs-no-slice scaling)
python synthetic_people/scripts/spikes/spike2_arrow_streaming.py \
    --workers 4
```

The test file is overwritten each run (Arrow's IPC writer creates
fresh files, unlike Spike 1 which reuses).

## Recording results

After each spike run, paste the script's stdout into a results
note alongside the script (e.g., `spike2_results_<date>.txt`)
so the evidence is preserved with the project. PR #30 (the 5d
plan update) can then cite the spikes' verdicts as the
pre-implementation gate. Spike 1 results from 2026-05-09 are
already recorded in `spike1_results_2026-05-09.txt`.
