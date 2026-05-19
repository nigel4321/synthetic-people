# synthetic_people — Performance Plan

Tracking checklist for the runtime / memory optimisations we agreed to
roll out across the `generate_people` pipeline. Each phase below is its
own pull request: implement only the items in one phase, open a PR,
land it, then move to the next phase. **Do not bundle phases.**

> **Standing rule for every PR in this plan:** when you implement
> anything in a section, also update the affected unit tests *and* the
> user-facing docs (`README.md`, `TUTORIAL.md`, `IMPLEMENTATION_PLAN.md`,
> in-script help text). New flags need help strings, new behaviours
> need a TUTORIAL recipe or note, and any change to determinism /
> seeding semantics must be called out in the docs.

---

## How to use this checklist

- Tick each box as you complete the item.
- When a phase's boxes are all ticked, open the PR for that phase and
  stop. Wait for review / merge before starting the next phase.
- If you're interrupted mid-phase, the unticked boxes plus the in-flight
  branch are enough to resume — no other context required.
- Branch naming convention: `perf/phase<N>-<short-slug>` (e.g.
  `perf/phase1-concurrency`).

## Workload assumptions (drives priority)

- **Primary scaling axis: `--n`** (cohort size). Treat per-person work
  as the dominant cost. The highest-ROI items are therefore Phase 1b
  (parallel per-person writes) and Phase 3a (numpy genotype matrix).
- **Target host: single machine.** `ProcessPoolExecutor` with the
  Linux fork start method is the right primitive — no need for a
  chunked work-queue or any inter-host coordination layer.
- **Stretch target: `n = 100 000+`.** A real-world request that the
  current in-RAM design cannot support: at `--n 30 --chromosomes 1-22
  --chr-length-mb 70` (≈ 2.3M cohort sites) the per-person fork pool
  already gets OOM-killed because cohort_sites alone consumes ~4-5 GB
  before fork. Phase 3's dense numpy matrix is one order of magnitude
  better but still tops out around `n ≈ 5 000` (a `n_sites × 2*n`
  uint8 matrix at `n=100k` is ≈ 460 GB). Reaching 100k requires
  *disk-backed intermediates* — Phases 5 and 6 below.
- **Strategic target: `n = 1 000 000`.** The realistic ceiling for
  a single-machine pipeline given 2026-era cloud hosts (1–4 TB NVMe,
  64–256 GB RAM). At this scale the in-RAM model breaks fundamentally,
  not just quantitatively: even Phase 5e's "parent holds the sites
  list, workers fork-share it" design fails because (a) the sites
  list at n=1M would be hundreds of GB of Python objects, and
  (b) CPython's reference-counting writes break fork-COW sharing
  in workers (memprof26 demonstrated this at n=3000 already — see
  Phase 5e Phase A post-mortem). The path to n=1M is therefore a
  **streaming, mmap-backed intermediate**: parent never materialises
  the cohort sites in Python at all; data flows
  `msprime.TreeSequence → streaming variant iterator → mmap'd
  Apache Arrow IPC file → workers (mmap-read their slice) → partial
  BCFs → bcftools merge`. The parent's RAM stays at O(1) per variant
  for the entire cohort phase. This is **Phase 5d** below — the
  scope expanded substantially after memprof26 surfaced that the
  n=1M ceiling is structural, not throughput-bound.
- **Storage and I/O become first-class constraints at n=1M.** A raw
  cohort BCF for `n=1M × 22 chroms × 70 Mb` is ~1.7 TB on disk
  (compressed). The intermediate Arrow file roughly doubles peak
  disk during a chromosome (drop-after-merge pattern). Total
  scratch disk needed: 200-300 GB during a single chromosome's
  processing, with ~80 GB final BCF per chromosome staying for
  the duration of the run. **Disk space and I/O bandwidth become
  the binding resources at this scale**, not RAM. NVMe at ~3 GB/s
  is comfortable; spinning disks at ~150 MB/s are 20× slower and
  dominate wall time. Cloud-instance disk choice matters more than
  CPU at n=1M.

## Strategy summary

The phases compose, not compete:

| Phase | What | Scale unlock | Primary cost |
|---|---|---|---|
| 1 ✓ | Parallel chrom + per-person | Faster, not larger | None |
| 2 ✓ | Overlap I/O loaders with simulation | Faster, not larger | None |
| 3 (partial) | Dense numpy genotype matrix | `n ≈ 5 000` | New in-RAM shape |
| 4 | Hot-path cleanups | Faster | None |
| **5a/b/c ✓** | **Disk-backed cohort BCF + sparse carriers + streamed run** | **`n = 100k+`** | **+disk I/O at phase boundaries** |
| **5f ✓** | **Chunked simulation within a chromosome** | Fits 16 GB hosts | Chunk-boundary LD |
| **5e** Phase A | Sample-slice parallel BCF write | Cohort phase ~2.7× faster | Removes parallel-chromosome path |
| **5e** Phase B | Workers walk tree directly (deferred) | Extraction parallelism | Larger refactor |
| **5f'** | Constant-term chunk-RAM calibration (open) | Auto-derate sees per-worker constant cost | None — bug fix |
| **5g.1/.2 ✓** | **Batched per-person fanout extraction** | Fanout ~10× fewer bcftools invocations | None — drop-in |
| **5g.3** | Disk-spilled fanout batch handoff (planned) | Decouples B from W RAM ceiling | +staging disk |
| **5d.1 ✓** | **Streaming-mmap cohort intermediate (Apache Arrow)** | **`n = 1M`** — parent RAM becomes O(1) per variant | New dep (`pyarrow`) + ~200-300 GB scratch disk per chrom |
| **bcftools merge ✓** | `--threads min(4, workers)` on cohort merge + index | Cohort merge wall ~1.3-2× faster at real-cohort scale | None — drop-in (PR #89) |
| **5d.2** | Direct binary BCF via pysam (optional) | Writer throughput at `n = 1M+` | New compiled dep |
| **6** | **SQLite for side-state (truth / ancestry / manifest)** | Filesystem at 100k | +schema |

Phase 3's `gts: list[str]` baseline measurement is recorded and still
useful (it confirms the underlying RAM problem is real). The dense
numpy refactor itself is **superseded by Phase 5** for the 100k goal —
once cohort state lives on disk in BCF form, the in-memory data shape
matters less. Two tasks from Phase 3 *survive* into Phase 5 and become
more useful there: the `person_records_from_cohort` generator pattern
and the `carriers[i]` pre-bucketing.

---

## Phase 1 — concurrency, low risk

**Goal:** parallelise the work that's already independent. No
data-shape changes. Visible end-to-end runtime win on multi-core hosts.

**Branch:** `perf/phase1-concurrency`

- [x] Parallelise per-chromosome msprime simulations
  - `coalescent.py:simulate_cohort` currently iterates chromosomes
    serially. Each chromosome is independent.
  - Use `concurrent.futures.ProcessPoolExecutor`,
    `max_workers=min(len(chromosomes), os.cpu_count())`.
  - Pre-derive one deterministic seed per chromosome from the master
    `--seed` *before* spawning workers — determinism must survive.
  - Add a `--workers` CLI flag (default: auto) to cap concurrency for
    constrained hosts.
- [x] Parallelise per-person VCF writes
  - The loop at `cli.py:552` is the hot path for large `--n`.
  - Use `ProcessPoolExecutor` with the default fork start method on
    Linux so `cohort_sites` is shared copy-on-write — do not pickle
    cohort_sites per task.
  - Each worker derives its own `random.Random` from a per-person seed
    so output stays bit-identical to the serial run when the same
    `--seed` is passed.
  - Reuse the `--workers` flag added above.
- [x] Stream straight into `bgzip -c`
  - Drop the plain `.vcf` intermediate in `writer.py:81-167`.
  - Open `Popen(["bgzip", "-c"], stdin=PIPE, stdout=open(out, "wb"))`,
    write records into stdin, close, then run `tabix` against the
    final `.vcf.gz`.
  - Saves one disk pass per person and removes a fork/exec.
- [x] **Tests** — extend `tests/test_writer.py` (or add a new module) to
  cover: deterministic output across `--workers=1` vs `--workers=N`
  given the same seed; correctness of the bgzip-pipe path (round-trip
  via `bcftools view`). Added in `tests/test_phase1_concurrency.py`.
- [x] **Docs** — `README.md` Performance section, `TUTORIAL.md` §9
  ("Performance and scaling"), and the `--workers` `--help` text.

> **Phase 1 caveat:** Phase 1 changed how the master rng is consumed
> (one `rng.randint` per chromosome and per person up front), so output
> at a given `--seed` differs from pre-Phase-1 runs at the same seed.
> Output remains deterministic for any choice of `--workers` once on
> the post-Phase-1 code, which is what the determinism tests verify.

---

## Phase 2 — overlap I/O with compute

**Goal:** make ClinVar / rsID / COSMIC loaders run while msprime is
still simulating, since they're I/O-bound on bcftools subprocesses.

**Branch:** `perf/phase2-overlap-loaders`

- [x] Run overlay loaders concurrently with simulation
  - `load_clinvar_index`, `load_rsid_pool`, and
    `load_cosmic_records` are all bcftools subprocess + I/O. They
    release the GIL.
  - Submitted via `submit_overlays(...)` to a 3-worker
    `ThreadPoolExecutor` *before* `simulate_cohort` runs. Blocked on
    each future at the point `cohort_sites` needs it (inside the
    existing overlay-injection block).
  - COSMIC submission is gated on `--somatic`, preserving the
    "registration-gated, never auto-fetch" guarantee.
  - `--somatic` / `--cosmic-vcf` validation moved up before the
    prefetch so a bad path fails fast (before the simulation).
- [x] **Tests** — `tests/test_phase2_prefetch.py` covers skip rules
  (default / `--rsid-density=0` / `--somatic`), loader-arg pass-
  through, and identity (resolved-future payloads match what calling
  the loaders directly returns).
- [x] **Docs** — README.md notes overlays prefetch in parallel with
  the simulation; PERFORMANCE_PLAN.md ticks the Phase 2 boxes. No
  user-visible flag change.

---

## Phase 3 — memory, medium risk

**Goal:** drop the per-site `list[str]` representation that dominates
RAM at large `n × n_sites`. Big peak-RAM drop, also speeds up
downstream loops.

**Branch:** `perf/phase3-genotype-matrix`

> **Status — skipped in favour of Phase 5.** The measure-first task
> is done and validated the underlying problem (gts share is 60–87%
> of peak RSS at `n = 200 / 500 / 1000`). The *dense numpy refactor*
> below was confirmed during the Phase 5 strategy review to **not be
> implemented** — Phase 5's disk-backed cohort BCF reaches 100k+
> directly, and a transitional in-RAM dense matrix would mean two big
> data-shape changes in a row. The `person_records_from_cohort`
> generator and `carriers[i]` pre-bucketing tasks listed below
> survive into Phase 5 where they map naturally onto the chromosome-
> streaming model; the dense-matrix tasks are kept here for the
> historical record only.

> **Measure before refactoring.** CPython interns short repetitive
> strings (`"0|0"`, `"0|1"`, etc.), so the projected ~1.25 GB number
> below is a worst case. The refactor only earns its complexity if a
> measured baseline confirms the win. The first task in this phase is
> to take that measurement; the rest of the phase only proceeds if the
> baseline justifies it.

### What we lose by switching representations

Be eyes-open about the tradeoffs before starting. The list-of-strings
representation gives us several things that the numpy uint8 path
either has to encode explicitly or drop:

- **Phasing as data, not assumption.** `"0|1"` (phased) vs `"0/1"`
  (unphased) lives in the string today. uint8 only encodes indices, so
  we'd be hardcoding "the cohort is always phased" into the data
  model. The pipeline does emit phased GTs uniformly, so this is
  mostly a documentation issue — but it means the data shape can no
  longer represent unphased calls without an extra flag. Mitigation:
  add a `cohort_phased: bool` flag at the top of the data structure
  to keep phasing explicit.
- **Per-call missingness.** `./.` for dropouts has no obvious uint8
  home — you'd reserve a sentinel (e.g. 255) or carry a parallel mask.
  Currently dodged because dropouts are layered on at write time in
  `writer.py` and never stored on cohort sites; only a problem if a
  future change wants missingness to live on `cohort_sites`.
- **Multi-allelic ceiling.** uint8 caps allele indices at 255. Trivial
  for human germline biallelic SNVs/indels; uint16 buys 65k alts at 2×
  the RAM if we ever need it.
- **Overlay-injection delicacy.** ClinVar / rsID / COSMIC injection
  (`cli.py:443-504`) currently mutates site dicts in place by index.
  With a cohort-wide matrix `(n_sites, 2*n_people)` plus parallel
  metadata, sorting / deduping / reordering becomes "two arrays moving
  together" instead of "one list of dicts." Mitigation: wrap the
  matrix + metadata in a small `CohortSites` dataclass with mutation
  methods rather than letting them float free.
- **Debugger / fixture ergonomics.** `site["gts"] == ["0|1", "1|1",
  ...]` is instantly readable in a debugger and trivially JSON-
  serialisable. Numpy arrays aren't. Every test that hand-constructs
  a site dict needs rewriting against the new shape — affects roughly
  10 test files including `test_cohort.py`, `test_overlays.py`,
  `test_truth.py`. This is dev-time cost, not runtime cost, but it's
  real.
- **Per-site heterogeneity (theoretical).** Each site dict can in
  principle carry its own length / ploidy / encoding. A cohort-wide
  matrix locks every site to exactly `2*n_people` slots. Not exploited
  today, but it removes a degree of freedom we have today.

If the baseline measurement comes back small (e.g. <30% RAM win at
`n=500`), reconsider the cost/benefit before pressing on. A cheaper
intermediate is to keep `list[str]` but call `sys.intern("0|1")` etc.
explicitly so the saving is documented, not incidental.

### Tasks

- [x] **Measure first** — establish the baseline before any refactor
  - `scripts/profile_memory.py` runs `simulate_cohort` for one `--n`
    at a time and reports peak RSS via
    `resource.getrusage(RUSAGE_SELF).ru_maxrss` plus the per-site
    `gts: list[str]` overhead measured with `sys.getsizeof` and
    string-id deduplication.
  - Numbers recorded in the **Baseline** subsection below, taken on
    `nigel-test2` (Linux x86_64, CPython 3.12.3) at
    `--chr-length-mb 5.0`.
  - **Decision gate result:** gts share is **60.0% / 79.6% / 86.7%**
    at `n=200/500/1000` — far above the 20% threshold. Refactor
    proceeds; `sys.intern` fallback not taken.
- [ ] Replace `gts: list[str]` with a numpy uint8 representation
  - Today: every site dict carries a list of `"0|1"` strings of length
    `n_people`. At `n=500, n_sites=50_000` this is ~25M Python strings
    (~50 B each → ~1.25 GB *worst case* — measure before relying on
    this number).
  - Switch to `numpy.ndarray[uint8]` of shape `(2*n_people,)` per site,
    or — preferred — a single cohort-wide matrix
    `(n_sites, 2*n_people)` that all downstream code indexes.
  - Wrap the matrix + per-site metadata in a `CohortSites` dataclass
    with explicit mutation methods (sort, dedupe, inject, annotate)
    so overlay code keeps the two arrays in lockstep.
  - Add a `cohort_phased: bool` flag on the dataclass to keep the
    phasing assumption explicit.
  - Format `"0|1"` strings only at the moment of writing, in
    `writer.py`.
  - Touches: `coalescent.py`, `cohort.py`, `admixture.py`, `writer.py`,
    plus all test files that synthesise site dicts.
- [ ] Make `person_records_from_cohort` a generator
  - Today it returns a full list per person, then the writer iterates
    it in genome order.
  - Yield records already in genome order so the writer streams into
    the bgzip pipe without holding the per-person list in RAM. Pairs
    naturally with the Phase 1 bgzip-pipe change.
- [ ] Pre-bucket per-person alt observations during cohort assembly
  - Build `carriers[i] = list_of_site_indices` so the per-person walk
    is O(records that person actually carries) instead of
    O(total cohort sites).
  - Win is linear in `n_people` and largest on singleton-heavy SFS
    (the default).
- [ ] **Re-measure** — re-run the baseline script after the refactor
  and record the realised RAM win in this file alongside the
  baseline. If the realised win is materially below the projection,
  capture why so we don't make the same mistake on a future
  optimisation.
- [ ] **Tests** — every test that constructs a site dict by hand needs
  updating; the genotype-matrix refactor is the riskiest item in this
  plan, so plan extra coverage of:
    - serial vs parallel (Phase 1 still passes deterministic check);
    - `person_records_from_cohort` generator yields the same records
      in the same order as the old list-returning version;
    - SFS histogram and overlay-stats numbers are byte-for-byte
      unchanged;
    - `cohort_phased=False` round-trips correctly (even if not used
      in production today, the flag should be exercised so it doesn't
      silently rot).
- [ ] **Docs** — note the new in-memory representation and the
  `cohort_phased` flag in `IMPLEMENTATION_PLAN.md` (architecture
  section); call out the realised (measured, not projected) RAM win
  in `README.md` Performance and `TUTORIAL.md` §9.

### Baseline (filled in by `scripts/profile_memory.py`)

Measured 2026-05-04 on `nigel-test2` (Linux 6.17.0-1011-azure x86_64,
glibc 2.39, CPython 3.12.3) at `--chr-length-mb 5.0 --chromosomes 22
--demo-model OutOfAfrica_3G09 --population CEU --seed 42`.

| n | n_sites | sim_secs | peak RSS (MB) | gts overhead (MB) | gts share | unique GT strings |
|---|--------:|---------:|--------------:|------------------:|----------:|------------------:|
| 200  | 13,610 |  2.67 |   226.3 |   135.7 | **60.0%** |  2,722,000 |
| 500  | 18,394 |  6.45 |   577.4 |   459.9 | **79.6%** |  9,197,000 |
| 1000 | 23,837 | 14.62 | 1,385.8 | 1,201.6 | **86.7%** | 23,837,000 |

Reproduce with:

```bash
.venv/bin/python synthetic_people/scripts/profile_memory.py --n 200  --chr-length-mb 5.0
.venv/bin/python synthetic_people/scripts/profile_memory.py --n 500  --chr-length-mb 5.0
.venv/bin/python synthetic_people/scripts/profile_memory.py --n 1000 --chr-length-mb 5.0
```

### Decision

**Proceed with the refactor.** Gate was: refactor only if gts overhead
≥20% of peak RSS at any of the three sample sizes. The realised
share is **60% / 79.6% / 86.7%** — three to four times the gate, and
*growing* with `n`. A finding from the measurement that should be
called out:

- **No interning.** `unique_gt_strings == n_sites × n` at all three
  cohort sizes, exactly. The `f"{a}|{b}"` formatter creates a fresh
  string object every call — CPython does not intern these. The
  worst-case "≈50 B per string" projection from the original Phase 3
  text is therefore the realised cost, not an upper bound.

The refactor described in the next subsection is therefore expected
to recover most of the gts overhead (substituting a packed
`numpy.uint8` matrix at ≈1 byte per haplotype slot for what is
currently ≈55 B per slot via PyObject + small-string overhead).
Projected post-refactor gts share: well under 5% at every cohort
size.

### After-refactor numbers (to be filled in by the re-measure task)

| n | n_sites | host | peak RSS | delta vs baseline |
|---|---------|------|----------|-------------------|
| 200  | TBD | TBD | TBD | TBD |
| 500  | TBD | TBD | TBD | TBD |
| 1000 | TBD | TBD | TBD | TBD |

---

## Phase 4 — small cleanups (do alongside Phase 1 if convenient)

**Goal:** tidy hot-path micro-issues uncovered while reading the code.
Low priority; ship as a follow-up if Phase 1 doesn't naturally absorb
them.

**Branch:** `perf/phase4-cleanups`

- [ ] Replace the `used_positions: set` collision-defence in
  `_tree_sequence_to_sites` with a single trailing dedupe pass after
  sorting. The set grows to `n_sites` entries; collisions are rare
  with biallelic msprime and integer rounding.
- [ ] Remove the per-record `"0|1"` parse step in `writer.py`
  (`alt_dosages`, `gq_from_ad`) once Phase 3 lands the numpy
  representation. They become array slices.
- [ ] **Tests** — confirm site-dedupe parity with a fixed seed; spot-
  check `alt_dosages` on the new path.
- [ ] **Docs** — none needed unless behaviour changes.

---

## Phase 5 — disk-backed cohort BCF (path to 100k+)

**Goal:** the cohort genotype state never has to fit in RAM. Simulate
chromosome-by-chromosome straight onto disk as BCF; per-person VCFs
become an opt-in derivation via `bcftools view -s SAMPLE` from the
cohort BCF. This is the architectural shift that unblocks the 100k
stretch target — the in-memory phases stop scaling around `n ≈ 5 000`
even with Phase 3's dense matrix.

The phase ships in two PRs:

- **Phase 5a — `perf/phase5-cohort-bcf`**: BCF writer module +
  `--mode {per-person, cohort, both}` flag. The cohort BCF gets
  written from the in-memory cohort_sites at the end of the overlay
  phase; `--mode cohort` skips the per-person fan-out. The streaming
  refactor and per-person-from-BCF derivation are deferred to 5b so
  this PR stays reviewable.
- **Phase 5b — `perf/phase5b-streaming` (planned)**: chromosome-by-
  chromosome streaming refactor (no in-memory cohort_sites accumulator),
  per-person derivation via `bcftools view -s`, resume contract via
  `cohort.meta.json`. Without 5b, 5a still bypasses the per-person
  fork-pool RAM amplification (`--mode cohort` skips fan-out
  entirely), which is a meaningful scale unlock on its own.

### Why BCF (not SQLite, not Parquet)

`htslib` is already a hard dep, BCF is the binary VCF spec, `bgzip +
tabix` gives indexed concurrent reads, and `bcftools view -s SAMPLE`
extracts one person in one command — that *is* the per-person
derivation step. SQLite for the same data would require designing a
sites/samples/genotypes schema that re-implements what BCF already
does at roughly 10–100× the I/O cost per query. SQLite is still the
right home for the *side-state* that isn't shaped like VCF — see
Phase 6.

### What we lose vs. an in-RAM cohort

- **Disk-I/O at every phase boundary.** Today everything stays in
  Python objects from `simulate_cohort` through to per-person writers.
  Phase 5 introduces a serialise / deserialise round-trip per
  chromosome. For sub-100 cohorts this is measurably slower; for 100k
  it's the only way the run completes at all.
- **Less freedom to mutate cohort state in-place.** Today's overlay
  injection (`cli.py:443-504`) walks the in-RAM list and rewrites
  fields by index. Once the canonical state is on disk, overlays
  must apply *before* writing each chromosome's BCF (one-pass), or
  via a separate "rewrite-the-BCF" pass. The one-pass model is the
  recommended shape — see the chromosome streaming task below.
- **Resume semantics need a contract.** If we let an existing
  `cohort.chr*.bcf` skip simulation, the seed + contig + n + demo
  model + length must all match. A small JSON sidecar
  (`cohort.meta.json`) listing those parameters lets us reject a
  cohort BCF that doesn't belong to this run.

### Phase 5a tasks (BCF deliverable + `--mode` flag)

- [x] **Pick the on-disk layout.** Resolved: single
  `out/cohort/cohort.bcf` with a CSI index. Per-chromosome split
  deferred to 5b once the streaming refactor lands and per-chromosome
  files become a natural unit; for 5a a single file is simpler and
  matches what `--mode cohort` produces in one pass.
- [x] **Add a BCF writer.** Resolved: `subprocess.Popen(["bcftools",
  "view", "-O", "b", "-"], stdin=PIPE)`. pysam isn't on the
  dependency tree — msprime / tskit / stdpopsim don't pull it in,
  and adding it is a meaningful install-time cost (its own htslib
  build with C extensions). `syntheticgen/bcf_writer.py` follows the
  same Popen pattern `writer.py` uses for `bgzip -c`.
- [x] **`--mode {per-person, cohort, both}` flag.** Default
  `per-person` (zero behaviour change for existing users). `cohort`
  writes the BCF and skips per-person fan-out; `both` writes both
  deliverables.
- [x] **Progress logging on long phases.** Throttled (~5 s cadence)
  heartbeat lines during the cohort BCF write loop and the per-person
  fan-out, so a multi-hour 100k-sample run has visible progress.
- [x] **Tests** — `tests/test_bcf_writer.py` (BCF round-trip,
  per-sample extraction, header shape, arg validation) and
  `tests/test_cli_modes.py` (each of the three modes lands the right
  artefacts and manifest fields). Both gate on bcftools/tabix/bgzip
  on PATH; the cli-mode tests additionally gate on msprime + stdpopsim.
- [x] **Docs** — README §Performance describes the flag and the
  manifest-output table; TUTORIAL §9.3 walks through the three modes
  and the `bcftools view -s` per-person derivation pattern; §9.4
  documents the progress-logging cadence.

### Phase 5a-completed manifest task

- [x] **Manifest extension.** `shape` (`per-person` / `cohort` /
  `both`) records which `--mode` produced the run; top-level
  `samples[]` always present so callers don't need a per-mode path
  for "list every sample"; `cohort_bcfs` is a list (singleton in 5a;
  populated with per-chromosome paths once 5b's streaming refactor
  lands). Per-person `people[]` entries only emitted when the
  per-person fan-out ran.

### Phase 5b1 tasks (streaming `--mode cohort` only)

- [x] **Stream chromosome-by-chromosome in `simulate_cohort`.** New
  `simulate_cohort_iter` generator yields one chunk per chromosome;
  the legacy flat-list `simulate_cohort` is now a thin wrapper for
  callers (admixture path, fixture builders) that need the full
  cohort materialised. Per-chromosome seeds are still pre-derived
  before any worker fan-out, so determinism within the streamed path
  is preserved at fixed `--seed`.
- [x] **Per-chromosome BCF layout** when streaming —
  `out/cohort/cohort.chr<N>.bcf` files. Pairs naturally with 1000G's
  per-chromosome shape. Manifest's `cohort_bcfs[]` (list shape from
  5a) populated with one entry per chromosome.
- [x] **`--mode cohort` only** in 5b1. Per-person and both modes
  keep today's in-memory path, unified in 5b2.

### Phase 5b1 measurements

Quick comparison at `n=500` × `chromosomes=20,21,22` × `--chr-length-mb=5`
(no overlays / no SVs / no errors), Linux x86_64, CPython 3.12.3:

| `--mode` | Peak RSS | Wall time | Notes |
|---|--------:|--------:|---|
| `per-person` | 1.9 GB | 2:27 | in-memory cohort_sites + per-person fan-out |
| `cohort` (streamed) | 0.85 GB | 0:19 | per-chrom BCFs, no per-person fan-out |

The peak-RSS halving comes from never materialising the whole cohort
in RAM at once; the wall-time speedup is partly that and partly
skipping the per-person fan-out entirely (cohort mode doesn't write
per-person VCFs by design — derive them later via `bcftools view -s`).
A full benchmark at `n = 1 000 / 10 000 / 100 000` lands alongside
Phase 5b2 once per-person derivation closes the loop.

### Phase 5b2 tasks (per-person derivation + resume) — completed
- [x] **Per-person derivation from cohort BCF.** New
  `syntheticgen/cohort_derivation.py` runs a two-step
  `bcftools view -s SAMPLE | bcftools view -e 'GT="ref"'` pipeline
  against each per-chrom cohort BCF and parses the records into the
  same dict shape `person_records_from_cohort` returns, so
  `write_person_vcf` consumes either source identically. The
  pipelined form is necessary because `bcftools view -s SAMPLE -e
  'GT="ref"'` evaluates the filter against the multi-sample GT
  before sample subset is applied — fixed in 5b2 by chaining two
  views.
- [x] **Streamed pipeline now drives all three modes.** The 5b1
  diversion check widened: any `--mode` value on the non-legacy
  non-admixture coalescent path now goes through
  `_run_cohort_streamed`, which after streaming branches on mode —
  cohort mode early-returns (today's 5b1 behaviour), per-person and
  both run a fan-out fed by `derive_person_records`.
- [x] **Resume contract via `cohort.meta.json`.** New
  `syntheticgen/resume.py` persists the resume-identity-relevant
  params (seed, n, build, chromosomes, chr_length_mb, demo_model,
  population, rec_rate, mu) plus the sample IDs, per-person seeds,
  and per-chromosome overlay seeds drawn at run start. Each
  chromosome's completion appends to `completed_chromosomes` via an
  atomic-rename write, so a SIGINT mid-flush leaves the prior
  version intact. On startup, a matching meta.json skips already-
  complete chromosomes; mismatched params raise `ResumeMismatch`
  with a clear `--no-resume` hint.
- [x] **Per-chrom overlay seeds.** Each chromosome's overlay rng is
  now seeded from `resume.overlay_seeds[chrom]` rather than the
  master rng. Required for resume to be deterministic — without it,
  the master rng's state at chrom_K depends on K's predecessors, so
  a resume after chrom K-1 would consume rng differently than a
  non-interrupted run did.
- [ ] **Full re-measure at large n.** A focused benchmark at
  `n = 1 000 / 10 000 / 100 000` with the streamed pipeline + the
  new BCF-based per-person derivation, measuring peak RSS and total
  wall time. Lands as a follow-up — the 5b2 PR ships the path
  itself; we don't gate the merge on a multi-hour benchmark run.
- [x] **Tests.**
  - Cohort-derivation parity:
    `tests/test_cohort_derivation.py` round-trips a cohort through
    the BCF writer + derivation module and asserts per-record
    equivalence vs `person_records_from_cohort`.
  - End-to-end: `tests/test_cli_modes.py` updated for the new
    streamed-with-derivation behaviour under all three `--mode`
    values; cohort BCFs land as intermediates under per-person too.
  - Resume: `tests/test_resume.py` covers the fresh-start, matching-
    params reuse, mismatched-params error, `--no-resume` wipe paths,
    plus an end-to-end test that deletes one chromosome's BCF mid-
    run and re-runs to confirm the surviving chrom's BCF is
    untouched (mtime preserved) while the deleted one regenerates.
- [x] **Docs.** README Performance / Output layout updated for the
  resume contract and `--no-resume` flag; TUTORIAL §9.3 walks
  through resuming an interrupted run.

### Open questions for review

- Should overlays still inject after-simulation, or fold into the
  per-chromosome streaming loop? Folding in is cleaner but means the
  ClinVar / rsID / COSMIC pools must all be loaded before the first
  chromosome simulates (they already are, via Phase 2's prefetch).
- Per-chromosome BCFs vs. one-big-BCF: are there downstream tools that
  prefer one shape? 1000G ships per-chromosome; gnomAD ships per-
  chromosome. Recommend per-chromosome.
- ~~Threshold for default `--mode`~~ resolved: default is
  `per-person` with no threshold; large-cohort users opt in via
  `--mode cohort`.

---

## Phase 5c — sparse in-memory genotype storage

**Goal:** drop the per-chromosome RAM ceiling that 5b's streaming
inherited. Phase 5b bounded peak at one chromosome's working set,
but each chromosome's working set is `n × n_sites × ~100 B` for the
dense `gts: list[str]` representation — at `n=3 000 × chr1×70 Mb`
that's already ~30 GB, OOMs a 32 GB host. 5c stores cohort-site
genotypes sparsely (carriers list of `(haplotype_idx, allele_idx)`
for non-zero entries only), so per-chromosome RAM scales with alt
observations rather than `n × n_sites`. SFS is singleton-dominated,
so total alt observations grow as `n_sites × log(n)` — flat enough
that n=100k+ fits comfortably on a workstation.

**Branch:** `perf/phase5c-sparse-carriers`

### Memory model

Per chromosome at full chr1 (~70 Mb, ~100 k cohort sites):

| `--n` | Dense `list[str]` (5b) | Sparse `carriers` (5c) |
|---|--------:|--------:|
| 500 | ~5 GB | ~10 MB |
| 3 000 | ~30 GB (OOMs 32 GB host) | ~50 MB |
| 100 000 | ~1 TB (infeasible) | ~600 MB |
| 1 000 000 | ~10 TB (infeasible) | ~6 GB |

Numbers reflect cohort-sites RAM only; tskit's own tree-sequence
overhead is independent and adds ~1-2 GB at 1M scale. Sparse storage
reaches n=100k cleanly and gets within reach of n=1M; the next
bottleneck at 1M is the writer side — see Phase 5d below.

### Tasks

- [ ] **Add `cohort_sites.py` helper module** with
  `carriers_from_dense_gts()`, `dense_gts_from_carriers()`,
  `gt_for_person()`. Tests can keep building site fixtures from
  dense GT lists; the helpers convert at the boundary.
- [ ] **Refactor `coalescent._tree_sequence_to_sites`** to emit
  sparse carriers from `var.genotypes` instead of dense
  `list[str]`. Per-record memory drops from `O(n_people)` strings to
  `O(non-zero entries)` int tuples.
- [ ] **Refactor `admixture` simulator** the same way (it walks
  tree sequences identically).
- [ ] **Refactor `cohort.draw_cohort_background`** (legacy 1000G
  path) to emit carriers from the slot array `assign_haplotypes`
  already produces. Drops `_gts_from_slots`.
- [ ] **Refactor `cohort.person_records_from_cohort`** to derive
  this person's GT by scanning the carriers list rather than
  indexing into a dense `gts` list. Per-site cost is O(carriers in
  this site) — typically a small constant under SFS-realistic
  cohorts.
- [ ] **Refactor `bcf_writer.CohortBcfWriter.write_site`** to expand
  carriers to dense GT strings at write time. Keep accepting
  `site["gts"]` as a fallback for incremental fixture migration.
- [ ] **Tests.** Update fixtures that hand-build cohort sites to
  go through the helper. Add a regression test that runs
  `simulate_cohort_iter` at moderately-large `n` and asserts
  per-chromosome RAM stays under a sparse bound.
- [x] **Quick measurement** at `n = 500 × chr22 × 70 Mb` (the
  full-chromosome case the user originally failed on at n=30):
  peak RSS **2.4 GB**, wall **1:19**. Cohort-sites RAM has dropped
  from ~9 GB pre-5c (dense `list[str]`) to ~50 MB (sparse
  carriers); the remaining 2 GB is tskit's tree sequence + ClinVar
  load (independent of the cohort-sites refactor).
- [ ] **Full benchmark at `n = 1 000 / 10 000 / 100 000`** —
  follow-up. The 5c PR ships the path itself; a multi-hour
  benchmark sweep lands separately.
- [x] **Docs** — README + TUTORIAL: the sparse refactor is
  internal, users see the same CLI and the same output. Plan
  table shows projected vs realised peak RSS.

---

## Phase 5f — chunked simulation within a chromosome

**Status:** planned, not started. **Highest priority of the
remaining 5 series** — directly addresses the user-reported OOM
that surfaced post-5c, and Phase 5e's parallel-extraction model
isn't useful until each tree sequence fits.

**Goal:** bound msprime's per-chromosome working memory by splitting
each chromosome into sub-chunks of fixed length (e.g. 5-10 Mb),
each simulated as an independent tree sequence. Targets the OOM
the user hit *during the simulation itself*, not the extraction
phase.

**Branch (when implemented):** `perf/phase5f-chunked-simulation`

### Why even one tree sequence doesn't fit

Phase 5c took cohort-sites RAM out of the picture; the user's
n=3000 × full chr1 70Mb run *now OOMs in msprime's working memory*
during simulation, before extraction even starts. The output trace
from the failing run:

```
  simulating chrom 1 (length 70.0 Mb, model=OutOfAfrica_3G09)...
Killed
```

stdpopsim's `OutOfAfrica_3G09` simulates three populations with
bottlenecks, migrations, and an out-of-Africa expansion. During
the coalescent backward-walk, the **active-lineage count peaks at
many multiples of n** through bottlenecks (likely 10-50× at
n=3000), and each lineage carries segment metadata that
recombination splits further. Combined with full chr1 length
(most recombination, biggest tree count), msprime's working memory
can be 8-16+ GB at n=3000 × 70Mb. That exceeds the budget on
workstation-class hosts even before Phase 5e gets a chance to
parallelise the extraction.

### Proposed structure

1. Split each chromosome into K sub-chunks of fixed length
   (e.g. `--chr-chunk-mb 10` produces 7 chunks for a 70 Mb chrom).
2. For each chunk, simulate independently using msprime/stdpopsim
   with the contig sliced to that chunk's range. Each chunk's
   tree sequence is roughly 1/K the size of the full-chromosome
   version.
3. Iterate variants out of each chunk's tree sequence into the
   per-chrom BCF in genome order — chunks are processed in order,
   positions adjusted by the chunk's start offset.
4. Free each chunk's tree sequence before starting the next.

### Memory model

| Path | Per-chrom peak | Status at n=3000 × 70Mb |
|---|---|---|
| Pre-5f, full-chrom simulation | 8-16+ GB | OOMs workstation-class hosts |
| 5f at `--chr-chunk-mb 10` (7 × 10Mb) | ~1-2 GB per chunk | Fits 16 GB host |
| 5f at `--chr-chunk-mb 5` (14 × 5Mb) | ~500 MB - 1 GB per chunk | Fits 8 GB host |

Total per-chrom wall time: roughly the same total work, with
overhead from K independent simulations. Each chunk pays msprime's
startup cost; for 5-10 Mb chunks at n=3000 that overhead is small
compared to the simulation itself.

### What we lose: cross-chunk LD

The big tradeoff. Recombination events that would have spanned
chunk boundaries are lost — chunks simulate independently, so
haplotypes are uncorrelated across chunk boundaries. In genetics
terms: linkage disequilibrium decays sharply at chunk boundaries.

Effect on common analyses:

| Analysis | Impact |
|---|---|
| Ti/Tv | unaffected |
| Allele frequency spectrum | unaffected |
| Per-person genotype lists | unaffected (each chunk's haplotypes stay consistent across the cohort) |
| ClinVar / dbSNP overlay placement | unaffected (overlays operate on positions within chunks) |
| Short-range LD (≤ chunk_length) | preserved within chunks |
| Long-range LD (> chunk_length) | NOT realistic — analyses requiring chr-scale haplotype block structure should not use chunked mode |

Documented as a 5f-specific caveat alongside Phase 1's
rng-consumption note.

### Tasks

- [x] Add `--chr-chunk-mb N` CLI flag. Default `0` = auto-pick
  from available RAM at run start; explicit `N > 0` overrides.
  Chunk-size selection logged so the user sees what was picked.
  `psutil>=5.9` added to `synthetic_people/requirements.txt`.
- [x] **Auto-pick logic** in `cli.py` calling
  `coalescent.auto_pick_chunk_size_mb(n, length, demo_model,
  available_bytes, workers)`. Estimate calibrated against the
  user's failing run: ~80 KiB per (sample × Mb) for OOA-class
  demography, ~16 KiB for constant-Ne. Auto-pick targets ≤ 50%
  of `psutil.virtual_memory().available / workers` and caps at
  the configured `--chr-length-mb`.
- [x] **Chunk overlap** at 10% of chunk size, clamped to
  `[0.5 Mb, 5 Mb]`. Each chunk simulates
  `chunk_size + overlap_margin` bp; variants past `chunk_size`
  in chunk-local coordinates are dropped at write time so the
  per-chrom site list stays duplicate-free. Documented as
  boundary smoothing rather than true cross-chunk LD recovery.
- [x] **Refactored `coalescent.simulate_chromosome`** to
  dispatch chunked vs single via `_simulate_chromosome_chunked`.
  Each sub-chunk runs `_simulate_one` (the shared msprime
  invocation) against `species.get_contig(chrom,
  right=chunk_size + overlap)` with a chunk-specific seed. The
  admixture path stays on the single-pass simulator for now —
  its per-person ancestry segments interact with the tree-
  sequence walk in ways that need a separate refactor.
- [x] **Per-chunk seeds** derived deterministically via a
  Knuth multiplicative mix:
  `(chrom_seed + chunk_index * 0x9E3779B9) & 0x7FFFFFFF`. Avoids
  rng-state-dependence so a resumed run sees the same chunk
  seeds regardless of which chunks were already simulated.
- [x] **Free per-chunk tree sequence** before the next chunk's
  simulation starts. `del ts` after each chunk's variant
  iteration; per-chrom site list accumulates the variants.
- [x] **Tests** in `tests/test_chunked_simulation.py` (18 cases):
  - RAM estimator linearity in n and chunk_mb, OOA vs constant-Ne
    rate selection, "none" string handled.
  - Auto-pick correctness (full-fits returns length, doesn't-fit
    returns smaller, workers divide budget, constant-Ne picks
    larger chunks than OOA).
  - Chunked output: positions sorted, unique, in range; record
    count within stochastic noise of unchunked; deterministic at
    fixed seed; no duplicates at chunk boundaries.
  - Overlap-bp clamping at floor and ceiling.
- [x] **Docs.** `README.md` Performance section + CLI reference
  table updated; `TUTORIAL.md` §9.5 added with a recipe and the
  cross-chunk LD caveat table.

### Resolved decisions

1. **Default chunk size: auto-pick based on host hardware.**
   `--chr-chunk-mb 0` (the default) detects available RAM at
   startup via `psutil.virtual_memory().available`, estimates
   per-chunk working memory from `(n_people, demo_model,
   chunk_size)`, and picks the largest chunk size that fits the
   available budget with comfortable margin (target: peak working
   set ≤ 50% of available RAM, leaving room for parent process,
   overlay loaders, ClinVar pool, OS cache, and the
   simulation-startup spike). `--chr-chunk-mb N` (with N > 0) is
   the explicit override for users who want to pin a specific
   chunk size — useful for reproducibility across heterogeneous
   hosts, or to force a smaller chunk for safety. Chunk-size
   selection is logged so the user sees what was picked and why
   (`auto-selected chunk size 5 Mb based on 14 GB available RAM
   at n=3000`).
2. **Chunk overlap: yes, implemented.** Adjacent chunks simulate
   with a configurable overlap margin (default ~5-10% of chunk
   size). The overlap regions are discarded at write time —
   variants in the overlap of chunk K and chunk K+1 are written
   only once, taken from chunk K. This doesn't fully recover
   cross-chunk LD (each chunk is still an independent simulation),
   but it gives the central region of each chunk a less abrupt
   boundary effect: LD decays naturally inside each chunk's
   simulated region rather than terminating sharply at the chunk
   boundary. Documented in 5f's caveat as "boundary smoothing"
   rather than "true LD recovery". Full cross-chunk LD recovery
   would require msprime's tree-sequence continuation API and is
   left to a future phase.
3. **stdpopsim contig slicing: use `right=`.** Each chunk
   simulates an *independent* small contig (length =
   chunk_size + overlap_margin) using its own seed; chunk K's
   variants then have positions offset by `K × chunk_size`
   when written to the per-chrom BCF. We don't try to "extract
   chunk K of chr1" using `right=chunk_end` (which would re-
   simulate the prefix at every chunk and balloon work
   quadratically). Instead each chunk is biologically equivalent
   to "the first chunk_size bp of chr1", with chunk-specific
   seeds making them independent. For the user-facing analyses
   we care about (Ti/Tv, AF spectrum, per-person genotypes,
   ClinVar overlay placement, short-range LD inside each chunk)
   this is identical; for long-range LD across chunks it isn't,
   which is the documented caveat.
4. **Sizing: detect host RAM, no hardcoded numbers.** The
   auto-pick formula uses `psutil.virtual_memory().available`
   at run start. Plan and TUTORIAL document the formula as
   "chunk size scales linearly with available RAM and inversely
   with cohort size; auto-picked chunks aim for peak working set
   ≤ 50% of free RAM". Users can override with
   `--chr-chunk-mb N` if their environment has surprising memory
   pressure (e.g. running alongside other big processes).

### Relationship to Phase 5e

5f and 5e compose. **5f reduces the size of each tree sequence**
so it fits in RAM at all. **5e parallelises extraction** across
workers consuming sample slices of that smaller tree sequence.

| Host RAM | Chunk size | Workers (5e) | What gets unlocked |
|---|---|---|---|
| 16 GB | 5 Mb chunks | 2-4 sample-slice workers | n=3000 completes |
| 32 GB | 10 Mb chunks | 4-8 sample-slice workers | n=3000+ at higher throughput |
| 64+ GB | 20 Mb chunks | 8-16 sample-slice workers | n=10k+ feasible |

**Implementation order: 5f first** (unblocks the immediate user
case), then 5e (squeezes within-chunk extraction throughput on
top). 5e standalone helps only when one full-chromosome tree
sequence fits in RAM — for the user's current workload it
doesn't.

### Phase 5f post-mortem — RAM-model recalibration after the user retest

The first 5f deployment auto-picked an 8.7 Mb chunk size for
``--n 3000 --chromosomes 1-22 --chr-length-mb 70`` on a 16 GB
host with auto ``--workers 4``. The user's
``--profile-memory`` trace showed:

- Children RSS climbing linearly from 0 to **~16 GB** over 100 s
  (4 workers × ~4 GB tree sequence each at the picked chunk
  size).
- The host's RAM ceiling hit at t≈120 s, kernel started
  swap-thrashing.
- **No ``chrom X sites yielded`` mark fired** in the entire
  48-minute run — workers stalled before completing any
  chromosome.

The first 5f calibration of `80 KiB/(sample × Mb)` came from a
single full-chromosome OOM observation; the new ratio'd
measurement at a known chunk size shows the actual cost is
**~153 KiB/(sample × Mb)** at OOA scale — almost exactly 2×
under-estimate. Three follow-ups landed:

- **Coefficient: 80 → 160 KiB/(sample × Mb)** for OOA-class
  demography. Pessimistically rounded up from 153. Constant-Ne
  scales proportionally (16 → 32 KiB).
- **Auto-pick safety target: 50% → 25%** of available RAM.
  Halves the chance that residual model error or unbudgeted
  overhead (parent process, ClinVar pool, bcftools subprocesses)
  pushes total over the host ceiling.
- **`auto_derate_workers` helper.** When the auto-picked chunk
  size would drop below ~2 Mb at the requested worker count,
  reduce workers instead. Below 2 Mb the per-chunk msprime
  startup cost dominates per-chunk simulation cost; better to
  trade parallelism for chunk size when RAM is the bound, not
  CPU. Only fires when ``--workers 0`` (auto) — explicit
  ``--workers N`` is honoured.

After the recalibration, the user's failing config at
``--workers 1`` should pick a ~10 Mb chunk size with peak ≈ 4 GB,
fitting comfortably in 16 GB. With auto ``--workers``, the
derate caps parallelism so the multiplied total stays in budget.

The post-mortem is also a reminder that the calibration constants
should be re-validated whenever the demographic model catalogue
changes: a heavier-than-OOA model would push the constant up;
extreme-cohort runs (n>10k) might need their own coefficient if
the lineage-tracking overhead scales super-linearly with n.

### Phase 5f' — calibration update: add a constant term to the chunk RAM model

**Status (2026-05-08):** open. Diagnosed empirically from a user
``--workers 5`` run that OOM'd during the cohort phase on a 32 GB
host (memprof24.tsv): per-worker peak ~5.6 GB even at the
auto-picked ~1.7 Mb chunk size, indicating a per-worker
*constant* cost (~3.6 GB at OOA n=3000) that doesn't shrink with
chunk size. The current `estimate_chunk_ram_bytes` is purely
linear (`rate × n × chunk_mb`), so `auto_pick_chunk_size_mb`
underestimates total cost and `auto_derate_workers` doesn't see
that 5 workers × 3.6 GB constant = 18 GB before any chunk-linear
contribution. Result: 5 workers × ~5.6 GB ≈ 28 GB peak → OOM-kill
within the first 150 s of cohort sim.

The post-mortem above flagged exactly this risk
("extreme-cohort runs might need their own coefficient if the
lineage-tracking overhead scales super-linearly"); 5f' is the
follow-up.

**Workaround until 5f' lands:** the user's manual two-phase
invocation works — `--workers 1` for the cohort phase, then
resume + `--workers 8` for the fanout. 5f' just automates this.

**Phase 5e Phase A makes 5f' largely moot for the cohort-write
path** (because `--workers` no longer multiplies tree-sequence
RAM there). 5f' still matters for the chunked-simulation path
within a single chromosome at very high `n` or with heavier
demographic models, where per-chunk peak can exceed budget even
serially. So 5f' is lower priority now than it would have been
pre-5e — but it's the correct fix for the model and worth
shipping if/when an extreme-`n` user surfaces it again.

**Branch (when implementation starts):**
`perf/phase5f-prime-constant-term-calibration`

#### Tasks

- [ ] **Add constant-per-sample term to `estimate_chunk_ram_bytes`:**

  ```python
  CHUNK_RAM_BYTES_BASE_PER_SAMPLE_OOA = 1.2 * 1024 * 1024
  CHUNK_RAM_BYTES_BASE_PER_SAMPLE_CONSTANT_NE = 200 * 1024

  def estimate_chunk_ram_bytes(n, chunk_mb, demo):
      if demo is None or str(demo).lower() == "none":
          base, rate = (CHUNK_RAM_BYTES_BASE_PER_SAMPLE_CONSTANT_NE,
                        CHUNK_RAM_BYTES_PER_SAMPLE_PER_MB_CONSTANT_NE)
      else:
          base, rate = (CHUNK_RAM_BYTES_BASE_PER_SAMPLE_OOA,
                        CHUNK_RAM_BYTES_PER_SAMPLE_PER_MB_OOA)
      return int(base * n + rate * n * chunk_mb)
  ```

  Per-sample constants come from the empirical 3.6 GB / 3000 ≈
  1.2 MiB ratio at OOA scale; constant-Ne is the ~5× cheaper case
  the existing rate ratio implies.
- [ ] **Update `auto_pick_chunk_size_mb`:** the linear-only
  `factor = per_worker_target / full_estimate` math no longer
  works once the constant dominates. Solve for chunk_mb where
  `base + linear × chunk_mb = budget`. If `budget < base`, the
  return value should signal "this n + this many workers cannot
  fit at any chunk size" — i.e. `auto_derate_workers` must reduce
  W before the auto-pick can succeed.
- [ ] **Update `auto_derate_workers`:** include the constant term
  in the per-worker cost check. Walk down from requested workers
  until `base × n + rate × n × floor_chunk_mb ≤ available × target / W`.
- [ ] **Tests** in `tests/test_chunked_simulation.py`:
  - At n=3000, OOA, 24 GB available: derate from W=8 → W=1 (or 2,
    depending on rounding).
  - At n=500, OOA, 24 GB available: no derate, pick large chunk.
  - At n=100k, OOA: explicit budget-too-small error message,
    pointing user at `--workers 1` or a larger host.
- [ ] **Docs:** post-mortem subsection above gets a "resolved"
  status update; `--workers` `--help` text mentions auto-derate
  considers per-sample constant cost.

---

## Phase 5e — within-chromosome parallel extraction over a shared tree sequence

**Status (2026-05-08):** Phase 5f has landed (PRs #19, #21, #23). 5e
is unblocked. **Scope decided**: split into Phase A (sample-slice
BCF write only — small structural change, captures most of the
practical wall-time win) and Phase B (workers walk the tree
sequence directly — bigger refactor, deferred). Phase A is the
work to start with; Phase B is opportunistic follow-up if msprime
sim parallelism turns out to be worth the complexity. **Branch in
flight (Phase A):** TBD when implementation starts.

The motivating bug — the user's n=3000 × full-chr1 OOM during
simulation of a single chromosome — was already addressed by 5f
(chunked simulation). 5e's parallel-extraction model is the next
step: *now* one tree sequence fits in RAM, parallelising extraction
across that smaller tree sequence is the structural fix that
removes the `workers × tree_sequence_size` RAM ceiling.

**Goal:** bound peak RAM at one msprime tree sequence per run
regardless of `--workers`. After Phase 5c made cohort-sites RAM
negligible, msprime's tree sequence is the central remaining RAM
cost at any non-trivial cohort size — and today's
"one-worker-per-chromosome" parallelism multiplies it by the worker
count, OOMing on workstation-class hardware at cohort sizes msprime
itself could simulate fine.

**Branch (when implemented):** `perf/phase5e-shared-tree-extraction`

### Why each worker holds its own tree today

The parallelism is sliced *across chromosomes*, not within one.
Worker 1 builds chr1's full tree sequence (all n samples, full
chromosome length); worker 2 builds chr2's full tree (also n
samples). Each worker is doing different work — there's no
redundancy — but each holds its own multi-GB tree sequence
concurrently, so peak system RAM is:

```
peak system RAM ≈ workers × tree_sequence_size + ~1-2 GB process overhead
```

For `n=3000 × 70 Mb`, tree-sequence size is ~3-5 GB; four parallel
workers OOM a 16 GB host. This was the user-reported failure at
`--n 3000 --chromosomes 1-22 --chr-length-mb 70` after Phase 5c
landed.

### Proposed structure

1. Parent process simulates **one** chromosome's tree sequence (n
   samples, full chromosome length) — msprime is single-threaded
   internally, so running it serially in the parent matches today's
   per-chromosome cost.
2. Parent forks N workers via Linux's copy-on-write semantics.
   Workers inherit the parent's tree sequence in shared memory —
   no per-worker copy.
3. Each worker handles a different sample slice (samples
   `0..n/W`, `n/W..2n/W`, …), iterating variants for its slice and
   writing its portion of the chromosome's records to a per-worker
   intermediate.
4. After workers finish, merge their outputs into the per-chrom BCF
   in genome order (`bcftools concat` or comparable join).
5. Parent frees the tree sequence; loop to next chromosome.

### Why it works

tskit's tree sequence is mostly C-level allocations — raw tables
(`NodeTable`, `EdgeTable`, `MutationTable`, `SiteTable`) and numpy
buffers. Fork-COW shares these cleanly because workers read but
don't write. Per-worker incremental RAM is the small Python
overhead plus per-variant wrapper objects (which dirty COW pages
but are bounded). New memory model:

```
peak system RAM ≈ 1 × tree_sequence_size
                + workers × ~100 MB
                + ~1-2 GB process overhead
```

At `n=3000 × 70 Mb` that's ~3-5 GB regardless of `--workers`, vs
~12-20 GB pre-5e.

### What we lose

Simulation goes serial across chromosomes — only one tree sequence
is ever being built at a time. Today's parallel-chromosome model
finishes faster *when RAM allows it* because msprime simulation is
the dominant per-chromosome cost (~5× the extraction phase). Wall
time tradeoff:

| Path | Wall time | RAM |
|---|---|---|
| Pre-5e parallel chromosomes (when RAM allows) | fastest | OOMs at n=3000+ on 16-32 GB hosts |
| Phase 5e (this) | ~3-5× slower than parallel-chromosome ideal | bounded at 1 × tree sequence |
| `--workers 1` (current OOM workaround) | ~22× slower than parallel-chromosome ideal | 1 × tree sequence |

5e is significantly faster than `--workers 1` (parallel extraction
inside one chromosome at a time) and bounded-RAM, at the cost of
being slower than parallel-chromosome when the host has the RAM.
There's no good way to recover the parallel-chromosome speed
without paying its RAM cost; that path stays available behind a
flag if a user has the RAM and wants the wall-time win.

### Phase A — sample-slice BCF write (scoped, ready to ship)

**Goal:** parallelise the per-chromosome cohort BCF write across W
workers, each writing a contiguous sample-slice partial BCF. Parent
process keeps the existing serial flow up to the BCF write; only
that final write step gets sample-sliced. Approximate diff size:
300–400 LOC + tests.

**Why this scope (vs. parallelising msprime sim too):** the BCF
write is the dominant per-chrom cost in measured runs (~962 s /
~75% of the per-chrom wall vs ~308 s simulation+extraction at
n=3000 × 70 Mb on a 32 GB host). The write phase is dominated by
Python text formatting + bgzip/bcftools encode — both parallelise
cleanly across sample slices with simple fork-COW sharing of the
parent's sites list. msprime's coalescent + mutation simulation is
single-threaded internally and would need workers walking the tree
sequence themselves to parallelise — that's the deeper refactor
deferred to Phase B.

**Architecture (per chromosome, all serial in parent):**

1. Parent simulates tree sequence (existing 5f code, unchanged).
2. Parent walks tree → sites list with sparse `carriers` (existing
   `_tree_sequence_to_sites`, unchanged).
3. Parent applies overlays + sorts (existing `cli.py` code).
4. **NEW:** parent forks W workers via fork-COW. Each worker gets
   a contiguous sample slice `[i × n/W, (i+1) × n/W)`. Sites list
   is shared read-only via COW, no per-worker copy.
5. **NEW:** each worker writes its slice's partial cohort BCF using
   a sample-slice variant of `CohortBcfWriter` (formats only its
   slice's columns of the per-site sample block).
6. **NEW:** parent runs `bcftools merge -O b` to combine partials
   into the final per-chrom cohort BCF, then `bcftools index`.
7. Parent frees sites list. Loop to next chromosome.

**Removed in this PR:** the `ProcessPoolExecutor` parallel-chromosome
path in `simulate_cohort_iter` (the path that OOMs the user with
`--workers 5`). Cohort simulation becomes serial across
chromosomes; `--workers` controls within-chromosome BCF-write
parallelism. Users who had enough RAM for parallel-chromosome sim
get a follow-up `--parallel-chromosomes` opt-in flag (deferred —
see Resolved decisions below).

**RAM bound after Phase A:**

```
peak RAM ≈ 1 × tree_sequence_size + 1 × sites_list
        + W × ~50 MB (per-worker bgzip + temp buffers)
        + ~1-2 GB process overhead
```

At `n=3000 × 70 Mb`: ~17 GB tree + 1 GB sites + 0.4 GB workers ≈
19 GB regardless of `--workers`. Fits 32 GB host comfortably; no
auto-derate needed.

**Wall-time estimate after Phase A** (n=3000, --workers 8):
- msprime sim per chrom: ~200 s (serial, unchanged)
- tree-walk to sites: ~100 s (serial, unchanged)
- Parallel BCF write: ~960 s ÷ 8 + merge ~30 s ≈ 150 s
- Per-chrom total: ~450 s vs ~1270 s pre-5e
- 22 chroms: **~2.7 h vs ~7.4 h pre-5e — ~2.7× speedup**

Not the ~3.5× speedup of 5e Full (which would also parallelise the
~308 s extraction); Phase B captures that remainder.

#### Phase A tasks

- [ ] **Branch:** `perf/phase5e-sample-slice-bcf-write` off main.
- [ ] **`bcf_writer.py`:** add a sample-slice mode to
  `CohortBcfWriter` so a worker can format only its
  `[lo, hi)` person range's GT block per site. Sparse-carrier
  helpers in `cohort_sites.py` need a slice-aware sibling
  (`dense_gts_from_carriers_slice(carriers, slice_lo, slice_hi)`).
- [ ] **`cli.py`:** replace the existing per-chrom
  `with CohortBcfWriter(chrom_bcf, ...) as bw: bw.write_sites(...)`
  block with a new helper that:
  - Forks W workers (fork mp_context); each gets `(slice_lo,
    slice_hi)` and writes a partial BCF to `cohort/.partials/
    cohort.chr<N>.slice<i>.bcf`.
  - After all partials done, runs `bcftools merge -O b -o
    cohort/cohort.chr<N>.bcf cohort/.partials/cohort.chr<N>.slice*.bcf`
    and indexes.
  - Cleans up the `.partials/` dir for that chrom.
  - Resume: `cohort.meta.json`'s completed-chromosomes list still
    drives the skip; partial BCFs only live mid-chromosome and are
    cleaned on success.
- [ ] **`coalescent.simulate_cohort_iter`:** drop the
  `ProcessPoolExecutor` parallel-chromosome path. Keep only the
  serial path (which already exists for `workers <= 1`). The
  generator now always yields `(chrom, sites)` serially regardless
  of `--workers`. Note the 5f auto-pick / auto-derate stays — chunk
  size and chunked simulation still apply per chromosome.
- [ ] **Tests** in `tests/test_cohort_parallel_write.py`:
  - **Determinism:** same `--seed` produces byte-identical
    `cohort.chr*.bcf` regardless of `--workers ∈ {1, 2, 4, 8}`.
    Diff via `bcftools view -H` md5 across two runs.
  - **End-to-end:** small cohort (n=20, chr22, 1 Mb,
    `--demo-model none`) round-trips through the parallel-write
    path; per-sample columns match a serial-write reference.
  - **Memory bound:** regression test asserting peak RSS at
    moderate `n` doesn't scale with `--workers` (within
    constant-factor noise). Marked optional / `@unittest.skipIf` on
    CI runners that don't expose `psutil`.
  - **Sample-slice helpers:** `dense_gts_from_carriers_slice`
    matches `dense_gts_from_carriers` on the relevant slice for a
    range of `slice_lo`/`slice_hi`.
- [ ] **Docs:** README + TUTORIAL — `--workers` now controls
  within-chromosome BCF-write parallelism in the cohort phase,
  not across-chromosome sim. Chromosomes simulate one at a time
  in the parent process.

#### Phase B — workers walk the tree sequence (required to make Phase A work at n=3000+)

**Status (2026-05-09):** promoted from "deferred / nice-to-have" to
"required follow-up after PR #28" by the memprof26 trace. At
`--n 3000 --workers 8`, Phase A OOM-kills on chrom 1: parent's
~17 GB sites list is fork-shared to 8 workers, but CPython's
reference-counting writes `ob_refcnt` on every PyObject access —
*including read-only iteration*. Workers iterating the shared
sites list trigger COW page faults across the whole structure.
Physical RAM grows ~linearly with worker count and exhausts on a
32 GB host within ~5 s of fork. Phase B sidesteps this entirely
because workers don't share a giant Python object graph — each
walks `ts.variants(samples=slice)` over the fork-shared tskit
TreeSequence, whose data is held in numpy arrays that *are*
fork-COW-safe (one Python wrapper, not one PyObject per element).

**Architecture:**

1. Parent simulates tree sequence (existing 5f code, unchanged).
2. Parent does *not* materialise the sites list. Instead, it forks
   W workers immediately after sim.
3. Each worker iterates `ts.variants(samples=[slice_lo*2 ..
   slice_hi*2])` for its sample slice, applies overlays per-site
   (overlay tables are ~500 MB and fork-COW-share cleanly because
   they're read-mostly Python dicts touched once per site, not per
   element), and writes its slice's partial BCF.
4. Parent runs `bcftools merge -O b` exactly as in Phase A.

**Wall-time estimate:** at n=3000 × 70 Mb × W=8: ~150 s/chrom (each
worker walks the tree once at ~100 s + writes its slice at ~50 s,
all in parallel). Compare to Phase A's measured 940 s/chrom serial
write or its theoretical ~150 s parallel write. Same wall time as
Phase A *would* have delivered, except actually achievable.

**What we lose vs Phase A's ideal (the tradeoffs):**

1. **CPU duplication.** Each worker re-walks the tree to materialise
   its slice's variants. At ~100 s tree-walk × W workers = 800
   worker-seconds of duplicated CPU per chromosome. Wall-time it's
   parallel-fast, but on metered cloud cycles you're paying ~3×
   more CPU than Phase A's "parent walks once" design. Workstation
   users feel this as fan noise, not budget.
2. **Overlay application moves into workers.** Today (Phase A)
   parent applies ClinVar (282k records) and rsID (440k records)
   overlays once before fork; everyone reads the same overlaid
   sites. In Phase B, each worker applies overlays per-site as it
   walks. **Constraint:** every overlay function must be
   deterministic and side-effect-free, so identical inputs yield
   identical outputs in every worker. Today's overlay code already
   is, but this becomes a hard contract for any future overlay
   logic.
3. **AC computation across slices.** Allele count (AC) for a site
   is the sum of carriers across all slices. Today (Phase A) parent
   computes AC from the full carriers list before workers see it.
   In Phase B, each worker only sees its slice's contribution to
   AC. Either: (a) workers each write their slice's partial AC and
   `bcftools merge` re-derives the merged AC (works — already
   tested in PR #28's parity tests), or (b) workers compute their
   slice's AC and write it explicitly into the partial BCF for
   later sum. (a) is what we already do; no change needed.
4. **Sparse `carriers` list as a debugging artifact disappears.**
   Today you can pickle the sites list and inspect carriers
   post-hoc for diagnostics or downstream stats. After Phase B,
   carriers exist only inside per-worker iteration scratch and are
   gone after the BCF write. Inspection requires re-walking the
   tree. Acceptable but worth noting.
5. **Admixture path needs separate refactoring.**
   `admixture.simulate_cohort` emits per-person ancestry segments
   from the same tree walk. Under one parent-side walker those
   segments come out alongside variants. Under W walkers, each
   worker emits its slice's ancestry segments — needs an
   ancestry-merge step (or each worker writes its own per-slice
   ancestry file and parent concatenates). Treat as a separate
   follow-up PR; landing Phase B for the non-admixture path first
   is fine.
6. **Site ordering source-of-truth moves.** Today parent sorts the
   sites list once. In Phase B, each worker walks in tree-position
   order; bcftools merge re-sorts on `(chrom, pos, ref, alt)`. This
   is functionally equivalent but it's a behaviour change that
   could surface ordering-sensitive bugs in overlay code or
   downstream consumers. Add a regression test: full pipeline
   produces byte-identical BCF (after `bcftools query` field
   normalisation) under W ∈ {1, 2, 4, 8}.
7. **Test refactoring surface.** Several tests build a sites list
   manually and feed it to `CohortBcfWriter`. Phase B either keeps
   that path supported (parent-side sites list as an *optional*
   input for tests) or migrates the affected tests to TreeSequence
   fixtures. ~150 LOC of test work.
8. **Apparent-RSS still inflates in monitoring.** Phase B reduces
   *physical* RSS to "1 × tree sequence + W × per-slice scratch",
   but each worker still inherits the parent's address space at
   fork. `/proc/<pid>/status` will still show each worker with the
   parent's RSS. Linux kernel OOM-scoring also still treats COW
   pages as resident in each child. Phase B fixes the *physical*
   exhaustion (which is what kills the run) but not the apparent
   inflation (which trips monitoring tools and cgroup limits).
   Operators on cgroup-bounded hosts will still need to size memory
   limits as `~W × tree_sequence_size` for the cgroup, even though
   physical use is ~1 ×.

**Bottom line.** Phase B is the right architecture for n=3000–30,000.
Above ~30k, even the per-worker tskit walk starts allocating
non-trivially (mutation count scales with n × sites), and the path
to n=1M shifts to **Phase 5d** below — streaming the tree sequence
to disk via Arrow IPC and never materialising the full per-slice
variant set in any single process. Phase B is the bridge; 5d is the
ceiling.

#### admixture mirror (separate follow-up)

`admixture.simulate_cohort` uses `BinaryMutationModel` and emits
per-person ancestry segments alongside cohort sites. Local-ancestry
tracking interacts with the tree sequence walk in ways the
non-admixture path doesn't, so the refactor is more delicate.
Sample-slice writes work cleanly for the cohort BCF side; ancestry
segments would either stay in the parent (parent emits them after
the tree walk) or get split per-worker. Separate PR after Phase A
lands.

### What this means for `--workers`

After Phase A, `--workers` semantics in the cohort phase change
from "parallel chromosomes" to "parallel sample-slice writers
within one chromosome at a time". For a fixed `--n`, doubling
workers no longer doubles peak RAM — it just speeds up the BCF
write. The 5f auto-derate-workers heuristic remains useful only
for pathological `--workers` choices in the chunked-simulation
path (which still runs in parent); it's a no-op for the cohort
write path because RAM doesn't scale with `--workers` there.

The fanout phase's `--workers` semantics are unchanged: still one
worker per person VCF write, batched per `--fanout-batch-size`.

### Resolved decisions

(Migrated from "Open questions for review" — answers locked in
during the 2026-05-08 design review.)

1. **Per-worker output format: per-slice partial BCF + `bcftools
   merge`.** Each worker writes its slice through the existing
   `CohortBcfWriter` (a thin sample-slice variant), so the
   subprocess pipeline is parity-identical with the current writer.
   Parent then runs `bcftools merge -O b` to join partials by
   `(chrom, pos, ref, alt)` key — every partial has the same site
   set in the same order, so merge collapses to a sample-column
   join. Pickle-and-rejoin-in-Python was rejected because it
   doesn't avoid disk I/O at this scale (the BCF still has to be
   written) and forces a custom binary format.
2. **Sample-slice (not position-slice).** Equal work per worker,
   deterministic split, falls out cleanly from `sample_ids[lo:hi]`.
   Position-slice would force per-worker variable workload (allele
   frequency varies with chromosome region) and wouldn't help the
   admixture follow-up.
3. **Parallel-chromosome opt-in flag deferred.** A
   `--parallel-chromosomes` flag (the pre-5e behaviour as a
   user-opt-in for 64+ GB hosts) was discussed but ruled out for
   Phase A. If a user explicitly asks for it later, add it as a
   small follow-up; it's a feature, not a regression to fix.
4. **`bcftools merge` (not `bcftools concat`).** The plan's earlier
   text loosely said "`bcftools concat` or comparable join" —
   that was wrong: concat is for region-disjoint VCFs (different
   genomic ranges per file), merge is for sample-disjoint VCFs
   (different sample columns per file). Sample-slice is the latter.

### Relationship to Phase 5d

Orthogonal axis. 5d (pysam-based direct binary BCF write) addresses
the *writer* throughput at n=1M+. 5e addresses the *simulation*
RAM at n=3000+. Both real follow-ups; 5e is more urgent because
it unblocks workstation-class users *now*, whereas 5d only matters
for cluster-class n.

---

## Phase 5g — batched per-person fanout

**Status (2026-05-08):** Phase 5g.1 + 5g.2 shipped (PRs #24, #25);
Phase 5g.3 (disk-spilled batch handoff) is the open follow-up.

**The bottleneck this addresses:** after Phase 5b made the cohort
BCF the canonical disk-backed handoff to per-person derivation,
the per-person fan-out itself became the dominant phase at
n=3000+. The original `derive_person_records` spawned
`bcftools view -s SID | bcftools view -e 'GT="ref"'` *once per
(person, chrom) pair* — at n=3000 × 22 chroms that was 66,000
bcftools subprocesses, each scanning a few hundred MB of
multi-sample BCF to keep one of 3,000 sample columns. Measured
wall: ~45 s/person with `--workers 1`, projecting to ~38 hours
fan-out for n=3000.

The plan didn't anticipate this — Phase 5b modelled per-person
derivation as a free disk read. The reality is that at high `n`
the multi-sample-BCF decode is the dominant cost on the
extraction side and the per-person VCF formatting is the
dominant cost on the write side.

### Phase 5g.1 — batched extraction (shipped, PR #24)

`derive_persons_batch(cohort_bcf_paths, sample_ids)` runs **one**
`bcftools query -s s1,...,sB -f '...[\t%GT]\n'` per chromosome
that emits all batch members' GT columns in a single decode pass;
the parser dispatches each row's GTs into per-person record lists
in the parent. The `_run_cohort_streamed` fan-out groups sample
IDs into batches and pipelines:

1. Parent calls `derive_persons_batch` — one bcftools subprocess
   per chrom for the batch's B sample IDs.
2. Parent stages `{sid: records}` in `_PERSON_WORKER_STATE`.
3. Parent forks a fresh `ProcessPoolExecutor`. Workers fork-inherit
   the records via copy-on-write.
4. Workers consume their `sid`'s records and write per-person
   VCFs.
5. Pool exits, parent drops the batch's references, next batch
   starts.

Bcftools subprocess count: 66,000 → `(n / B) × n_chroms` (66k →
1,320 at B=50, 16,500 at B=4).

### Phase 5g.2 — safer default + memprofile marks (shipped, PR #25)

The PR-24 default of B=50 OOM'd a worker on a 32 GB host with
`--workers 8`: parent extracts 14 GB of records (per-person ~280
MB at n=3000), 8 workers fork from a 14-GB-RSS parent, kernel's
per-process OOM scoring counts COW pages as resident in each
worker, total apparent RSS hits 8 × 14 = 114 GB → reap.

The binding ceiling is `(parent_baseline + B × per_person) ×
workers`, not `B × per_person`. PR #25 dropped the default to
B=4 and added per-batch memprofile marks
(`batch N stage A start / extracted / pool spawned / done`) so
the next trace pinpoints per-stage peak with no guesswork.

### Measured pre-#28 baseline (memprof25, 2026-05-09)

User trace at `--n 3000 --chromosomes 1-22 --chr-length-mb 70
--workers 1 --fanout-batch-size 7`, run on a 32 GB host, completed
successfully. Captured before PR #28 (Phase 5e Phase A) merged, so
cohort BCF write was still serial; serves as the canonical
reference point for measuring future improvements.

| Phase | Wall | RSS plateau |
|---|---|---|
| Setup (clinvar / candidates / rsid pool) | 0 → 16 s | 567 MB |
| Cohort (sim + serial BCF write × 22 chroms) | 16 → 27,038 s ≈ **7.5 h** | ~17 GB per chrom |
| Cohort → fanout boundary | (instant) | 17 GB → **11.9 GB** (~5 GB freed) |
| Per-person fanout (W=1, B=7) | 27,041 → 92,665 s ≈ **18.2 h** | 11.89 GB rock-steady |
| **Total** | **25.7 h** | — |

Per-chrom split inside the cohort phase: msprime sim + tree-walk
≈ 290 s, serial BCF write ≈ 940 s. The BCF write is **76 % of
per-chrom wall** — exactly what Phase 5e Phase A (PR #28)
parallelises. Per-batch fanout: ~150 s for B=7 = **~21 s/person**
single-threaded (vs ~45 s/person pre-5g.1, confirming 5g.1+5g.2
held).

**Course-of-action implications captured at the time:**

- **Fanout is now the dominant cost** — 18.2 h / 25.7 h ≈ 71 % of
  total runtime. Validates 5g.3 as the next high-impact unit of
  work after PR #28 lands. Linear extrapolation at W=8 (if RAM
  permitted): 18.2 h → ~2.3 h; with B=50 W=8 after 5g.3 the
  earlier ~1 h estimate looks consistent.
- **Per-chrom cohort plateau at W=1 is ~17 GB** — matches the
  Phase A RAM-bound estimate above (line ~1112). This exceeds the
  5f' constant-term hypothesis (~3.6 GB OOA n=3000) because that
  was framed around per-chunk constants, not full-chromosome
  state. Worth flagging as additional motivation for Phase B (see
  next subsection) since per-process OOM scoring on Linux counts
  COW pages as resident.
- **The cohort → fanout 5 GB drop** is sites-list /
  tree-sequence / overlay state being released. A ~12 GB parent
  baseline at fanout-fork time is what 5g.3 inherits — see the
  5g.3 RAM-bound discussion below.

### Phase 5g.3 — disk-spilled batch handoff (planned)

**Goal:** decouple parent's RSS at fork time from the batch size.
After 5g.3, B and W can both be large simultaneously without
multiplying apparent RSS.

**Branch (when implementation starts):**
`perf/phase5g-disk-spill-fanout`

**Architecture:**

1. Stage A (parent): `derive_persons_batch_to_disk` writes per-
   person records to per-person tempfiles in
   `out/.fanout-staging/person_<i>.tsv` (one row per record,
   tab-separated CHROM/POS/ID/REF/ALT/INFO/GT). Parent RSS stays
   at baseline because records are streamed straight to disk, not
   buffered in a Python dict.
2. Stage B: parent forks workers from a thin parent (~600 MB
   baseline). Each worker reads its assigned tempfile, parses
   records into the dict shape `write_person_vcf` expects, runs
   the existing per-person VCF write.
3. Worker deletes its tempfile on success (truth-bed and per-
   person VCF are the canonical outputs; staging is intermediate).

**RAM bound after 5g.3:**

```
peak parent RSS ≈ ~600 MB (baseline only, no batch records held)
peak per-worker RSS ≈ ~600 MB (parent shared) + ~280 MB private records
peak system RSS  ≈ baseline + W × ~280 MB
```

For `n=3000` × `W=8`: ~2.8 GB instead of the current 14 GB at B=4.
Headroom for `B=50+` × `W=8` simultaneously.

**Pre-fanout state release is part of the 600 MB target.**
memprof25 shows the parent's RSS at fanout-fork is ~12 GB
(11.89 GB), only ~5 GB shed at the cohort → fanout boundary out
of the ~17 GB cohort-phase peak. The remaining ~11 GB is a mix of
sites-list references, ClinVar / rsID overlay tables, and Python
allocator fragmentation that survives the post-cohort GC. To hit
the 600 MB parent baseline assumed by the RAM bound above, 5g.3
also needs to drop those references explicitly before
`derive_persons_batch_to_disk` runs (`del sites; del overlays;
gc.collect()` plus consider `malloc_trim` to release allocator
arenas). If we don't shed those, parent baseline at fork is ~12 GB
and the apparent-RSS-multiplier problem is unchanged from today —
the disk-spill only helps with the *batch records* term, not the
parent-baseline term. Treat the cohort-state-release path as a
required sub-task of 5g.3, not an optimisation.

**Disk cost:** ~280 MB × n_records-per-person × n_persons. At
n=3000 that's ~840 MB total staging at any one time (one batch
spilled at a time); at n=100k it's ~28 GB peak staging. Cleanup
is per-person on success; staging dir gets `rm -rf`'d at the end
of the fan-out. Acceptable on cloud instances; surface a clear
disk-space requirement note in `--help`.

**Wall-time impact:** the per-person VCF write itself (~70 s/
person, dominated by `draw_site_quality` Python RNG calls) is
unchanged. 5g.3 just removes the `B × W` ceiling so workers stay
fully busy. With `B=50, W=8`: fanout for n=3000 should drop to
~1 hour (vs ~8 h at B=7,W=8 currently, vs ~16 h at B=4,W=8).

**Out of scope for 5g.3 — the deeper win waiting after this:**
numpy-vectorise `draw_site_quality` (replace the per-record
Knuth-poisson Python loop with `np.random.poisson(lam, n)` etc.).
~5–7× speedup on the per-person VCF write; brings fanout from
~1 h to ~10–15 minutes at n=3000 × W=8. **Trade-off: breaks bit-
exact reproducibility** because numpy's RNG stream differs from
Python's `random.Random`. Same statistical properties, different
output bytes. Worth it for wall-clock-sensitive workflows; not if
downstream tooling pins golden-hash equality. Tracked here as a
candidate future PR rather than a confirmed task — needs a
product call from the user first.

#### Phase 5g.3 tasks

- [ ] **`cohort_derivation.py`:** add
  `derive_persons_batch_to_disk(cohort_bcf_paths, sample_ids,
  staging_dir)` mirroring `derive_persons_batch` but writing each
  parsed record line to `staging_dir/person_<sid>.tsv` instead of
  appending to an in-memory list. Return the dict
  `{sid: tsv_path}` for workers.
- [ ] **`cli.py`:** swap the in-memory `batch_backgrounds` dict
  for `batch_staging_paths`. Worker reads its tempfile + parses
  records back to dicts. Tempfile gets deleted on worker success.
  Staging-dir cleanup at fan-out end.
- [ ] **`--fanout-batch-size` default revisit:** with the
  parent-RSS ceiling gone, default can rise. 50 was the PR-24
  pre-OOM target; pick a default that balances bcftools
  invocation count (smaller batches = more) vs disk staging
  footprint (smaller batches = less peak staging). Tentative
  default: 50 once the disk-spill is in place.
- [ ] **Tests** in `tests/test_cohort_derivation.py`:
  - `derive_persons_batch_to_disk` parity with the in-memory
    `derive_persons_batch` — same per-person records, just on
    disk.
  - Tempfile cleanup on worker success and on worker failure
    (failure mode shouldn't leave stale staging).
  - Disk-space check that the staging dir size scales as
    `O(B × per_person)` not `O(n × per_person)`.
- [ ] **Docs:** README + TUTORIAL notes on staging-dir disk
  requirement; troubleshooting section for "fanout phase running
  out of disk".

### Phase 5g.4 — single-pass cohort-to-person dispatch (planned)

**Goal:** break the per-batch full-cohort-BCF scan that makes the
current fan-out quadratic in `n`. After 5g.4, the cohort BCFs are
scanned exactly **once total** during the per-person phase
(regardless of `n`), and per-person VCFs are emitted from local
staging streams.

**Branch (when implementation starts):**
`perf/phase5g.4-dispatch-prototype`

#### Why this supersedes 5g.3

Phase 5g.3 spills the *batch* dict to disk so that B and W can
grow on small-RAM hosts. That's useful at n ≤ 10k where the
per-batch scan cost is small relative to wall time. **It does not
help at n ≥ 50k** because the dominant cost is the cohort BCF
scan itself, which is repeated `n / B` times. At n=100k the
per-batch wall time is ~60–80 min regardless of B (it's I/O- and
decode-bound on the full cohort); making B bigger reduces batch
*count*, not batch *cost*. 5g.4 eliminates the batch-and-rescan
loop entirely.

If engineering bandwidth is constrained, **skip 5g.3 in favor of
5g.4.** 5g.3's disk-spill becomes architecturally redundant once
5g.4 lands (the dispatch path never holds a batch dict in
parent RAM in the first place). 5g.3 was sized for n=3k–10k
users; 5g.4 is the n=100k+ path.

#### The user observation that triggered this

**Reported run (2026-05-18):** `--n 100000 --chr-length-mb 70
--chromosomes 1-22,X --cohort-mode arrow-streaming` on a 32 GB
host. Cohort phase completed in ~4 days. Per-person fan-out
projected to **~4 years** at the default `--fanout-batch-size 4`:

```
person VCFs: 1/100,000 written
  (0.0/s, elapsed 1h 3m 58s, eta 106601h 12m 43s)
```

User retry with `--fanout-batch-size 100 --workers 12` was OOM-
killed during stage B pool spawn. Memprofile and dmesg
confirmed the failure mode:

| Event | Value |
|---|---|
| Parent RSS after stage A extract (B=100) | 29.4 GB |
| Apparent total RSS at pool spawn (12 fork-children) | 382 GB |
| Per-OOM-killed-worker `anon-rss` | 30.1 GB |
| dmesg `constraint` | `CONSTRAINT_NONE, global_oom` |

The fork-COW pages did not stay shared. CPython's per-object
refcount header is in the same page as the object data, so every
read access bumps the refcount and writes the page, which COW
promotes to private within seconds of worker startup. Realistic
fork-pool budget on a small-RAM host is therefore:

```
RAM_GB ≳ (1 + n_workers) × per_person_GB × B + ~4 GB system
```

**At the per-person bundle size measured in this run (~0.29 GB
at n=100k × 70 Mb × 23 chroms):**

| Workers | Safe B on 32 GB host | Batches | ETA |
|---|---|---|---|
| 1 | ~40 | 2,500 | ~111 days |
| 4 | ~15 | 6,667 | ~296 days |
| 12 | ~5 | 20,000 | ~890 days |

Maximum B is what controls wall time (per-batch scan time is
~constant in B), so **`workers=1, B=as-high-as-RAM-allows`** is
the optimal tuning *for the current algorithm.* On a 32 GB host
that yields ~4 months — the algorithmic floor today.

#### Cheap parallel quick-win — `gc.freeze()` (Phase 5g.4.0)

Before the dispatch refactor lands, a 3-line patch to
`cli.py:_run_cohort_streamed` can recover a meaningful chunk of
the lost COW sharing. Pattern:

```python
import gc
# parent: just before pool spawn
gc.collect()
gc.freeze()
with ProcessPoolExecutor(...) as ex:
    ...
```

`gc.freeze()` (Python 3.7+) marks all currently-tracked objects
as permanent so the garbage collector stops traversing them.
Most of the refcount writes that defeat fork-COW come from GC
generational sweeps, not from data access. Realistic gain on
this workload: **~30–50 % more effective B for the same RAM**,
i.e. roughly **2× speedup on small-RAM hosts** by enabling
larger batches.

This is a low-risk PR that can ship in days and is
architecturally independent of 5g.4. **Ship it as
`perf/gc-freeze-fanout` separately.** Caveats: workers still
write refcounts on data access (e.g. `dict[key]` lookup bumps
the looked-up value's refcount), so the win is partial; and any
mutation of frozen objects post-freeze defeats it, so the
`batch_backgrounds` build must complete before the freeze call.

#### Architecture (Phase 5g.4 main work)

1. **Dispatch pass (single scan).**
   `dispatch_cohort_to_staging(cohort_bcf_paths, sample_ids,
   staging_dir, chunk_chrom_at_a_time=True)` runs:
   ```
   bcftools query -s s1,s2,...,sN \
     -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO[\t%GT]\n' \
     cohort.chr<X>.bcf
   ```
   once **per chromosome BCF, with all N sample IDs in
   `-s`**, not split into batches. Each row is parsed and the
   per-sample GTs are dispatched into per-person append-mode
   writers. Records the same per-person record shape that the
   current `derive_persons_batch` parser emits.

2. **Staging format.** One file per (person, chrom) for the
   chunked variant; or one file per person for the non-chunked
   variant. Compact packed binary preferred over TSV — per-
   record cost ~12 bytes (4-byte pos + 1-byte GT + small
   per-record payload) instead of the ~30 bytes a TSV line costs.
   At n=100k × ~830k records/person that's ~10 GB total staging
   (chunked) vs ~2.5 TB without chunking.

3. **Per-chrom dispatch chunking.** Default mode: dispatch chr1,
   emit chr1 fragment for every person, delete chr1 staging, then
   chr2, etc. Each per-person final VCF is assembled from per-
   chrom bgzip fragments via `bcftools concat` at the end. Peak
   staging disk: ~per-chrom slice ≈ 400–800 GB at n=100k WGS.
   With aggressive compression (zstd or bgzip on staging
   fragments) that drops to ~200–400 GB.

4. **FD-limit handling.** Opening 100k staging file handles
   simultaneously exceeds typical `ulimit -n` of 1024–4096.
   Dispatch in rotating windows of ~1000 persons at a time per
   chrom-scan-pass; on a window flush, close the window's
   writers and open the next 1000. Each chrom is rescanned per
   window, so for a 1000-wide window at n=100k that's 100
   re-scans per chrom (~25 chroms × 100 = 2500 scans total).
   Even this is **vastly better** than today's `n / B` ≈ 25,000
   re-scans at B=4. If the user has bumped `ulimit -n` to 1M,
   collapse to a single window and a single scan per chrom.

5. **Per-person VCF emit (stage B).** Workers each read their
   own staging file(s), parse into the record dict shape that
   `write_person_vcf` already consumes, and emit the per-person
   VCF with DP/GQ/AD noise model layered identically to today.
   Staging files are deleted on worker success.

#### Expected wall time at n=100k

| Phase | Approximate wall time | Bound |
|---|---|---|
| Cohort scan (single pass, all chroms) | ~5–15 min | NVMe sequential I/O |
| Dispatch + staging write | ~1–3 hours | Sequential write of ~10 GB (chunked + compressed) |
| Per-person VCF emit (parallel) | ~12–48 hours @ 8–16 workers | CPU on `write_person_vcf` |
| **Total** | **~15–60 hours** | (vs ~3 years today) |

**Order-of-magnitude: ~3 years → ~1–3 days. Roughly 500–1000×
speedup.** The new bottleneck after 5g.4 is per-person VCF
emit CPU (DP/GQ/AD draws + bgzip), which is exactly the right
shape — local to each worker, parallelisable, doesn't scale
with `n` per person.

#### Disk cost

| Mode | Peak staging disk | Notes |
|---|---|---|
| Per-chrom chunked + compressed | ~200–400 GB at n=100k WGS | Recommended default |
| Per-chrom chunked, uncompressed | ~400–800 GB | Trade CPU for disk |
| Whole-genome staging | ~2.5 TB at n=100k WGS | Avoid; only viable with very large scratch |
| Whole-genome compressed | ~1 TB | Acceptable if disk is cheaper than per-chrom concat CPU |

At n=1M the chunked-compressed numbers scale linearly: ~2–4 TB
peak. Surface this clearly in `--help` and the docs.

#### Data quality

**Zero loss if implemented correctly. Output should be byte-
identical to the current path.** The simulator's randomness
(msprime, person_seeds, DP/GQ/AD noise, MT lineage carrier
draws) is entirely determined by inputs that don't change.
5g.4 only changes how records are *routed* between cohort BCFs
and per-person VCFs.

Implementation hazards that **must** be preserved or quality
degrades:

1. **Record ordering** — per-person records must reach
   `write_person_vcf` in the same chrom/pos order as today.
   Genome-sorted output from `bcftools query` per-chrom ensures
   this trivially.
2. **M13.5 MT clonal-inheritance contract** — every MT record
   must reach every person regardless of original simulator GT
   (`cohort_derivation.py:230-239`). The dispatch path must
   carve MT out of the hom-ref skip filter exactly the same way.
3. **`afs=[None]` semantics** — `write_person_vcf`'s MT lineage
   carrier fallback (AF=0.1) depends on the `afs=[None]` shape
   that today's `derive_persons_batch` emits. Staging must
   preserve this field.
4. **INFO field carriage** — staging must carry the full
   `_INFO_FIELDS_TO_CARRY` set (`CLNSIG`, `CLNDN`, `COSMIC_ID`,
   `COSMIC_GENE`, `SVTYPE`, `SVLEN`, `END`, `CIPOS`). Dropping
   any of these silently degrades ClinVar / COSMIC / SV
   overlay fidelity.
5. **Resume determinism** — partial-then-restart must reproduce
   identical record ordering. Keyed on (chrom, pos, alt_idx)
   yields trivially identical content; checkpoint after each
   chrom's dispatch completes.

**Validation harness** is the gating criterion for ship:
n=100 and n=1000 runs with `--cohort-mode arrow-streaming` and
`--cohort-mode dispatch` produce **byte-identical** per-person
VCFs (or, less strictly, identical (chrom, pos, ref, alt, GT,
INFO) tuples per person). PR-B does not merge until the diff
is empty across both n sizes.

#### Phase 5g.4 tasks

##### PR-A: `perf/phase5g.4-dispatch-prototype`

- [ ] **`syntheticgen/cohort_dispatch.py` (new module):**
  - `dispatch_cohort_to_staging(cohort_bcf_paths, sample_ids,
    staging_dir, *, chunked=True, window_size=None)` — single-
    pass scan and dispatch to per-person staging.
  - `read_person_staging(staging_path) -> list[dict]` —
    consumer that returns the same record-dict shape today's
    `derive_persons_batch` returns.
  - Compact packed staging format with `_INFO_FIELDS_TO_CARRY`
    preserved.
- [ ] **`syntheticgen/cli.py`:** new `--cohort-mode dispatch`
  choice. When selected, replace the `_run_cohort_streamed`
  per-batch loop with the dispatch+consume pipeline.
- [ ] **`syntheticgen/resume.py`:** bump `_SCHEMA_VERSION` to 4.
  Add `dispatch_state: dict[str, str]` (per-chrom dispatch
  status: `pending` / `staged` / `consumed`). Preserve all v3
  fields. Migration from v3 = "treat unknown dispatch_state as
  empty, restart dispatch from scratch."
- [ ] **`tests/test_cohort_dispatch.py`:**
  - byte-equivalence against `derive_persons_batch` at n=10
    and n=100 (the load-bearing test for "zero quality loss")
  - per-chrom resume mid-dispatch
  - disk-full mid-dispatch surfaces a clear error, no partial
    staging silently left behind
  - FD-limit handling: forced low ulimit (e.g. `resource.setrlimit`)
    completes via window rotation
  - MT carve-out preserved (same-lineage record-set contract
    pinned by existing M13.5 tests)
  - `afs=[None]` carry-through preserved (MT lineage fallback
    still fires)
  - `_INFO_FIELDS_TO_CARRY` full set carried (ClinVar / COSMIC /
    SV overlay metadata round-trips)

##### PR-B: `perf/phase5g.4-validation-and-promote` (stacked on PR-A)

- [ ] **Larger-scale validation:** n=1000 byte-diff run vs
  `arrow-streaming`; record the result in
  `DATA_QUALITY_ASSESSMENT.md`.
- [ ] **Benchmark numbers in this doc:** measured wall time at
  n=10k and n=100k on the project's reference hardware.
- [ ] **`--cohort-mode auto` promotion:** at `n ≥ 50000` (or
  configurable threshold), auto picks `dispatch` over
  `arrow-streaming`. Lower thresholds stay on `arrow-streaming`
  until the dispatch path has more bake time.
- [ ] **Docs:** README + TUTORIAL updates for the new mode and
  staging-disk requirement; troubleshooting section for "dispatch
  phase out of disk" and "FD limit too low."
- [ ] **`PERFORMANCE_PLAN.md` revision:** flip 5g.4 status to
  "shipped" with measured numbers; mark 5g.3 as superseded.

##### PR-C: `perf/gc-freeze-fanout` (independent, can ship first)

- [ ] **`syntheticgen/cli.py`:** add `gc.collect(); gc.freeze()`
  before the `ProcessPoolExecutor` spawn in
  `_run_cohort_streamed`'s stage B.
- [ ] **`tests/test_cli_modes.py`:** smoke test that fan-out
  still completes at n=10 with the freeze in place. Behavioral
  parity, no determinism break.
- [ ] **Note in this doc:** measured RAM savings on a small
  reference run; updated effective-B table on a 32 GB host.

#### Out of scope for 5g.4

- **Per-person VCF emit speedup** (the new bottleneck after
  5g.4). The `draw_site_quality` Knuth-poisson Python loop is
  the next target — `numpy.random.poisson(lam, n)` is ~5–7×
  faster, but breaks bit-exact reproducibility against today's
  golden hashes (numpy's RNG stream differs from `random.Random`).
  Tracked as a candidate future PR; needs a product call on
  whether bit-exact replay matters more than wall-clock.
- **In-memory Arrow transpose alternative** (the "Option 2"
  we considered). Avoids disk staging but reintroduces the
  memory-pressure failure mode the dispatch path was designed
  to escape. **Explicitly ruled out** for the n=100k+ workload
  this phase targets, per user direction (2026-05-19).
- **Single-pass dispatch directly from msprime tree sequence**
  (skipping cohort BCFs entirely). Would cut another ~4 days
  off the n=100k path but requires re-plumbing the cohort →
  overlay → per-person pipeline. Track as a future Phase 6+
  candidate after 5g.4 has settled.

---

## Phase 5d — streaming-mmap cohort intermediate (path to n=1M)

**Status (2026-05-15):** **Phase 5d.1 shipped.** `--cohort-mode
arrow-streaming` is in production; the cli + cohort writer code
paths under `syntheticgen/cohort_arrow*.py` + `bcf_writer.py`
implement the streamed mmap'd Arrow IPC handoff described in this
section. Phase 5d.2 (pysam direct BCF writer) remains optional —
ship only after measuring whether the BCF-write phase is still
>25 % of cohort wall at n=100k+ in production.

**Original status note (2026-05-09):** viability spikes complete;
both PASSED. Phase 5d.1 implementation green-lit.

Scope expanded substantially after memprof26 demonstrated that
Phase 5e Phase A's "parent holds the sites list, workers fork-share
it" design fails at n=3000 due to CPython refcount-COW divergence.
Phase B (workers walk the tskit TreeSequence directly) bridges to
~n=30,000. Beyond that, even the per-worker tskit walk allocates
non-trivially: mutation count grows with `n × sites`, and per-
iteration scratch state per worker starts hitting many-GB
territory at n=100k+. The design that actually scales to n=1M is
**streaming the cohort representation to disk and mmap-ing it from
workers** — parent never materialises a full in-memory site set,
workers never copy more than their slice.

**Goal:** make parent's RAM in the cohort phase O(1) per variant,
independent of `n`. Workers consume their slice's columns via mmap
without copying. Parent and worker memory both stay bounded as
`n` grows; only disk and I/O scale.

### Pre-implementation viability spikes (gate cleared)

Before committing to the ~600–800 LOC build, two standalone scripts
under [`scripts/spikes/`](scripts/spikes/) validated the load-
bearing assumptions in isolation. Both passed cleanly on the user's
32 GB workstation on 2026-05-09.

| Spike | What it tested | Result |
|---|---|---|
| 1 — [`spike1_mmap_fork_smoke.py`](scripts/spikes/spike1_mmap_fork_smoke.py) | OS-level: does `np.memmap` + `multiprocessing.fork` share physical pages across W workers via the page cache? | **PASS** — see [`spike1_results_2026-05-09.txt`](scripts/spikes/spike1_results_2026-05-09.txt). Total system RAM delta during the 8-worker mmap-read phase was 0 GB; apparent per-worker RSS sum (4.21 GB) tracked file size, not n_workers × file_size. |
| 2 — [`spike2_arrow_streaming.py`](scripts/spikes/spike2_arrow_streaming.py) | Arrow-specific: does `pyarrow.ipc.new_file` truly stream batches without buffering? Does `pa.memory_map` + `FixedSizeListArray.values.to_numpy(zero_copy_only=True)` give workers true zero-copy reads under a realistic per-site Python loop? | **PASS** — see [`spike2_results_2026-05-09.txt`](scripts/spikes/spike2_results_2026-05-09.txt). All four checks cleared: parent peak RSS during streaming write 184 MB (< 500 MB threshold), system RAM delta during workers +0.01 GB against a 7.83 GB apparent-RSS sum (apparent-to-physical ratio ~700:1), write throughput 245 MB/s (≥ 200 MB/s target), aggregate read throughput 1,355 MB/s (≥ 500 MB/s target). |

The Spike 2 result is the decisive evidence. The realistic per-site
Python loop — the same shape that killed Phase 5e Phase A through
refcount-COW divergence at n=3000 — ran cleanly on all 8 workers in
0.69 s with zero physical RAM growth. **Numpy holds one PyObject
wrapper for the whole mmap'd array, not one per element**, so per-
element index access doesn't trigger refcount-COW divergence. This
is the core architectural difference that makes Phase 5d safe where
Phase 5e Phase A wasn't.

The spikes did not exhaustively test every n=1M-relevant dimension
(only the default ~1 GB Arrow file size; no admixture path; no real
msprime input). Treat the green-light as "the mechanism works on
this host class" rather than "every implementation detail is
solved." Implementation-time discoveries are still possible; budget
for them.

A follow-up sub-spike ([`spike2b_pool_batch_matrix.py`](scripts/spikes/spike2b_pool_batch_matrix.py)
/ [results](scripts/spikes/spike2b_results_2026-05-10.txt))
characterised parent peak RSS during the streaming write as a
function of `batch_size`, memory pool, and `release_unused()`-
between-writes. Findings: `batch_size` is the only knob that
materially affects parent peak; pool choice (jemalloc vs system) and
`release_unused()` calls are noise. **Decision for 5d.1:**
`batch_size = 256`, default memory pool, no `release_unused()`
calls. The 5d task list and architecture pseudocode below reflect
this default.

### Why Apache Arrow IPC

Several mmap-able binary formats could work (custom struct format,
HDF5, Parquet, Arrow). **Apache Arrow's IPC ("Feather v2") format
is the right pick** for this problem:

1. **Columnar layout matches our access pattern.** Workers want to
   read one slice (range of sample columns) per site. Arrow stores
   each column contiguously, so a worker's read is a contiguous
   mmap range — no scatter, no per-site decode. Row-oriented
   formats (HDF5, custom struct) force every worker to scan every
   site's full row to find its slice columns.
2. **Zero-copy mmap is a first-class operation.**
   `pyarrow.ipc.open_file(path).read_all()` returns a Table whose
   columns are zero-copy views into the mmap'd buffer. Workers see
   numpy arrays without ever materialising Python objects. No
   refcount-COW issue (numpy arrays hold one PyObject wrapper for
   the whole array, not one per element).
3. **Variable-length carriers fit naturally.** Arrow's `ListArray`
   handles per-site variable-length child arrays via offsets — the
   exact shape our sparse `carriers` representation needs. No
   manual offsets bookkeeping.
4. **Streaming write API.** `pyarrow.ipc.new_file(...)` lets parent
   write one record batch at a time as it iterates
   `ts.variants()`. No requirement to hold the whole table in RAM
   to write it. Parent's RAM stays at one record batch's worth
   (~1k variants).
5. **Parquet is overkill, HDF5 is non-columnar, custom format is
   maintenance burden.** Arrow gives us the columnar mmap pattern
   without inventing or maintaining a binary format.

**New dependency: `pyarrow`.** Pulls in a ~70 MB wheel (Arrow C++
core + Python bindings). Heavier than the existing pure-Python
deps but lighter than pysam's compiled htslib build. Worth the
weight at the n=1M ask; conditional-import-and-degrade-gracefully
keeps it optional for n≤30k users who only need Phase B.

### Architecture

```
Parent (per chromosome):
  1. msprime.sim_ancestry(...) → TreeSequence ts
  2. ts.dump("scratch/chrom_N.trees")        # disk-backed tskit
     del ts                                   # free Python ref
  3. ts = tskit.load("scratch/chrom_N.trees") # mmap-backed reload
  4. with pyarrow.ipc.new_file(
         "scratch/chrom_N.cohort.arrow",
         schema=cohort_schema(n_samples)) as writer:
         for batch in stream_variants_to_arrow_batches(
                          ts, batch_size=256):
             writer.write_batch(batch)
     # Parent peak RSS during write ≈ ~9.5 bytes/element × batch_size
     # × n_samples (Spike 2b empirical fit). At batch_size=256:
     # n=10k → ~25 MB,  n=100k → ~250 MB,  n=1M → ~2.5 GB.
  5. Signal workers to start.

Workers (each, in parallel):
  1. table = pyarrow.ipc.open_file(
                "scratch/chrom_N.cohort.arrow").read_all()
     # Zero-copy mmap; table.columns are numpy views
  2. Slice columns to assigned [sample_lo, sample_hi]:
        slice_table = table.select(
            ["pos", "ref", "alt", "info_ac", ...,
             f"gt_{sample_lo}", ..., f"gt_{sample_hi-1}"])
  3. Stream-write partial BCF using existing CohortBcfWriter
     (or 5d.2's pysam direct writer if shipped) over the
     numpy-array view.
  4. Close mmap, exit.

Parent: bcftools merge partials → cohort/cohort.chr_N.bcf
        rm scratch/chrom_N.trees scratch/chrom_N.cohort.arrow
```

### RAM bound after 5d

```
peak parent RSS ≈ ~500 MB baseline + ~20 MB Arrow batch buffer
peak per-worker RSS ≈ ~500 MB baseline + ~50 MB BCF write buffer
peak system RSS  ≈ baseline + W × ~50 MB
```

**Independent of `n`.** At n=1M × W=8 the working set is ~5 GB
total — fits comfortably on any modern host. The cost moves to
disk.

### Storage and I/O — first-class constraints at this scale

**Per-chromosome scratch disk during processing:**
- Tree sequence dump: ~10-50 GB at n=1M × 70 Mb (tskit's
  table-store format; compresses well but still substantial).
- Arrow IPC intermediate: ~80 GB at n=1M × 70 Mb (uncompressed
  GT columns dominate; 1M samples × 500k sites × 1 byte = 500 GB
  raw, halved by Arrow's RLE / dictionary encoding for typical
  allele frequencies).
- Partial BCFs (W workers): ~10 GB per slice × 8 slices = 80 GB.
- Final merged BCF: ~80 GB.

**Peak per-chromosome scratch: ~250 GB.** Drop-after-merge pattern
keeps only one chromosome's intermediates in flight at a time; the
final BCFs accumulate at ~80 GB × 22 = ~1.7 TB.

**Total disk requirement at n=1M × 22 chroms:**
- Final outputs (cohort BCFs + per-person VCFs): ~2 TB
- Peak scratch (one chromosome in flight): ~250 GB
- Headroom buffer: ~10%
- **Recommended free disk: 2.5-3 TB**

**I/O bandwidth dominates wall time at n=1M.** Per-chromosome data
moves through disk twice (write Arrow, read mmap) plus tree
sequence dump and partial BCF writes. Total per-chrom I/O:
~400 GB. On NVMe (~3 GB/s) that's ~140 s — tolerable. On a
spinning disk (~150 MB/s) it's ~45 minutes per chromosome,
~16 hours just for I/O across 22 chromosomes. **Cloud-instance
disk type matters more than CPU at n=1M.**

**Surface this in `--help` and TUTORIAL:** users should know
upfront that running n≥100k requires a fast SSD/NVMe and 2-3 TB
free. Add a pre-flight check that warns if the output filesystem
is rotational or has insufficient free space.

### Phase 5d.1 vs 5d.2 split

**5d.1 — streaming-mmap intermediate (the load-bearing change):**
Replace parent's in-memory sites list with an Arrow IPC file.
Workers consume via mmap. This is what unblocks n=1M. Most of the
~600-800 LOC sits here:

- New `cohort_arrow.py`: schema definition + `stream_variants_to_
  arrow_batches(ts) -> Iterator[pyarrow.RecordBatch]` + reader
  helpers.
- `cli.py`: orchestration — dump tree, stream Arrow, fan out
  workers, merge.
- `bcf_writer.py`: a new mode that consumes Arrow column views
  rather than a Python list-of-dicts.
- Conditional `pyarrow` import; if missing, fall back to Phase B
  with a clear "n=1M requires pyarrow; install with `pip install
  pyarrow` or pass `--n ≤ 30000`" message.
- Pre-flight disk-space check + clear error.
- Tests: parity with Phase B output at small n; large-n smoke
  test gated on `RUN_LARGE_N=1`.

**5d.2 — direct binary BCF write via pysam (optional optimisation
on top of 5d.1):**
At n=1M the BCF write itself becomes a wall-time concern: piping
through `bcftools view -O b -` from a Python `subprocess.PIPE`
caps at ~50 MB/s (Python GIL + subprocess overhead). pysam's
`VariantFile.write()` writes BCF binary directly through htslib at
~500 MB/s, an order of magnitude faster.

This is **independent of 5d.1** in design: 5d.1's worker can use
either the existing subprocess pipeline or pysam. Ship 5d.1 first;
add 5d.2 only if the BCF write phase is measured to bottleneck
real n=1M runs (it might not — at NVMe disk speeds the partial
BCF writes already saturate disk before they saturate bcftools).

**Why pysam is cheaper to add now than at 5a evaluation time:**
back when 5a evaluated pysam, the codebase had no compiled deps.
5d.1 already ships pyarrow (a compiled wheel of similar size), so
adding pysam shifts the dep model from "pure Python" to "two
compiled wheels" — marginal cost. The 5a rejection no longer
applies cleanly.

### Out of scope for 5d entirely

- **Per-person fanout at n=1M.** The fanout phase has its own
  scaling story (5g.3 + the disk-spill design). At n=1M, fanout
  produces 1M per-person VCFs at ~50 MB each = ~50 TB. That's
  obviously a separate disk-space ceiling and a separate wall-time
  problem. 5d covers the cohort phase up to and including the
  cohort BCF; per-person derivation at n=1M needs its own design
  pass — likely chunked persons, each chunk sharing the cohort
  Arrow intermediate via mmap.
- **Distributing across machines.** Out of scope per the
  workload-assumptions section; if anyone needs to go beyond
  single-machine, the disk-backed Arrow + BCF intermediates make
  the data layer transferable, but the orchestration is a separate
  product.

### Phase 5d tasks

- [ ] **Branch:** `perf/phase5d-streaming-mmap-arrow` off main.
- [ ] **`cohort_arrow.py` (new):**
  - `cohort_schema(n_samples) -> pyarrow.Schema`: pos/ref/alt/qual/
    filter/info_ac/info_an/info_af + carriers as `ListArray<int32>`
    or per-sample GT columns as `int8[n_samples]`. Decision point:
    columns-per-sample is denser-on-disk but reads cleanly per
    slice; carriers-as-list is sparser-on-disk but workers must
    re-densify. **Tentative choice: per-sample int8 columns** —
    matches Arrow's mmap-slice strength.
  - `stream_variants_to_arrow_batches(ts, batch_size=256) ->
    Iterator[RecordBatch]`: walks `ts.variants()` once, emits
    one batch every 256 variants. Default chosen empirically from
    Spike 2b ([results](scripts/spikes/spike2b_results_2026-05-10.txt)):
    `batch_size` is the only knob that materially affects parent
    peak RSS, and 256 keeps the predicted peak at n=1M to ~2.5 GB
    vs ~9.4 GB at 1024 / ~37 GB at 4096, for a ~7% throughput cost.
    Memory bound (empirical fit): ~9.5 bytes/element × batch_size
    × n_samples per batch.
  - `read_slice(arrow_path, sample_lo, sample_hi) -> Table`:
    mmap-open + project columns. Workers' entry point.
- [ ] **`cli.py`:** new orchestration path gated on `--cohort-mode
  arrow` (or auto-pick for n≥100k):
  - Dump tree to scratch.
  - Stream Arrow IPC.
  - Fork workers; each writes partial BCF from its Arrow slice.
  - Merge + cleanup.
  - Pre-flight disk-space check; clear error if insufficient.
- [ ] **`bcf_writer.py`:** add `CohortBcfWriter.from_arrow_slice(
  arrow_table, slice_lo, slice_hi)` constructor that reads
  Arrow columns directly without a Python list-of-dicts
  intermediate.
- [ ] **Conditional pyarrow import:** if missing, raise
  `ImportError` with install hint; degrade to Phase B for
  n ≤ 30000.
- [ ] **Tests:**
  - Parity: output BCF byte-identical to Phase B's output at
    n=100, n=500, n=2000 across W ∈ {1, 2, 4, 8}.
  - Streaming bound: assert parent RSS during cohort phase stays
    under 1 GB independent of n (smoke test at n=10k).
  - Disk-space pre-flight: simulate insufficient space, assert
    clear error before any work begins.
  - Large-n smoke test gated on `RUN_LARGE_N=1` env: n=10k full
    pipeline completes; record wall time + peak disk.
- [ ] **Docs:** README + TUTORIAL — n=100k+ requires `pyarrow` +
  fast disk + 2-3 TB free; cloud-instance disk-type guidance.
  `--help` updates for `--cohort-mode` if surfaced.

### Phase 5d persistent regressions (track and verify post-merge)

5d eliminates most of the Phase B regressions for free — parent walks
the tree exactly once during the Arrow stream-write, so workers
never re-walk. That removes Phase B trade-offs **#1 (CPU
duplication)**, **#2 (overlay-in-workers determinism contract)**,
**#3 (AC-across-slices reconstruction)**, and **#8 (apparent-RSS
inflation per worker)**. The list below is what *does* persist into
5d and needs deliberate handling — captured here so we don't lose
them when 5d.1 implementation begins.

1. **Site-ordering source-of-truth moves.** Phase A: parent sorts
   the in-memory sites list once. 5d.1: parent writes Arrow in
   tree-position order; `bcftools merge` re-sorts on
   `(chrom, pos, ref, alt)`. Functionally equivalent, but a
   behaviour change vs Phase A that could surface ordering-sensitive
   bugs in overlay code or downstream consumers.
   - **Action:** add a regression test asserting full pipeline
     produces byte-identical BCF (after `bcftools query` field
     normalisation) under W ∈ {1, 2, 4, 8} on a small fixture.
2. **Sparse `carriers` list as a pickle-able debugging artifact
   disappears.** Phase A keeps a Python-pickleable in-memory
   carriers list that downstream diagnostics can introspect. 5d.1
   never materialises it — the on-disk Arrow file is the
   source-of-truth instead. Strictly *better* observability than
   Phase B (full GTs are recoverable per-site by re-mmap-ing the
   Arrow file), but anything that called `pickle.dump(sites_list)`
   for diagnostics needs to migrate to the Arrow file.
   - **Action:** audit any debugging / diagnostic code paths that
     consume the in-memory sites list. Provide a small
     `read_arrow_carriers(arrow_path, pos)` helper if the audit
     finds real consumers.
3. **Test refactoring surface (larger than Phase B's).** Several
   tests build a sites list manually and feed it to
   `CohortBcfWriter`. Phase B was estimated at ~150 LOC of test
   work; 5d.1 is bigger because we're introducing
   `CohortBcfWriter.from_arrow_slice(...)` as a new constructor.
   Either keep the manual-list path supported as an *optional* test
   input or migrate the affected tests to write a small Arrow
   fixture.
   - **Action:** during 5d.1 implementation, decide between
     dual-path support and full migration; document the call in the
     PR.
4. **Admixture path needs separate refactoring.**
   `admixture.simulate_cohort` emits per-person ancestry segments
   alongside cohort sites. Under the 5d.1 single-walker streaming
   design, ancestry segments come out alongside the variant stream
   in the parent — but the schema and downstream merge for ancestry
   is not yet designed. Land 5d.1 for the non-admixture path first,
   then mirror.
   - **Action:** track as a separate follow-up PR after 5d.1 lands.
     Same disposition as Phase B's admixture follow-up.

**Out-of-scope reminder (not a 5d regression — separate):** the
optional 5g.3 follow-up that numpy-vectorises `draw_site_quality`
*does* break bit-exact reproducibility (numpy RNG ≠ Python
`random.Random`). That regression is independent of which cohort-
phase architecture ships and lives or dies on its own product call.
Tracked in §5g.3.

### Phase 5d.2 tasks (optional, after 5d.1 measurements)

- [ ] **Branch:** `perf/phase5d2-pysam-bcf-write` off main (after
  5d.1 lands).
- [ ] **`bcf_writer.py`:** add a pysam-backed `write_site` path;
  benchmark vs subprocess pipeline.
- [ ] **Tests:** parity with the subprocess writer at modest n.
- [ ] **Decision gate:** ship only if measured 5d.1 BCF write
  phase is >25% of cohort wall at n=100k. Otherwise close as
  "not worth the dep".

---

## Phase 6 — SQLite for side-state at scale

**Goal:** consolidate the per-person side files
(`out/truth/person_NNNN.golden.bed`,
`out/truth/person_NNNN.noise.bed`, `out/ancestry/person_NNNN.bed`,
`out/manifest.json` with N entries) into a single cohort SQLite
database. At `n = 100 000`, the per-person-files layout becomes a
filesystem-stress problem (300k+ small files in a single directory)
that's slow to enumerate, slow to back up, and awkward to query.
SQLite handles a few hundred MB of structured per-sample rows
trivially and supports concurrent WAL-mode `SELECT`s for any
downstream "grade caller per-sample" workflow.

**Branch:** `perf/phase6-sqlite-sidestate`

### Why SQLite (not Parquet, not JSONL)

- WAL mode supports concurrent readers without a server.
- Schema is small and stable (samples / events / segments / manifest
  KV). Parquet is overkill and forces column-store thinking on
  small heterogeneous payloads.
- `sqlite3` is in the stdlib — zero new deps.
- Easy ad-hoc inspection during runs (`sqlite3 cohort.db "SELECT ..."`).
- Per-person BED derivation stays a one-liner via `bcftools`-friendly
  bedfile output: `dump_truth_bed(sample_id, db) → person.bed`.

### Tasks

- [ ] **Schema.** One DB at `out/cohort.db` with:
  ```sql
  CREATE TABLE samples (
    sample_id INTEGER PRIMARY KEY,
    sample_name TEXT NOT NULL UNIQUE,
    ancestry_summary_json TEXT
  );
  CREATE TABLE truth_events (
    sample_id INTEGER NOT NULL REFERENCES samples(sample_id),
    contig TEXT NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    kind TEXT NOT NULL,        -- HIGHLIGHTED / CLINVAR / COSMIC / SV / RSID / FLIP / DROPOUT
    payload_json TEXT NOT NULL
  );
  CREATE INDEX truth_events_sample_pos
    ON truth_events(sample_id, contig, start);
  CREATE TABLE ancestry_segments (
    sample_id INTEGER NOT NULL REFERENCES samples(sample_id),
    contig TEXT NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    hap1_pop TEXT,
    hap2_pop TEXT
  );
  CREATE INDEX ancestry_segments_sample_pos
    ON ancestry_segments(sample_id, contig, start);
  CREATE TABLE manifest_kv (
    key TEXT PRIMARY KEY,
    value_json TEXT NOT NULL
  );
  ```
- [ ] **Extend `truth.py` with a parallel DB writer.** Keep
  `TruthBedWriter` as the canonical path; add a `TruthDBWriter` that
  accepts the same events and inserts into `truth_events`. Run both
  in tandem when `--cohort-db` is set. Batches inserts in-memory and
  flushes per chromosome (matching Phase 5's chromosome-streaming
  model).
- [ ] **Extend admixture ancestry writer similarly.** BED files stay;
  the DB gets a parallel feed gated by `--cohort-db`.
- [ ] **Manifest.** Top-level fields stay in `out/manifest.json`. The
  DB's `manifest_kv` is a mirror for SQL-friendly access, populated
  at run end from the same data.
- [ ] **Variant-scan integration (optional).** Extend
  `nextflow_pipeline/bin/scan_variant.py` so it can read truth events
  for a sample directly from `cohort.db` instead of opening a per-
  person BED. Useful for caller-grading workflows at scale.
- [ ] **Tests.** Round-trip a small synthetic cohort through SQLite
  and assert `dump_truth_bed` output matches what the M11 BED writer
  produced directly. Concurrent-read test: 8 parallel processes
  query different sample IDs, assert no `database is locked` errors
  in WAL mode.
- [ ] **Docs.** TUTORIAL §9 (when each path activates),
  IMPLEMENTATION_PLAN architecture diagram, README Output layout.

### Open questions for review

- ~~DB primary or complement to BED files~~ resolved: **DB is a
  complement, not canonical** (locked in during plan review). BED
  files remain the default deliverable; the DB is an additional
  output that downstream tooling can opt into. Per-person BEDs at
  100k samples remain a filesystem-stress problem, but that's a
  scale tradeoff users opt into rather than a default behaviour
  change.
- Do we ship SQL helper views for common queries (e.g.
  `golden_per_sample_chrom`)? Probably yes, in a separate
  `synthetic_people/scripts/cohort_db.sql` reference file.

---

## Resuming after an interruption

1. `git branch --list 'perf/*'` — see which phase has an in-flight
   branch.
2. Open the matching section above and continue from the first
   unticked box.
3. If no `perf/*` branch exists, the next phase to start is the first
   one with any unticked boxes.
