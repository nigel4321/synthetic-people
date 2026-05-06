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

## Strategy summary

The phases compose, not compete:

| Phase | What | Scale unlock | Primary cost |
|---|---|---|---|
| 1 ✓ | Parallel chrom + per-person | Faster, not larger | None |
| 2 ✓ | Overlap I/O loaders with simulation | Faster, not larger | None |
| 3 (partial) | Dense numpy genotype matrix | `n ≈ 5 000` | New in-RAM shape |
| 4 | Hot-path cleanups | Faster | None |
| **5** | **Disk-backed cohort BCF** | **`n = 100k+`** | **+disk I/O at phase boundaries** |
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

### Tasks (when implementation starts)

- [ ] Add `--chr-chunk-mb` CLI flag (default `0` = unchunked,
  preserves today's behaviour for users with the RAM headroom).
- [ ] Refactor `coalescent.simulate_chromosome` (and its admixture
  twin) to optionally split into sub-chunks. Each sub-chunk runs
  the existing simulation pipeline against a contig slice.
- [ ] Per-chunk seeds derive deterministically from the chromosome's
  seed + chunk index, so `--seed` reproducibility holds.
- [ ] Each chunk's variants flow into the per-chrom BCF with
  positions adjusted by the chunk's start offset. Tree sequences
  are freed before the next chunk's simulation starts.
- [ ] **Tests.**
  - Chunked vs unchunked at small n produces equivalent per-record
    summary statistics (counts, Ti/Tv, AF distribution) within
    stochastic noise.
  - Memory bound: regression test runs `--n 1000 --chr-length-mb 70
    --chr-chunk-mb 10` and asserts peak RSS scales with one chunk's
    tree sequence (a few GB, not the unchunked 8-16 GB).
  - Chunk-boundary determinism: same seed + same chunk size →
    byte-identical BCFs across runs.
- [ ] **Docs.** README + TUTORIAL: explain `--chr-chunk-mb`
  semantics, document the cross-chunk LD caveat, recommend chunk
  sizes per host RAM budget.

### Open questions for review

1. **Default chunk size.** `0` (unchunked, today's behaviour) is
   the safest default — users opt in via `--chr-chunk-mb 10` when
   they hit OOM. An auto-pick based on n × demo_model + available
   RAM is friendlier but adds complexity; lean toward explicit-
   opt-in for the first cut.
2. **Chunk overlap.** Simulate adjacent chunks with overlapping
   margins to recover some cross-chunk LD? Adds 10-20% redundant
   simulation but recovers LD across boundaries up to overlap
   length. Worth it for downstream LD-aware analyses; out of scope
   for the first 5f PR.
3. **stdpopsim contig slicing API.** `species.get_contig(chrom,
   right=slice_end)` is documented, but `left=` may or may not be
   supported; check at implementation time. If `left` is missing,
   we either simulate `[0, chunk_end]` and discard `[0,
   chunk_start]` (wastes work) or use msprime's lower-level API
   directly.
4. **Sizing.** What's the host RAM the user is running on? At
   16 GB, even after 5f a chunk size ≤5 Mb is needed at n=3000 ×
   OOA_3G09; at 32 GB, 10 Mb chunks should fit comfortably.
   Knowing the target host helps pick a sensible default chunk
   size for the docs / TUTORIAL recipe.

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

---

## Phase 5e — within-chromosome parallel extraction over a shared tree sequence

**Status:** planned, **gated on Phase 5f landing first**. The user's
n=3000 × full-chr1 OOM happens *during simulation* of a single
chromosome — Phase 5e's parallel-extraction model only helps once
one tree sequence fits in RAM. 5f (chunked simulation) gets each
tree sequence to fit; 5e then parallelises extraction across that
smaller tree sequence.

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

### Tasks (when implementation starts)

- [ ] **Refactor `coalescent.simulate_cohort_iter`** along the
  structure above — simulate in parent, fork extractor pool, merge
  per-worker outputs into the per-chrom BCF.
- [ ] **`admixture.simulate_cohort` mirror.** The admixture path
  uses `BinaryMutationModel` and emits per-person ancestry segments
  alongside cohort sites; local-ancestry tracking interacts with
  the tree sequence walk, so the refactor is more delicate but
  follows the same principle.
- [ ] **Tests.**
  - Determinism: same seed → same per-chrom BCF byte-for-byte
    regardless of `--workers`. Sample slices split deterministically;
    merge order is fixed.
  - Memory bound: regression test that runs at moderate n and
    asserts peak RSS scales with one tree sequence (not workers ×
    tree sequence).
  - End-to-end: the user's failing command (`--n 3000 --chromosomes
    1-22 --chr-length-mb 70`) completes on a 32 GB host.
- [ ] **Re-measure** at the configurations 5c was measured against
  (`n=500 / 3000 / 10 000`) so the comparison is direct.
- [ ] **Docs** — README + TUTORIAL: clarify that `--workers` now
  controls within-chromosome parallelism (sample-slice extraction),
  not across-chromosome simulation; chromosomes simulate one at a
  time. Note that hosts with enough RAM for parallel chromosome
  simulation can opt into the old behaviour via a separate flag
  if we keep it.

### What this means for `--workers`

After 5e, `--workers` semantics change from "parallel chromosomes"
to "parallel sample-slice extractors within one chromosome at a
time". For a fixed `--n`, doubling workers no longer doubles peak
RAM — it just speeds up the (already-cheap) extraction phase. The
auto-derate-workers idea that surfaced as a stopgap during the 5c
discussion falls away because RAM stops scaling with `--workers`.

### Open questions for review

- **Per-worker output format.** Per-worker text-VCF chunks merged
  via `bcftools concat`, or per-worker pickle of records merged in
  the parent? Concat is simpler and matches the existing BCF
  writer; pickle keeps everything in Python and avoids the extra
  disk write. Lean toward concat for parity with 5b's pipeline.
- **Sample-slice vs position-slice.** Sample-slice gives equal
  work per worker assuming roughly uniform allele frequencies;
  position-slice gives genome-ordered output but unequal work.
  Recommend sample-slice because msprime's variant iteration is
  deterministic and equal slices fall out cleanly from the tree
  sequence's sample list. Sample-slice also fits admixture's
  per-person ancestry segments naturally (each worker emits its own
  slice's ancestry).
- **Keep parallel-chromosome path as opt-in?** Users with 64+ GB
  hosts and large chromosome counts get the old wall-time win for
  free today. Consider a `--parallel-chromosomes` flag that retains
  the pre-5e behaviour for those users; default off.

### Relationship to Phase 5d

Orthogonal axis. 5d (pysam-based direct binary BCF write) addresses
the *writer* throughput at n=1M+. 5e addresses the *simulation*
RAM at n=3000+. Both real follow-ups; 5e is more urgent because
it unblocks workstation-class users *now*, whereas 5d only matters
for cluster-class n.

---

## Phase 5d — direct binary BCF writes (path to n=1M+)

**Goal:** get past the writer-side bottleneck that emerges once
sparse storage takes RAM out of the picture. At `n=1M × ~500 k
sites` per chromosome, the text-VCF representation we currently
pipe through `bcftools view -O b` is roughly 2 TB per chromosome
(GT block is `n × ~3 B` per record). Bgzip throughput caps the
write at ~5 hours per chromosome, ~5 days total — slow even on
a fat workstation.

**Status:** stub only. Keep on the plan so we can revisit when the
1M scale becomes a real ask.

**Approach:** swap the `bcftools view -O b -` subprocess for a
direct binary BCF writer. Two viable libraries:

- `pysam` — the obvious choice. Adds a ~30 MB compiled wheel to
  the dep tree (its own htslib build with C extensions). Already
  evaluated and rejected in 5a as a transitive dep cost; revisit
  the cost/benefit at 1M scale.
- A hand-rolled BCF encoder using `htslib` via `ctypes`. Lower
  install cost, much higher implementation cost.

**Out of scope for the current PR.** Phase 5c gets to n=100k
cleanly; n=1M waits for 5d when it's actually needed.

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
