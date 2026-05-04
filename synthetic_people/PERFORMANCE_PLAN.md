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

## Resuming after an interruption

1. `git branch --list 'perf/*'` — see which phase has an in-flight
   branch.
2. Open the matching section above and continue from the first
   unticked box.
3. If no `perf/*` branch exists, the next phase to start is the first
   one with any unticked boxes.
