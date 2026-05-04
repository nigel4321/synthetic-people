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

---

## Phase 1 — concurrency, low risk

**Goal:** parallelise the work that's already independent. No
data-shape changes. Visible end-to-end runtime win on multi-core hosts.

**Branch:** `perf/phase1-concurrency`

- [ ] Parallelise per-chromosome msprime simulations
  - `coalescent.py:simulate_cohort` currently iterates chromosomes
    serially. Each chromosome is independent.
  - Use `concurrent.futures.ProcessPoolExecutor`,
    `max_workers=min(len(chromosomes), os.cpu_count())`.
  - Pre-derive one deterministic seed per chromosome from the master
    `--seed` *before* spawning workers — determinism must survive.
  - Add a `--workers` CLI flag (default: auto) to cap concurrency for
    constrained hosts.
- [ ] Parallelise per-person VCF writes
  - The loop at `cli.py:552` is the hot path for large `--n`.
  - Use `ProcessPoolExecutor` with the default fork start method on
    Linux so `cohort_sites` is shared copy-on-write — do not pickle
    cohort_sites per task.
  - Each worker derives its own `random.Random` from a per-person seed
    so output stays bit-identical to the serial run when the same
    `--seed` is passed.
  - Reuse the `--workers` flag added above.
- [ ] Stream straight into `bgzip -c`
  - Drop the plain `.vcf` intermediate in `writer.py:81-167`.
  - Open `Popen(["bgzip", "-c"], stdin=PIPE, stdout=open(out, "wb"))`,
    write records into stdin, close, then run `tabix` against the
    final `.vcf.gz`.
  - Saves one disk pass per person and removes a fork/exec.
- [ ] **Tests** — extend `tests/test_writer.py` (or add a new module) to
  cover: deterministic output across `--workers=1` vs `--workers=N`
  given the same seed; correctness of the bgzip-pipe path (round-trip
  via `bcftools view`).
- [ ] **Docs** — `README.md` Performance section, `TUTORIAL.md` §9
  ("Performance and scaling"), and the `--workers` `--help` text.

---

## Phase 2 — overlap I/O with compute

**Goal:** make ClinVar / rsID / COSMIC loaders run while msprime is
still simulating, since they're I/O-bound on bcftools subprocesses.

**Branch:** `perf/phase2-overlap-loaders`

- [ ] Run overlay loaders concurrently with simulation
  - `load_clinvar_index` (`cli.py:431`),
    `load_rsid_pool` (`cli.py:466`), and `load_cosmic_records`
    (`cli.py:492`) are all bcftools subprocess + I/O. They release
    the GIL.
  - Submit all three to a `ThreadPoolExecutor` *before*
    `simulate_cohort` runs. Block on the futures only at the point
    `cohort_sites` actually needs them.
  - Skip COSMIC submission unless `--somatic` is set, to preserve the
    "registration-gated, never auto-fetch" guarantee.
- [ ] **Tests** — add a regression test that the overlay-stats output
  is identical with and without the threaded prefetch, given a fixed
  seed.
- [ ] **Docs** — note in `IMPLEMENTATION_PLAN.md` and `README.md` that
  overlays prefetch in parallel; no user-visible flag change.

---

## Phase 3 — memory, medium risk

**Goal:** drop the per-site `list[str]` representation that dominates
RAM at large `n × n_sites`. Big peak-RAM drop, also speeds up
downstream loops.

**Branch:** `perf/phase3-genotype-matrix`

- [ ] Replace `gts: list[str]` with a numpy uint8 representation
  - Today: every site dict carries a list of `"0|1"` strings of length
    `n_people`. At `n=500, n_sites=50_000` this is ~25M Python strings
    (~50 B each → ~1.25 GB before measuring overhead).
  - Switch to `numpy.ndarray[uint8]` of shape `(2*n_people,)` per site,
    or — preferred — a single cohort-wide matrix
    `(n_sites, 2*n_people)` that all downstream code indexes.
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
- [ ] **Tests** — every test that constructs a site dict by hand needs
  updating; the genotype-matrix refactor is the riskiest item in this
  plan, so plan extra coverage of:
    - serial vs parallel (Phase 1 still passes deterministic check);
    - `person_records_from_cohort` generator yields the same records
      in the same order as the old list-returning version;
    - SFS histogram and overlay-stats numbers are byte-for-byte
      unchanged.
- [ ] **Docs** — note the new in-memory representation in
  `IMPLEMENTATION_PLAN.md` (architecture section), and call out the
  RAM win in `README.md` Performance and `TUTORIAL.md` §9.

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

## Open questions to resolve before starting

These don't block Phase 1, but the answers shape Phase 3 and the
choice of mechanism in Phase 1b.

- Which dimension does the user typically scale — `--n` or
  `--chromosomes`/`--chr-length-mb`?
- Single-machine only, or do we ever need to distribute across hosts?
  (Multi-host changes Phase 1b away from `ProcessPoolExecutor` toward
  a chunked work-queue.)

---

## Resuming after an interruption

1. `git branch --list 'perf/*'` — see which phase has an in-flight
   branch.
2. Open the matching section above and continue from the first
   unticked box.
3. If no `perf/*` branch exists, the next phase to start is the first
   one with any unticked boxes.
