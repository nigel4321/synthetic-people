# Performance budgets

Performance & memory guardrails for `synthetic_people`. The point
of this document is to make the targets visible and reviewable, not
buried in spike results — a number you can read in a code review is
a number you can defend against in a code review.

If you're changing a budget here, you should also be changing the
matching constant in `tests/test_performance_budgets.py` in the
same commit, and you should be able to explain *why* in the commit
message. A regression that doubles RAM is not a problem you fix by
doubling the budget.

---

## Why this exists

Phase 5d.1 (streaming-cohort + Arrow intermediate + chunked
simulation) spent ten PRs designing around the memory ceiling that
WGS-scale runs hit on 32 GB hosts. That work is currently
unprotected: no test fails today if a future PR silently
re-materialises a cohort-wide list, and the regression would only
surface as an OOM in production weeks later.

The tests in `tests/test_performance_budgets.py` codify the budgets
below and run on every PR.

---

## Canary scenario

| Property | Value |
|---|---|
| `--n` | 20 people |
| `--chromosomes` | 22 |
| `--chr-length-mb` | 1.0 |
| `--demo-model` | none (constant-Ne msprime path; no stdpopsim catalogue) |
| `--workers` | 1 |
| overlays | all disabled (rsids, clinvar, cosmic, SVs, error model) |
| `--seed` | 4242 |

Chosen so each per-mode test runs in ~12 s while still exercising
the full cohort backbone. Overlays are disabled because the
guardrail is about the cohort-write pipeline, not about overlay
density tradeoffs.

---

## Per-mode peak-RSS budgets

Absolute peak resident-set size of a fresh Python subprocess
invoking the cli on the canary scenario. Subprocess isolation is
required: in-process measurement is contaminated by sticky
allocations from prior tests in the same interpreter
(arrow-streaming measured 190 MB delta in isolation but 315 MB
when run after sites_list and arrow in the same process).

| `--cohort-mode` | Observed (32 GB Linux host) | Budget | Headroom |
|---|---|---|---|
| `sites_list` | ~312 MB | **450 MB** | 44 % |
| `arrow` | ~341 MB | **500 MB** | 47 % |
| `arrow-streaming` | ~339 MB | **500 MB** | 47 % |

### Why are the three budgets so close?

At canary scale the dominant cost is **imports** (pyarrow + msprime
+ stdpopsim) — not the cohort-state intermediate the streaming
architecture optimises. The streaming advantage is real but only
visible at WGS scale (Phase 5d.1's memprof28 measured ~17 GB
materialised-parent peak at n=3000 × 70 Mb; the streaming path
keeps that well below 4 GB at the same scale).

The PR-time canary intentionally trades visibility of the
streaming-vs-materialised ratio for fast CI execution. A nightly
WGS-scale memprof would surface that ratio — tracked as deferred
work below.

### Headroom rationale

40–47 % headroom over the observed peak is what we trade between:

- **Too tight** — false positives when CI runners are under
  memory pressure or when a transient import order changes
  baseline by ~30 MB.
- **Too loose** — fails to catch a ~100 MB regression that would
  matter at scale.

If false positives become routine, the next move is *not* to widen
the budget but to switch to relative budgets via a stored baseline
(pytest-benchmark style). The current absolute-budget approach is
deliberately simple — invest in the more elaborate infra only when
it pays for itself.

---

## Streaming-shape invariants

These are structural guards independent of measured RSS. They
catch the most common quick-fix regression pattern: "the streaming
iterator's test is failing, let me just collect it into a list."

| Invariant | Test |
|---|---|
| `stream_cohort_sites` is a generator function | `tests/test_performance_budgets.py::StreamingShapeInvariantsTest::test_stream_cohort_sites_is_generator_function` |
| `simulate_cohort_ts_iter` is a generator function | same file, `test_simulate_cohort_ts_iter_is_generator_function` |
| Streaming yields incrementally (partial-consume pattern works) | same file, `test_stream_cohort_sites_yields_incrementally` |

If `inspect.isgeneratorfunction` returns False on either of those
two entry points, **the streaming guarantee has been
silently undone** — anything that calls them and expects a
generator is now collecting a full list in memory. The fix is to
revert the change that converted the `yield` to a `return`, not to
delete the test.

---

## What this catches

- A new dependency import that adds tens of MB of constant
  allocation.
- A forgotten cohort-wide list (sites accumulator, person records,
  bcf row buffer) re-introduced into the parent process.
- A leaked buffer that doesn't drop between chromosomes.
- A `yield` accidentally changed to `return list(...)` during a
  refactor.

## What this does NOT catch

- **Subtle O(n)-scaling regressions** whose effect at n=20 is
  below the ~100 MB noise floor. A 50 KB-per-sample regression
  costs ~1 MB at canary scale but 150 MB at WGS-n3000 and 50 GB at
  n=1M. Needs a nightly WGS-scale canary. *(Tracked as deferred
  work — see below.)*
- **Wall-clock regressions.** Memory ≠ time; a 10× slowdown that
  doesn't allocate more is invisible to RSS budgets. Needs
  pytest-benchmark or ASV. *(Not a priority unless a runtime
  regression actually hits us.)*
- **Auto-picker calibration drift.** If `_estimate_materialised_parent_peak_bytes`
  drifts from reality, the picker silently mis-routes (e.g.,
  chooses arrow-streaming when sites_list would fit). A separate
  calibration test would assert predicted ≈ observed peak within
  tolerance. *(Tracked as deferred work.)*
- **Per-worker peak RSS.** `--workers 1` is what the canary tests;
  worker fan-out has its own per-process budget separate from
  parent peak. Today only the parent is gated.

---

## Known scaling ceiling — `arrow-streaming` at n ≥ ~1M

**Empirical data point, 2026-05-12.** A user-driven run at
**n=1 000 000**, full WGS (chr 1-22 + X @ 70 Mb per chrom),
`--cohort-mode arrow-streaming`, `--mode per-person` was killed
at **63 GB parent RSS** ~115 min into chrom 1's streaming-write
pass — before any worker fan-out (`children_rss_mb` was 0
throughout). memprof phase trail:

| t | parent RSS | phase |
|---|---|---|
| 16 s | 718 MB | streaming sim start (23 chroms) |
| 901 s | 2.6 GB | chrom 1 TS ready (1 081 579 sites) |
| 6905 s | **63.3 GB** | killed (no further phase mark reached) |

### Root cause

The streaming heap in `_stream_cohort_pass2` (see
`coalescent.py:611-619` for the design-time memory model) holds
**full site dicts** whose `carriers` field is a Python
`list[tuple[int, int]]` — one tuple per non-zero haplotype.
Per-site cost scales linearly with n:

| AF | n=3 000 | n=100 000 | n=1 000 000 |
|---|---|---|---|
| singleton | ~80 B | ~80 B | ~80 B |
| 5 % common | ~24 KB | ~800 KB | **~8 MB** |
| 30 % common | ~144 KB | ~4.8 MB | **~48 MB** |

The buffer depth is bounded by overlay-injection-position pressure
(~`O(sqrt(N_inject))` if positions are uniform, larger when they
cluster). At the canonical 0.2 rsid density × real dbSNP positions
on chr1, the buffer plausibly holds hundreds of common-AF sites
concurrently — fine at n=3000 (hundreds of MB), fatal at n=1M
(hundreds of GB).

### Why this isn't caught by the PR-time canary

The canary at n=20 has a per-site cost of ~1.6 KB even for common
AF — six orders of magnitude smaller than n=1M. The architectural
scaling problem is invisible until n is large enough for
per-site cost to dominate, which is exactly why a nightly
WGS-scale memprof CI is on the deferred list.

### Supported envelope (today)

- **n ≤ ~100 000**, full WGS, `--mode cohort` — comfortable.
- **n ≤ ~500 000**, full WGS, `--mode cohort` — feasible with
  reduced overlay density (rsid ≤ 0.05) on a 64 GB host.
- **n ≥ ~1 000 000** — needs the carriers-packing fix below.

### Fix sketch (not yet implemented)

Pack `carriers` as a numpy array of haplotype indices rather than
a list of `(int, int)` tuples. For biallelic sites (the vast
majority) the allele field is always 1 and can be elided.
`np.array(hap_indices, dtype=np.int32)` ≈ 4 B per carrier, vs
~80 B for the Python tuple representation — a **~20× memory
reduction** across all carriers-holding code paths (streaming
heap, Arrow writer batches, materialised sites list). Should
push the streaming ceiling from ~n=500K to ~n=5M without
architectural change. Tracked in deferred-guardrails below as
**carriers packing**.

---

## Updating a budget

The right reasons to widen a budget:

- A legitimately new feature has been added that genuinely needs
  the memory (e.g., loading the GRCh38 reference FASTA — adds
  ~50 MB across all modes, expected and intentional).
- The test runner platform genuinely changed (e.g., we moved to a
  bigger CI runner and Python imports themselves got heavier).

The wrong reason to widen a budget:

- A PR introduces a regression that pushes RSS over the budget,
  and the easy fix is to bump the budget rather than understand
  what got allocated.

When you do bump it, update the constant in
`tests/test_performance_budgets.py::PEAK_RSS_BUDGET_MB` AND this
document's per-mode table in the same commit. The commit message
should explain what the new code legitimately needs.

---

## Deferred guardrails

These would tighten the safety net further but each carries its
own infrastructure cost; ship when the cost is justified by an
actual regression they would have caught.

| Guardrail / fix | Catches | Why not yet / justification |
|---|---|---|
| **Carriers packing** (numpy array, drop allele field for biallelics) | The n=1M scaling ceiling documented above — direct ~20× memory reduction across streaming heap + Arrow writer + materialised sites | We now have a concrete regression case (2026-05-12 user run). Right next move after this doc update. Estimate: 1 day implementation + 1 day touching consumers (Arrow writer, BCF emission). |
| Nightly WGS-scale memprof CI | Subtle O(n) re-materialisation, slow-bleed drift, scaling-ceiling regressions | Needs a dedicated runner with predictable RAM; ongoing cost. The 2026-05-12 incident is exactly the case that justifies this. Worth ~1 day of plumbing now that we have a regression case. |
| Auto-picker calibration test | Drift between `_estimate_materialised_parent_peak_bytes` and reality | The estimator is itself an approximation; calibration would need a tolerance band that's hard to set without observing real drift first. |
| pytest-benchmark wall-clock regression | Time-only regressions (no RSS impact) | Variance on shared CI runners makes time-based gates noisy; only worth wiring up if a real time regression hits us. |
| Per-worker peak RSS | Worker-side leaks under multi-process fan-out | Today only the parent is gated; per-worker is a smaller scaling concern (workers spawn fresh and reap promptly). |

---

## References

- [`PERFORMANCE_PLAN.md`](PERFORMANCE_PLAN.md) — the design document
  behind the streaming-cohort architecture this doc is gating.
- [`tests/test_performance_budgets.py`](tests/test_performance_budgets.py)
  — the tests this document is the source-of-truth for.
- [`DATA_QUALITY_ASSESSMENT.md`](DATA_QUALITY_ASSESSMENT.md) §6 — what
  performance compromises Phase 5d.1 made and why removing them is
  now feasible.
