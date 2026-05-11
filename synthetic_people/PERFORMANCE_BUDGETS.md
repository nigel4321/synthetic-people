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

| Guardrail | Catches | Why not yet |
|---|---|---|
| Nightly WGS-scale memprof CI | Subtle O(n) re-materialisation, slow-bleed drift | Needs a dedicated runner with predictable RAM; ongoing cost. Worth ~1 day of plumbing once we have a regression case that justifies it. |
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
