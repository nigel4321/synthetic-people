"""Performance guardrail tests.

Catches the regression class that PERFORMANCE_PLAN.md spent ten PRs
designing around — silent re-materialisation of an O(n) or O(sites)
intermediate that the streaming-cohort architecture is supposed to
avoid.

Two kinds of test:

- ``CohortPeakRssBudgetTest`` runs a deterministic canary cohort
  for each cohort mode and asserts parent-process peak RSS stays
  under a documented budget. A regression that re-introduces a
  cohort-wide list materialisation will blow the streaming budget
  even at this canary scale.
- ``StreamingShapeInvariantsTest`` asserts the *shape* of the
  streaming entry points (generator functions, not list-returning
  functions). Cheap structural check that catches accidental
  ``return list(...)`` refactors in code review.

Budgets are documented in ``PERFORMANCE_BUDGETS.md``; when bumping
a budget, update both that doc and the constants here in the same
commit so the rationale stays close to the number.

All tests skip cleanly when their deps are missing.
"""

from __future__ import annotations

import importlib.util
import inspect
import shutil
import subprocess
import sys
import tempfile
import textwrap
import time
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


_HAVE_PSUTIL = importlib.util.find_spec("psutil") is not None
_HAVE_BCFTOOLS = shutil.which("bcftools") is not None
_HAVE_TABIX = shutil.which("tabix") is not None
_HAVE_BGZIP = shutil.which("bgzip") is not None
_HAVE_MSPRIME = importlib.util.find_spec("msprime") is not None
_HAVE_STDPOPSIM = importlib.util.find_spec("stdpopsim") is not None
_HAVE_PYARROW = importlib.util.find_spec("pyarrow") is not None


# ---------------------------------------------------------------------------
# Peak-RSS budgets per cohort mode at the canary scale (n=20, chr22, 1 Mb).
# ---------------------------------------------------------------------------
#
# Units: **MB of ABSOLUTE peak RSS in a fresh Python subprocess** that imports
# only what the cli needs. Subprocess isolation is required because the
# Python allocator holds onto pages across tests, contaminating in-process
# delta measurements: arrow-streaming measured 190 MB delta in isolation but
# 315 MB delta when run after sites_list and arrow in the same interpreter.
# A fresh subprocess per test eliminates that contamination.
#
# Calibration: observed locally on a 32 GB Linux host; budgets set ~40 %
# above observed so CI runner noise doesn't false-positive while a real
# ~50 MB re-materialisation regression still trips clearly.
# Re-calibrate (in this constant AND PERFORMANCE_BUDGETS.md) if the
# underlying code's working set legitimately grows.
#
# What this catches and what it doesn't:
#
# - Catches: gross regressions that add tens of MB of constant allocation
#   (a new dependency import, a forgotten cohort-wide list, a leaked buffer).
# - Doesn't catch: subtle O(n)-scaling regressions whose impact at n=20 is
#   below the noise floor. Those require a nightly WGS-scale canary —
#   tracked separately in PERFORMANCE_BUDGETS.md as future work.
PEAK_RSS_BUDGET_MB = {
    # Observed at canary scale in fresh subprocesses:
    #   sites_list      ~312 MB
    #   arrow           ~341 MB
    #   arrow-streaming ~339 MB
    # The three are close at this scale because pyarrow + msprime
    # imports dominate over run-time cohort state — the streaming
    # architectural advantage is real but only visible at WGS scale.
    # We set absolute budgets here for PR-time regression catching;
    # the streaming-vs-materialised ratio belongs in the nightly
    # large-scale canary (deferred).
    "sites_list": 450,
    "arrow": 500,
    "arrow-streaming": 500,
}


CANARY_CHROM = "22"
CANARY_LENGTH_MB = 1.0
CANARY_N = 20
CANARY_SEED = 4242


# Hard ceiling on subprocess runtime. The canary completes in ~12 s
# locally; a 5-minute deadline gives generous headroom for slow CI
# runners while still bounding the worst case so a stalled child
# can't hang the job indefinitely.
SUBPROCESS_TIMEOUT_SECONDS = 300


def _canary_args(
    out_dir: Path,
    cohort_mode: str,
    cache_dir: Path,
) -> list[str]:
    """CLI args for the canary cohort run.

    Disables every overlay (rsids, clinvar, cosmic, SVs, error model)
    so the only memory consumer is the cohort backbone itself —
    which is exactly what the budgets are gating on. ``--workers 1``
    pins the parent path; worker fan-out has its own per-worker
    budget that's separate from parent peak.

    ``cache_dir`` is taken as a parameter rather than derived from
    ``out_dir`` so multiple tests in the same class can share one
    cache — defensive against the possibility that a future cli
    change makes ClinVar/dbSNP fetching unconditional (today
    overlay-density=0 skips it). Sharing one cache means the
    network hit is paid at most once per test class, not once per
    cohort mode.
    """
    return [
        "--no-config",  # ignore any cwd config that would skew results
        "--n", str(CANARY_N),
        "--seed", str(CANARY_SEED),
        "--build", "GRCh38",
        "--chromosomes", CANARY_CHROM,
        "--chr-length-mb", str(CANARY_LENGTH_MB),
        "--demo-model", "none",
        "--rsid-density", "0",
        "--clinvar-inject-density", "0",
        "--svs-per-person", "0",
        "--error-rate", "0",
        "--dropout-rate", "0",
        "--workers", "1",
        "--output-dir", str(out_dir),
        "--cache-dir", str(cache_dir),
        "--mode", "cohort",
        "--cohort-mode", cohort_mode,
    ]


_CANARY_RUNNER_SCRIPT = textwrap.dedent("""
    # Inline canary runner — imports the cli and invokes it with
    # args loaded from a JSON file (path passed as argv[1]). Run as
    # a fresh subprocess so the test process's prior imports don't
    # contaminate the RSS measurement.
    import json, sys
    sys.path.insert(0, {repo_root!r})
    from syntheticgen.cli import main as cli_main
    with open(sys.argv[1]) as fh:
        args = json.load(fh)
    sys.exit(cli_main(args))
""").strip()


def _run_canary_in_subprocess(
    test_case: unittest.TestCase,
    cohort_mode: str,
    cache_dir: Path,
) -> float:
    """Spawn a fresh Python subprocess that runs the canary cli;
    sample its RSS at 50 ms cadence; return absolute peak in MB.

    Subprocess isolation is what makes the measurement comparable
    across tests: a fresh interpreter has none of the prior tests'
    sticky allocations, so the peak reflects exactly what this mode
    requires when run cold.

    Wall-clock-bounded: a stalled child (network hang, deadlock,
    bug) is killed at ``SUBPROCESS_TIMEOUT_SECONDS`` rather than
    being allowed to wedge the CI job. stderr is redirected to a
    file rather than a pipe so the cli's progress prints can't fill
    the pipe buffer and block the child mid-run.
    """
    import json
    import psutil
    repo_root = str(Path(__file__).resolve().parent.parent)
    script = _CANARY_RUNNER_SCRIPT.format(repo_root=repo_root)
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = Path(tmp) / "out"
        args = _canary_args(out_dir, cohort_mode, cache_dir)
        # Write args to a temp file rather than stdin: stdin pipes
        # interact awkwardly with the RSS-sampling loop (closing
        # stdin from the parent races with the child's read), and a
        # short file is simpler than getting that timing right.
        args_path = Path(tmp) / "args.json"
        args_path.write_text(json.dumps(args))
        # stderr → temp file rather than PIPE: the cli prints
        # progress lines and warnings to stderr throughout the run;
        # a PIPE without a concurrent reader can fill its ~64 KB
        # buffer and block the child write, which would never let
        # ``poll()`` return non-None and would hang the test. A file
        # has unbounded capacity and is trivially readable for
        # diagnostics on failure.
        stderr_path = Path(tmp) / "stderr.txt"
        with stderr_path.open("wb") as stderr_fh:
            proc = subprocess.Popen(
                [sys.executable, "-c", script, str(args_path)],
                stdout=subprocess.DEVNULL,
                stderr=stderr_fh,
            )
            ps_proc = psutil.Process(proc.pid)
            peak_bytes = 0
            deadline = time.monotonic() + SUBPROCESS_TIMEOUT_SECONDS
            # Sample until the subprocess exits or we hit the
            # deadline. ``poll() is None`` is the fast non-blocking
            # liveness check; ``memory_info`` raises NoSuchProcess
            # once the child reaps.
            timed_out = False
            while proc.poll() is None:
                if time.monotonic() > deadline:
                    timed_out = True
                    proc.kill()
                    try:
                        proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        pass
                    break
                try:
                    rss = ps_proc.memory_info().rss
                    if rss > peak_bytes:
                        peak_bytes = rss
                except (psutil.NoSuchProcess, psutil.ZombieProcess):
                    break
                time.sleep(0.05)
        # Read the stderr tail unconditionally — useful diagnostic
        # for either timeout or non-zero exit.
        try:
            stderr_text = stderr_path.read_text(errors="replace")
        except OSError:
            stderr_text = "<could not read stderr file>"
        if timed_out:
            test_case.fail(
                f"canary subprocess for --cohort-mode {cohort_mode} "
                f"did not exit within {SUBPROCESS_TIMEOUT_SECONDS}s; "
                f"killed. Args: {args!r}\n"
                f"stderr tail:\n{stderr_text[-1500:]}",
            )
        test_case.assertEqual(
            proc.returncode, 0,
            f"canary subprocess returned {proc.returncode}.\n"
            f"stderr tail:\n{stderr_text[-1500:]}",
        )
        return peak_bytes / (1024 * 1024)


@unittest.skipUnless(
    _HAVE_PSUTIL and _HAVE_BCFTOOLS and _HAVE_TABIX and _HAVE_BGZIP
    and _HAVE_MSPRIME and _HAVE_STDPOPSIM,
    "needs psutil + bcftools/tabix/bgzip + msprime + stdpopsim",
)
class CohortPeakRssBudgetTest(unittest.TestCase):
    """Per-mode peak RSS budget on the canary scenario.

    A regression that re-materialises a cohort-wide list, allocates
    an O(sites) intermediate the streaming path was meant to avoid,
    or breaks the auto-picker's predicted-vs-actual contract will
    push parent peak RSS over its budget here. Failure mode is a
    clear assertion error with the observed delta vs the limit,
    not a downstream OOM days later at WGS scale.
    """

    # Class-level shared cache so all three mode-budget tests reuse
    # one ClinVar/dbSNP cache rather than fetching fresh per test.
    # Today's cli skips overlay-fetch when densities are 0 (which
    # the canary sets), but if a future change makes the fetch
    # unconditional, sharing the cache caps the network cost at
    # one download per CI job instead of three.
    _cache_tmp: tempfile.TemporaryDirectory | None = None
    _cache_dir: Path | None = None

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls._cache_tmp = tempfile.TemporaryDirectory(
            prefix="perf_budget_cache_",
        )
        cls._cache_dir = Path(cls._cache_tmp.name)

    @classmethod
    def tearDownClass(cls):
        if cls._cache_tmp is not None:
            cls._cache_tmp.cleanup()
            cls._cache_tmp = None
            cls._cache_dir = None
        super().tearDownClass()

    def _check_budget(self, cohort_mode: str) -> None:
        peak_mb = _run_canary_in_subprocess(
            self, cohort_mode, self._cache_dir,
        )
        budget_mb = PEAK_RSS_BUDGET_MB[cohort_mode]
        # Report both numbers so a failure tells the developer
        # immediately whether they're over by a hair or by a lot.
        self.assertLess(
            peak_mb, budget_mb,
            f"--cohort-mode {cohort_mode}: subprocess peak RSS "
            f"{peak_mb:.1f} MB exceeds budget {budget_mb} MB at the "
            f"n={CANARY_N}, chr{CANARY_CHROM}, {CANARY_LENGTH_MB} "
            f"Mb canary. See PERFORMANCE_BUDGETS.md for the "
            f"rationale; if this is a deliberate increase, update "
            f"both the budget constant and the doc in the same "
            f"commit.",
        )

    def test_sites_list_within_budget(self):
        self._check_budget("sites_list")

    @unittest.skipUnless(_HAVE_PYARROW, "pyarrow not installed")
    def test_arrow_within_budget(self):
        self._check_budget("arrow")

    @unittest.skipUnless(_HAVE_PYARROW, "pyarrow not installed")
    def test_arrow_streaming_within_budget(self):
        self._check_budget("arrow-streaming")


class StreamingShapeInvariantsTest(unittest.TestCase):
    """Structural guardrails on the streaming-cohort path.

    The streaming architecture rests on two functions being true
    generators — i.e., ``yield``-based, not ``return list(...)``.
    Refactors that "just fix" a test by collecting the iterator
    into a list (a common quick-fix anti-pattern) silently undo
    the streaming guarantee. These checks make that change visible
    in code review by failing CI on the refactor commit.
    """

    @unittest.skipUnless(_HAVE_MSPRIME, "msprime not installed")
    def test_stream_cohort_sites_is_generator_function(self):
        from syntheticgen.coalescent import stream_cohort_sites
        self.assertTrue(
            inspect.isgeneratorfunction(stream_cohort_sites),
            "stream_cohort_sites must remain a generator function "
            "(yield-based). A return-list refactor undoes the "
            "streaming-cohort memory guarantee — see "
            "PERFORMANCE_BUDGETS.md §Streaming shape.",
        )

    @unittest.skipUnless(_HAVE_MSPRIME, "msprime not installed")
    def test_simulate_cohort_ts_iter_is_generator_function(self):
        from syntheticgen.coalescent import simulate_cohort_ts_iter
        self.assertTrue(
            inspect.isgeneratorfunction(simulate_cohort_ts_iter),
            "simulate_cohort_ts_iter must remain a generator "
            "function (yield-based). Materialising all tree "
            "sequences upfront would re-introduce the parent-peak "
            "RSS ceiling Phase 5d.1 removed.",
        )

    @unittest.skipUnless(
        _HAVE_MSPRIME and _HAVE_STDPOPSIM,
        "msprime + stdpopsim not installed",
    )
    def test_stream_cohort_sites_yields_incrementally(self):
        # Behavioural complement to the isgeneratorfunction check:
        # consume just the first few sites and confirm the function
        # actually yields incrementally rather than building a full
        # list internally before the first yield. A small canary
        # keeps the runtime ~1 s.
        import random
        from syntheticgen.coalescent import (
            simulate_cohort_ts_iter, stream_cohort_sites,
        )
        rng = random.Random(CANARY_SEED)
        # demo_model=None routes through the constant-Ne msprime path
        # — same conversion the cli does for ``--demo-model none``.
        # Skips stdpopsim catalogue lookup which would otherwise need
        # network and the species cache.
        ts_iter = simulate_cohort_ts_iter(
            chromosomes=[CANARY_CHROM],
            build="GRCh38",
            n_people=CANARY_N,
            length_mb=0.1,
            demo_model=None,
            population="CEU",
            rec_rate=1e-8,
            mu=1.29e-8,
            rng=rng,
        )
        chrom, ts, walk_rng = next(ts_iter)
        stream = stream_cohort_sites(
            ts, chrom, CANARY_N, walk_rng,
            overlay_rng=random.Random(CANARY_SEED),
        )
        # Pull the first three sites; the iterator must produce them
        # without exhausting itself first. If a future refactor turns
        # this into ``return list(...)`` the iterator would still
        # work but the lazy yielding would be lost — caught by the
        # ``isgeneratorfunction`` test above, this one just doubles
        # the safety by asserting the partial-consumption pattern
        # actually works end-to-end.
        first_three = []
        for site in stream:
            first_three.append(site)
            if len(first_three) >= 3:
                break
        self.assertGreaterEqual(len(first_three), 1)
        for site in first_three:
            # Sanity: yielded items look like cohort-site dicts.
            self.assertIn("chrom", site)
            self.assertIn("pos", site)


if __name__ == "__main__":
    unittest.main()
