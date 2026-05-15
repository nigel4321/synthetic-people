"""Phase 5b2 resume contract tests.

The cohort.meta.json schema records the resume-relevant params plus
the sample IDs, per-person seeds, per-chromosome overlay seeds, and
the list of chromosomes whose cohort BCF has finished writing. This
file covers the contract:

- A fresh run writes a meta.json before simulating, with empty
  ``completed_chromosomes`` and freshly-derived seeds.
- Each chromosome's completion appends to ``completed_chromosomes``
  atomically (rename-based write so a SIGINT mid-flush leaves the
  prior version intact).
- A resume that finds matching params skips already-complete
  chromosomes and re-uses the persisted seeds.
- A resume with mismatched params surfaces a clear error (and
  ``--no-resume`` overrides it).

Heavy paths (full streamed pipeline) gate on bcftools/tabix/bgzip +
msprime + stdpopsim. Pure-Python tests of the resume module itself
run on any host.
"""

from __future__ import annotations

import importlib.util
import json
import os
import random
import shutil
import sys
import tempfile
import types
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen import cli as cli_module
from syntheticgen.resume import (
    Resume,
    ResumeMismatch,
    load_or_create_meta,
)
from tests._shared_cache import SHARED_TEST_CACHE_DIR


_HAVE_BCFTOOLS = shutil.which("bcftools") is not None
_HAVE_MSPRIME = importlib.util.find_spec("msprime") is not None
_HAVE_STDPOPSIM = importlib.util.find_spec("stdpopsim") is not None


def _args_ns(**overrides):
    """Build a minimal argparse-style Namespace for the resume helper."""
    base = dict(
        seed=42, n=3, build="GRCh38",
        chr_length_mb=0.2, demo_model="none", population="CEU",
        rec_rate=1e-8, mu=1.29e-8,
        # M13.1: load_or_create_meta uses ``args.male_fraction`` to
        # draw per-person sexes alongside person_seeds. Match the
        # production default so the fixture stays representative.
        male_fraction=0.5,
    )
    base.update(overrides)
    return types.SimpleNamespace(**base)


class ResumeFreshStartTest(unittest.TestCase):
    """A fresh run with no existing cohort.meta.json should derive
    new seeds from the rng and persist them on disk."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="resume_fresh_"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_fresh_writes_meta_json(self):
        rng = random.Random(42)
        r = load_or_create_meta(
            _args_ns(), ["22"], self.tmp, rng,
        )
        meta = self.tmp / "cohort.meta.json"
        self.assertTrue(meta.is_file(),
                        f"expected meta at {meta}")
        payload = json.loads(meta.read_text())
        self.assertEqual(payload["completed_chromosomes"], [])
        self.assertEqual(len(payload["samples"]), 3)
        self.assertEqual(len(payload["person_seeds"]), 3)
        # M13.1: sexes persisted alongside samples, same cardinality.
        self.assertEqual(len(payload["sexes"]), 3)
        self.assertTrue(all(s in ("m", "f") for s in payload["sexes"]))
        self.assertEqual(set(payload["overlay_seeds"]), {"22"})
        # Returned object reflects the persisted state.
        self.assertEqual(r.samples, payload["samples"])
        self.assertEqual(r.sexes, payload["sexes"])

    def test_seeds_are_deterministic_at_fixed_seed(self):
        a = load_or_create_meta(_args_ns(), ["22"],
                                self.tmp / "a", random.Random(42))
        b = load_or_create_meta(_args_ns(), ["22"],
                                self.tmp / "b", random.Random(42))
        self.assertEqual(a.samples, b.samples)
        self.assertEqual(a.person_seeds, b.person_seeds)
        # M13.1: sexes must also be deterministic at fixed seed —
        # they come from the same master rng as samples/person_seeds.
        self.assertEqual(a.sexes, b.sexes)
        self.assertEqual(a.overlay_seeds, b.overlay_seeds)

    def test_male_fraction_zero_yields_all_female(self):
        # M13.1 sanity: male_fraction=0.0 must produce all "f".
        # Defends against a polarity flip (e.g. someone interpreting
        # the field as P(female) instead of P(male)) — a real risk
        # because "sex ratio" is ambiguous in nature; the field name
        # `male_fraction` and these tests pin the convention.
        r = load_or_create_meta(
            _args_ns(male_fraction=0.0), ["22"],
            self.tmp, random.Random(42),
        )
        self.assertEqual(r.sexes, ["f", "f", "f"])

    def test_male_fraction_one_yields_all_male(self):
        r = load_or_create_meta(
            _args_ns(male_fraction=1.0), ["22"],
            self.tmp, random.Random(42),
        )
        self.assertEqual(r.sexes, ["m", "m", "m"])

    def test_sex_draws_do_not_advance_master_rng(self):
        # PR #95 review (Copilot): M13.1's load-bearing claim is that
        # adding sex assignment doesn't change the simulator's output
        # at a fixed seed. That contract holds iff sex draws DON'T
        # consume the master rng — otherwise overlay_seeds and every
        # downstream rng consumer would shift relative to pre-M13.1.
        #
        # This test fingerprints the master rng's state after a
        # load_or_create_meta call. The state must be byte-identical
        # to drawing only samples + person_seeds + overlay_seeds
        # (i.e. NOT advanced by the sex draws).
        rng_with_sex = random.Random(42)
        load_or_create_meta(
            _args_ns(male_fraction=0.5), ["22"],
            self.tmp / "a", rng_with_sex,
        )

        # Replay the exact same rng draws that load_or_create_meta
        # makes from the master rng, EXCLUDING sex draws.
        from syntheticgen.background import draw_sample_ids
        rng_no_sex = random.Random(42)
        _ = draw_sample_ids(3, rng_no_sex)
        _ = [rng_no_sex.randint(1, 2**31 - 1) for _ in range(3)]
        _ = rng_no_sex.randint(1, 2**31 - 1)  # overlay_seeds["22"]

        # Both rngs must now produce the same next value — proves the
        # sex draws didn't consume from the master rng.
        self.assertEqual(rng_with_sex.random(), rng_no_sex.random())


class ResumeMatchingParamsTest(unittest.TestCase):
    """An existing meta.json whose params match should be re-used and
    its completed_chromosomes carried through."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="resume_match_"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_resume_matches_returns_persisted_state(self):
        # Initial run: fresh meta + mark chr22 as complete.
        first = load_or_create_meta(
            _args_ns(), ["22", "21"], self.tmp, random.Random(42))
        first.mark_chromosome_done("22")
        # Resume: a different rng must not be consumed (seeds come
        # from disk). We pass a one-trick rng that explodes if drawn
        # from to enforce this.
        class ExplodingRng:
            def randint(self, a, b):
                raise AssertionError(
                    "resume must not redraw seeds from rng")
        second = load_or_create_meta(
            _args_ns(), ["22", "21"], self.tmp, ExplodingRng())
        self.assertEqual(second.samples, first.samples)
        self.assertEqual(second.person_seeds, first.person_seeds)
        # M13.1: sexes also round-trip via the persisted meta.
        self.assertEqual(second.sexes, first.sexes)
        self.assertEqual(second.overlay_seeds, first.overlay_seeds)
        self.assertEqual(second.completed_chromosomes, ["22"])


class ResumeMismatchedParamsTest(unittest.TestCase):
    """Params that don't match the persisted record should surface a
    clear ``ResumeMismatch`` error and ``--no-resume`` should
    override."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="resume_mismatch_"))
        load_or_create_meta(
            _args_ns(seed=42), ["22"], self.tmp, random.Random(42))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_seed_mismatch_raises(self):
        with self.assertRaises(ResumeMismatch) as cm:
            load_or_create_meta(
                _args_ns(seed=999), ["22"], self.tmp, random.Random(0))
        self.assertIn("--no-resume", str(cm.exception))

    def test_chromosomes_mismatch_raises(self):
        with self.assertRaises(ResumeMismatch):
            load_or_create_meta(
                _args_ns(seed=42), ["22", "21"],
                self.tmp, random.Random(0))

    def test_male_fraction_mismatch_raises(self):
        # PR #95 review (Copilot): male_fraction is now part of the
        # resume identity. Changing it between runs must trigger
        # ResumeMismatch — without this, the persisted ``sexes`` would
        # silently be reused even though the new run intends a
        # different sex composition.
        with self.assertRaises(ResumeMismatch):
            load_or_create_meta(
                _args_ns(seed=42, male_fraction=0.8), ["22"],
                self.tmp, random.Random(0))

    def test_force_fresh_wipes_existing_state(self):
        # Pre-populate cohort/ with a sentinel file; --no-resume
        # should remove it alongside the old meta.json.
        sentinel = self.tmp / "cohort.chr22.bcf"
        sentinel.write_text("stale")
        r = load_or_create_meta(
            _args_ns(seed=999), ["22"], self.tmp,
            random.Random(0), force_fresh=True,
        )
        self.assertFalse(sentinel.exists(),
                         "force_fresh must wipe stale BCF residue")
        self.assertEqual(r.params["seed"], 999)


def _streamed_args(out_dir: Path, no_resume: bool = False) -> list:
    args = [
        "--n", "3",
        "--seed", "42",
        "--build", "GRCh38",
        "--chromosomes", "21,22",
        "--chr-length-mb", "0.2",
        "--demo-model", "none",
        "--rsid-density", "0",
        "--clinvar-inject-density", "0",
        "--svs-per-person", "0",
        "--error-rate", "0",
        "--dropout-rate", "0",
        "--workers", "1",
        "--output-dir", str(out_dir),
        # Share the ClinVar download across tests in this process —
        # see tests/_shared_cache.py.
        "--cache-dir", str(SHARED_TEST_CACHE_DIR),
        "--mode", "cohort",
        # M12: opt out of the auto-fetch (default-on behaviour
        # would download a 3 GB FASTA into the per-test cache_dir).
        "--no-reference-fasta",
    ]
    if no_resume:
        args.append("--no-resume")
    return args


@unittest.skipUnless(
    _HAVE_BCFTOOLS and _HAVE_MSPRIME and _HAVE_STDPOPSIM,
    "needs bcftools + msprime + stdpopsim",
)
class ResumeEndToEndTest(unittest.TestCase):
    """Run the streamed pipeline; delete one chrom's BCF; re-run; the
    surviving chrom's BCF stays untouched (resume worked) and the
    deleted one gets re-generated."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(prefix="resume_e2e_"))
        rc = cli_module.main(_streamed_args(cls.tmpdir))
        if rc != 0:
            raise RuntimeError(f"first run failed (rc={rc})")
        # Snapshot mtimes of both per-chrom BCFs.
        cls.chr21 = cls.tmpdir / "cohort" / "cohort.chr21.bcf"
        cls.chr22 = cls.tmpdir / "cohort" / "cohort.chr22.bcf"
        cls.mtime_chr22_before = cls.chr22.stat().st_mtime_ns
        # Simulate "interrupted before chr21 finished" by removing
        # chr21 from disk and dropping it from completed_chromosomes.
        cls.chr21.unlink()
        Path(str(cls.chr21) + ".csi").unlink(missing_ok=True)
        meta_path = cls.tmpdir / "cohort" / "cohort.meta.json"
        meta = json.loads(meta_path.read_text())
        meta["completed_chromosomes"] = [
            c for c in meta["completed_chromosomes"] if c != "21"
        ]
        meta_path.write_text(json.dumps(meta))
        # Re-run: chr22 should survive, chr21 should regenerate.
        rc = cli_module.main(_streamed_args(cls.tmpdir))
        if rc != 0:
            raise RuntimeError(f"resume run failed (rc={rc})")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_chr22_was_not_resimulated(self):
        # mtime preserved means resume reused the existing BCF rather
        # than re-running simulate_cohort_iter for chr22.
        self.assertEqual(
            self.chr22.stat().st_mtime_ns,
            self.mtime_chr22_before,
            msg="chr22 BCF was rewritten — resume should have skipped it",
        )

    def test_chr21_was_regenerated(self):
        self.assertTrue(self.chr21.is_file())
        self.assertTrue(Path(str(self.chr21) + ".csi").is_file())

    def test_no_resume_wipes_and_redoes_everything(self):
        before = self.chr22.stat().st_mtime_ns
        rc = cli_module.main(_streamed_args(self.tmpdir, no_resume=True))
        self.assertEqual(rc, 0)
        # chr22 should have been rewritten under --no-resume.
        self.assertNotEqual(
            self.chr22.stat().st_mtime_ns, before,
            msg="--no-resume should have rewritten chr22 BCF",
        )


if __name__ == "__main__":
    unittest.main()
