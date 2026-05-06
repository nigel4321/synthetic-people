"""End-to-end --mode {per-person, cohort, both} smoke tests.

These tests run the actual CLI in-process against a tiny chr22
cohort and verify the right deliverables land on disk for each mode.
The runs are sized for CI: 3 people × 0.2 Mb × chr22 with the
``--demo-model none`` constant-Ne path (no stdpopsim model lookup,
no overlay downloads), so each invocation finishes in a few
seconds. Tests gate on bcftools/tabix/bgzip on PATH and on
msprime/stdpopsim importability — same gating pattern as
``test_coalescent.py``.

What we're guarding against:

- ``--mode per-person`` (the default) must remain byte-equivalent to
  today: no cohort/ directory should appear, manifest should not
  carry a ``cohort_bcfs`` field. A regression here would surprise
  every existing user.
- ``--mode cohort`` must skip per-person fan-out entirely and emit a
  cohort BCF whose per-sample columns match the cohort the
  simulator produced.
- ``--mode both`` must emit both deliverables in the same run, and
  the cohort BCF's record set must align with the per-person VCFs
  by sample id.
"""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen import cli as cli_module


_HAVE_BCFTOOLS = shutil.which("bcftools") is not None
_HAVE_TABIX = shutil.which("tabix") is not None
_HAVE_BGZIP = shutil.which("bgzip") is not None
_HAVE_MSPRIME = importlib.util.find_spec("msprime") is not None
_HAVE_STDPOPSIM = importlib.util.find_spec("stdpopsim") is not None


def _common_args(out_dir: Path, mode: str, n: int = 3) -> list:
    return [
        "--n", str(n),
        "--seed", "42",
        "--build", "GRCh38",
        "--chromosomes", "22",
        "--chr-length-mb", "0.2",
        "--demo-model", "none",   # skip stdpopsim model lookup
        "--rsid-density", "0",    # skip the dbSNP overlay (no cache)
        "--clinvar-inject-density", "0",
        "--svs-per-person", "0",  # SV emission is per-person; not
                                  # exercised in cohort-only mode
        "--error-rate", "0",
        "--dropout-rate", "0",
        "--workers", "1",
        "--output-dir", str(out_dir),
        "--cache-dir", str(out_dir / "cache"),
        "--mode", mode,
    ]


@unittest.skipUnless(
    _HAVE_BCFTOOLS and _HAVE_TABIX and _HAVE_BGZIP
    and _HAVE_MSPRIME and _HAVE_STDPOPSIM,
    "needs bcftools/tabix/bgzip + msprime + stdpopsim",
)
class CliModePerPersonTest(unittest.TestCase):
    """Default behaviour. Should land exactly the artefacts a pre-Phase-5
    user is used to and nothing new."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(prefix="cli_mode_pp_"))
        rc = cli_module.main(_common_args(cls.tmpdir, "per-person"))
        if rc != 0:
            raise RuntimeError(f"cli.main exited {rc}")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_per_person_vcfs_written(self):
        for i in range(1, 4):
            vcf = self.tmpdir / f"person_{i:04d}.vcf.gz"
            self.assertTrue(vcf.is_file(), f"missing {vcf}")
            self.assertTrue((vcf.with_suffix(".gz.tbi")).is_file()
                            or Path(str(vcf) + ".tbi").is_file())

    def test_no_cohort_bcf_written(self):
        # --mode per-person must not silently add a cohort BCF — that
        # would change the deliverable shape for every existing user.
        cohort_dir = self.tmpdir / "cohort"
        self.assertFalse(cohort_dir.exists(),
                         f"unexpected cohort dir at {cohort_dir}")

    def test_manifest_marks_shape(self):
        manifest = json.loads(
            (self.tmpdir / "manifest.json").read_text())
        self.assertEqual(manifest.get("shape"), "per-person")
        self.assertNotIn("cohort_bcfs", manifest)
        # `people` list still populated as before.
        self.assertEqual(len(manifest["people"]), 3)
        # Top-level `samples` list now lands in every mode so callers
        # get one code path for "list of sample IDs" lookup.
        self.assertEqual(
            manifest["samples"],
            [p["sample_id"] for p in manifest["people"]],
        )


@unittest.skipUnless(
    _HAVE_BCFTOOLS and _HAVE_TABIX and _HAVE_BGZIP
    and _HAVE_MSPRIME and _HAVE_STDPOPSIM,
    "needs bcftools/tabix/bgzip + msprime + stdpopsim",
)
class CliModeCohortTest(unittest.TestCase):
    """Cohort-only mode skips per-person fan-out and lands the BCF as
    the sole genotype deliverable."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(prefix="cli_mode_co_"))
        rc = cli_module.main(_common_args(cls.tmpdir, "cohort"))
        if rc != 0:
            raise RuntimeError(f"cli.main exited {rc}")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_cohort_bcfs_and_indices_exist(self):
        # Phase 5b1 streams cohort mode chromosome-by-chromosome and
        # lands one BCF per chromosome rather than a single combined
        # file. The fixture only simulates chr22 so we expect exactly
        # one cohort.chr22.bcf + its CSI sidecar.
        bcf = self.tmpdir / "cohort" / "cohort.chr22.bcf"
        self.assertTrue(bcf.is_file(), f"missing {bcf}")
        self.assertTrue(Path(str(bcf) + ".csi").is_file(),
                        f"missing CSI index for {bcf}")
        bcfs = sorted((self.tmpdir / "cohort").glob("cohort.chr*.bcf"))
        self.assertEqual(len(bcfs), 1)

    def test_no_per_person_vcfs(self):
        # The point of cohort mode is to skip per-person fan-out
        # entirely. If any person_NNNN.vcf.gz lands here, the early-
        # return path got bypassed.
        for i in range(1, 4):
            vcf = self.tmpdir / f"person_{i:04d}.vcf.gz"
            self.assertFalse(vcf.exists(),
                             f"unexpected per-person VCF: {vcf}")

    def test_manifest_marks_cohort_mode(self):
        manifest = json.loads(
            (self.tmpdir / "manifest.json").read_text())
        self.assertEqual(manifest.get("shape"), "cohort")
        # Phase 5b1 — streamed cohort mode lands per-chromosome BCFs.
        # The cohort_bcfs list shape (already in place from 5a) carries
        # one entry per chromosome simulated.
        self.assertEqual(
            manifest.get("cohort_bcfs"),
            ["cohort/cohort.chr22.bcf"],
        )
        # No per-person list in cohort mode.
        self.assertNotIn("people", manifest)
        # Sample IDs are recorded so a downstream `bcftools view -s`
        # caller knows what to ask for.
        self.assertEqual(len(manifest["samples"]), 3)

    def test_bcf_has_three_sample_columns(self):
        bcf = self.tmpdir / "cohort" / "cohort.chr22.bcf"
        proc = subprocess.run(
            ["bcftools", "query", "-l", str(bcf)],
            capture_output=True, text=True, check=True,
        )
        samples = [s for s in proc.stdout.splitlines() if s.strip()]
        self.assertEqual(len(samples), 3)


@unittest.skipUnless(
    _HAVE_BCFTOOLS and _HAVE_TABIX and _HAVE_BGZIP
    and _HAVE_MSPRIME and _HAVE_STDPOPSIM,
    "needs bcftools/tabix/bgzip + msprime + stdpopsim",
)
class CliModeBothTest(unittest.TestCase):
    """Both deliverables. Manifest carries cohort_bcfs + people list;
    the BCF and per-person VCFs were generated from the same cohort
    so their sample IDs and record counts line up."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(prefix="cli_mode_bo_"))
        rc = cli_module.main(_common_args(cls.tmpdir, "both"))
        if rc != 0:
            raise RuntimeError(f"cli.main exited {rc}")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_both_deliverables_present(self):
        self.assertTrue(
            (self.tmpdir / "cohort" / "cohort.bcf").is_file())
        for i in range(1, 4):
            self.assertTrue(
                (self.tmpdir / f"person_{i:04d}.vcf.gz").is_file())

    def test_manifest_carries_both(self):
        manifest = json.loads(
            (self.tmpdir / "manifest.json").read_text())
        self.assertEqual(manifest.get("shape"), "both")
        self.assertEqual(
            manifest.get("cohort_bcfs"), ["cohort/cohort.bcf"],
        )
        self.assertEqual(len(manifest["people"]), 3)
        # Top-level samples mirrors the per-person list's IDs in order.
        self.assertEqual(
            manifest["samples"],
            [p["sample_id"] for p in manifest["people"]],
        )

    def test_bcf_samples_match_per_person_files(self):
        bcf = self.tmpdir / "cohort" / "cohort.bcf"
        bcf_samples = subprocess.run(
            ["bcftools", "query", "-l", str(bcf)],
            capture_output=True, text=True, check=True,
        ).stdout.splitlines()
        manifest = json.loads(
            (self.tmpdir / "manifest.json").read_text())
        per_person_samples = [p["sample_id"] for p in manifest["people"]]
        # Same IDs in the same order — they were drawn from the same
        # rng pass against the same seed.
        self.assertEqual(bcf_samples, per_person_samples)


if __name__ == "__main__":
    unittest.main()
