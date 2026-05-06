"""Phase 5b1 — streaming cohort flow tests.

Phase 5b1 adds a chromosome-by-chromosome streaming variant of the
``--mode cohort`` path: each chromosome's sites are simulated,
overlaid, written to its own ``cohort.chr<N>.bcf``, and freed before
the next chromosome is simulated. This file covers the contract:

- The streaming generator ``simulate_cohort_iter`` yields one
  chunk per requested chromosome, with the right shape.
- ``--mode cohort`` over multiple chromosomes lands one BCF per
  chromosome, all carrying the same sample columns.
- Determinism: same seed → byte-identical streamed BCFs across
  re-runs (this is the determinism contract Phase 5b1 actually
  guarantees; cross-path identity vs the in-memory cohort flow
  is *not* preserved because rng consumption order differs).

Heavy paths gate on bcftools/tabix/bgzip + msprime + stdpopsim,
matching ``test_cli_modes.py``.
"""

from __future__ import annotations

import importlib.util
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen import cli as cli_module
from syntheticgen.coalescent import simulate_cohort_iter


_HAVE_BCFTOOLS = shutil.which("bcftools") is not None
_HAVE_TABIX = shutil.which("tabix") is not None
_HAVE_BGZIP = shutil.which("bgzip") is not None
_HAVE_MSPRIME = importlib.util.find_spec("msprime") is not None
_HAVE_STDPOPSIM = importlib.util.find_spec("stdpopsim") is not None


@unittest.skipUnless(_HAVE_MSPRIME, "msprime not available")
class SimulateCohortIterTest(unittest.TestCase):
    """The streaming generator must yield exactly one chunk per
    chromosome, in the requested order, with each chunk shaped like
    the existing flat-list output."""

    def test_yields_one_chunk_per_chromosome(self):
        rng = random.Random(42)
        chunks = list(simulate_cohort_iter(
            chromosomes=["22", "21"],
            build="GRCh38", n_people=2, length_mb=0.05,
            demo_model=None, population="CEU",
            rec_rate=1e-8, mu=1.29e-8, rng=rng, workers=1,
        ))
        self.assertEqual([c[0] for c in chunks], ["22", "21"])
        for chrom, sites in chunks:
            self.assertGreater(len(sites), 0)
            for s in sites:
                # Same site dict shape the flat-list path produces.
                for key in ("chrom", "pos", "ref", "alts", "gts"):
                    self.assertIn(key, s)
                self.assertEqual(s["chrom"], chrom)


def _common_args(out_dir: Path, n: int = 3, chroms: str = "22") -> list:
    return [
        "--n", str(n),
        "--seed", "42",
        "--build", "GRCh38",
        "--chromosomes", chroms,
        "--chr-length-mb", "0.2",
        "--demo-model", "none",
        "--rsid-density", "0",
        "--clinvar-inject-density", "0",
        "--svs-per-person", "0",
        "--error-rate", "0",
        "--dropout-rate", "0",
        "--workers", "1",
        "--output-dir", str(out_dir),
        "--cache-dir", str(out_dir / "cache"),
        "--mode", "cohort",
    ]


@unittest.skipUnless(
    _HAVE_BCFTOOLS and _HAVE_TABIX and _HAVE_BGZIP
    and _HAVE_MSPRIME and _HAVE_STDPOPSIM,
    "needs bcftools/tabix/bgzip + msprime + stdpopsim",
)
class CohortStreamedMultiChromTest(unittest.TestCase):
    """Multi-chromosome streamed cohort runs land one BCF per chrom,
    all carrying the same sample columns. The manifest's cohort_bcfs
    list reflects every chromosome simulated."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(prefix="cohort_stream_multi_"))
        rc = cli_module.main(_common_args(cls.tmpdir, chroms="20,21,22"))
        if rc != 0:
            raise RuntimeError(f"cli.main exited {rc}")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_one_bcf_per_chromosome(self):
        bcfs = sorted((self.tmpdir / "cohort").glob("cohort.chr*.bcf"))
        names = [p.name for p in bcfs]
        self.assertEqual(
            names,
            ["cohort.chr20.bcf", "cohort.chr21.bcf", "cohort.chr22.bcf"],
        )
        for bcf in bcfs:
            self.assertTrue(Path(str(bcf) + ".csi").is_file(),
                            f"missing CSI for {bcf}")

    def test_manifest_lists_every_chromosome_bcf(self):
        manifest = json.loads(
            (self.tmpdir / "manifest.json").read_text())
        self.assertEqual(
            manifest["cohort_bcfs"],
            [
                "cohort/cohort.chr20.bcf",
                "cohort/cohort.chr21.bcf",
                "cohort/cohort.chr22.bcf",
            ],
        )

    def test_all_bcfs_share_sample_columns(self):
        manifest = json.loads(
            (self.tmpdir / "manifest.json").read_text())
        expected = manifest["samples"]
        for path in manifest["cohort_bcfs"]:
            proc = subprocess.run(
                ["bcftools", "query", "-l", str(self.tmpdir / path)],
                capture_output=True, text=True, check=True,
            )
            actual = [s for s in proc.stdout.splitlines() if s.strip()]
            self.assertEqual(actual, expected,
                             f"sample mismatch in {path}")

    def test_each_bcf_holds_only_its_chromosome(self):
        for path in (self.tmpdir / "cohort").glob("cohort.chr*.bcf"):
            chrom_in_name = path.name.replace("cohort.chr", "").replace(
                ".bcf", "")
            proc = subprocess.run(
                ["bcftools", "view", "-H", str(path)],
                capture_output=True, text=True, check=True,
            )
            seen_chroms = set()
            for line in proc.stdout.splitlines():
                if not line.strip():
                    continue
                seen_chroms.add(line.split("\t", 1)[0])
            self.assertEqual(seen_chroms, {chrom_in_name},
                             f"{path.name} holds wrong chrom(s): "
                             f"{seen_chroms}")


@unittest.skipUnless(
    _HAVE_BCFTOOLS and _HAVE_TABIX and _HAVE_BGZIP
    and _HAVE_MSPRIME and _HAVE_STDPOPSIM,
    "needs bcftools/tabix/bgzip + msprime + stdpopsim",
)
class CohortStreamedDeterminismTest(unittest.TestCase):
    """Same seed → identical streamed BCFs across two re-runs.

    This is the determinism contract 5b1 actually guarantees. The
    streamed path's rng consumption order differs from the 5a in-
    memory path (overlays apply per-chunk vs globally), so cross-path
    byte-identity is *not* preserved — that's a separate concern
    documented in the plan.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir_a = Path(tempfile.mkdtemp(prefix="cohort_det_a_"))
        cls.tmpdir_b = Path(tempfile.mkdtemp(prefix="cohort_det_b_"))
        for td in (cls.tmpdir_a, cls.tmpdir_b):
            rc = cli_module.main(_common_args(td, chroms="22"))
            if rc != 0:
                raise RuntimeError(f"cli.main exited {rc} for {td}")

    @classmethod
    def tearDownClass(cls):
        for td in (cls.tmpdir_a, cls.tmpdir_b):
            shutil.rmtree(td, ignore_errors=True)

    def test_per_record_equivalence_across_runs(self):
        # bcftools view normalises BCF→VCF formatting in a few small
        # ways (timestamp in the BCF header, occasionally trailing
        # whitespace) so we compare the de-headered record streams
        # rather than raw bytes.
        def _records(td: Path) -> list:
            bcf = td / "cohort" / "cohort.chr22.bcf"
            proc = subprocess.run(
                ["bcftools", "view", "-H", str(bcf)],
                capture_output=True, text=True, check=True,
            )
            return proc.stdout.splitlines()
        rec_a = _records(self.tmpdir_a)
        rec_b = _records(self.tmpdir_b)
        self.assertEqual(rec_a, rec_b,
                         f"streamed cohort records diverged across runs "
                         f"at the same seed: "
                         f"{len(rec_a)} vs {len(rec_b)} records")


if __name__ == "__main__":
    unittest.main()
