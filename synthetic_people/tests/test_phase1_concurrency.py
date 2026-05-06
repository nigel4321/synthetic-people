"""Tests for the Phase 1 concurrency / streaming writer changes.

Three things are verified here:

* ``resolve_workers`` honours the documented contract (auto / explicit /
  non-Linux fallback / negative input).
* ``simulate_cohort`` is deterministic for a given master ``--seed``
  irrespective of ``workers`` — i.e. running with ``workers=1`` and
  ``workers=2`` produces byte-identical site sequences when seeded the
  same way.
* ``write_person_vcf`` produces a `bgzip`-compressed VCF that
  round-trips through ``bcftools view`` after being streamed straight
  through ``bgzip -c`` (no plain ``.vcf`` intermediate).
"""

from __future__ import annotations

import random
import shutil
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.cli import resolve_workers

try:
    import msprime  # noqa: F401
    import stdpopsim  # noqa: F401
    _HAVE_SIM_DEPS = True
except ImportError:
    _HAVE_SIM_DEPS = False


_HAVE_BGZIP = shutil.which("bgzip") is not None
_HAVE_TABIX = shutil.which("tabix") is not None
_HAVE_BCFTOOLS = shutil.which("bcftools") is not None


class TestResolveWorkers(unittest.TestCase):
    def test_negative_rejected(self):
        with self.assertRaises(ValueError):
            resolve_workers(-1)

    def test_explicit_one_returns_one(self):
        # On Linux this is unchanged; on other platforms the floor is
        # also 1, so this assertion holds everywhere.
        self.assertEqual(resolve_workers(1), 1)

    def test_non_linux_caps_to_one(self):
        # Skip on Linux (would just be redundant with test_explicit*).
        if sys.platform == "linux":
            self.skipTest("Linux: parallelism enabled, no cap to test")
        self.assertEqual(resolve_workers(0), 1)
        self.assertEqual(resolve_workers(8), 1)

    def test_explicit_value_passes_through_on_linux(self):
        if sys.platform != "linux":
            self.skipTest("non-Linux always caps to 1")
        self.assertEqual(resolve_workers(3), 3)


@unittest.skipUnless(_HAVE_SIM_DEPS, "msprime/stdpopsim not installed")
@unittest.skipUnless(sys.platform == "linux",
                     "fork-based parallelism is Linux-only")
class TestSimulateCohortDeterminism(unittest.TestCase):
    """workers=1 must match workers=N for the same master --seed."""

    SEED = 12345
    CHROMS = ["21", "22"]
    LENGTH_MB = 0.15
    N_PEOPLE = 4

    def _run(self, workers: int) -> list:
        from syntheticgen.coalescent import simulate_cohort
        rng = random.Random(self.SEED)
        return simulate_cohort(
            chromosomes=self.CHROMS, build="GRCh38",
            n_people=self.N_PEOPLE, length_mb=self.LENGTH_MB,
            demo_model=None, population="CEU",
            rec_rate=1e-8, mu=1.29e-8, rng=rng,
            verbose=False, workers=workers,
        )

    def test_serial_reproducible(self):
        a = self._run(workers=1)
        b = self._run(workers=1)
        self.assertEqual(_summary(a), _summary(b))

    def test_parallel_matches_serial(self):
        serial = self._run(workers=1)
        parallel = self._run(workers=2)
        self.assertEqual(_summary(serial), _summary(parallel))


def _summary(sites: list) -> list:
    """Project a sites list to its serialisable, comparable contents.

    Phase 5c: cohort sites carry sparse carriers rather than dense
    GT lists; sort the carrier tuples for stable comparison.
    """
    return [
        (s["chrom"], s["pos"], s["ref"], tuple(s["alts"]),
         tuple(s["acs"]), tuple(sorted(s["carriers"])))
        for s in sites
    ]


@unittest.skipUnless(_HAVE_BGZIP and _HAVE_TABIX,
                     "bgzip/tabix not on PATH")
class TestWriterBgzipPipe(unittest.TestCase):
    """The writer streams through `bgzip -c`; round-trip via gunzip / tabix."""

    def setUp(self):
        self.tmp = Path("/tmp/test_phase1_writer")
        if self.tmp.exists():
            shutil.rmtree(self.tmp)
        self.tmp.mkdir(parents=True)

    def tearDown(self):
        if self.tmp.exists():
            shutil.rmtree(self.tmp)

    def _person(self) -> dict:
        return {
            "sample_id": "PHASE1_SMOKE",
            "highlighted": {
                "id": "rs1", "chrom": "22", "pos": 1000,
                "ref": "A", "alts": ["G"], "gt": "0|1",
            },
            "background": [
                {"id": ".", "chrom": "22", "pos": 2000, "ref": "C",
                 "alts": ["T"], "gt": "0|1"},
                {"id": ".", "chrom": "22", "pos": 3000, "ref": "G",
                 "alts": ["A"], "gt": "1|1"},
            ],
        }

    def test_roundtrip_via_bgzip_decompression(self):
        from syntheticgen.writer import write_person_vcf
        out = self.tmp / "person_0001.vcf.gz"
        write_person_vcf(out, self._person(), "GRCh38",
                         random.Random(0))
        self.assertTrue(out.exists())
        # Tabix index landed too.
        self.assertTrue((self.tmp / "person_0001.vcf.gz.tbi").exists())
        # bgzip header magic — first two bytes are gzip magic 1f 8b.
        with open(out, "rb") as f:
            self.assertEqual(f.read(2), b"\x1f\x8b")
        # Decompress via bgzip -d -c and check the variant rows are
        # present.
        decoded = subprocess.run(
            ["bgzip", "-d", "-c", str(out)],
            check=True, capture_output=True,
        ).stdout.decode("utf-8")
        # Three variants expected.
        body_lines = [ln for ln in decoded.splitlines()
                      if ln and not ln.startswith("#")]
        self.assertEqual(len(body_lines), 3)
        self.assertIn("HIGHLIGHT", body_lines[0])

    @unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
    def test_roundtrip_via_bcftools(self):
        from syntheticgen.writer import write_person_vcf
        out = self.tmp / "person_0001.vcf.gz"
        write_person_vcf(out, self._person(), "GRCh38",
                         random.Random(0))
        result = subprocess.run(
            ["bcftools", "view", "-H", str(out)],
            check=True, capture_output=True,
        )
        rows = [ln for ln in result.stdout.decode("utf-8").splitlines()
                if ln.strip()]
        self.assertEqual(len(rows), 3)
        # Sample column carries the GT we wrote.
        self.assertIn("0|1", rows[0])
        self.assertIn("1|1", rows[2])

    def test_no_plain_vcf_left_behind(self):
        """The pipe-into-bgzip path must not create a sibling .vcf."""
        from syntheticgen.writer import write_person_vcf
        out = self.tmp / "person_0001.vcf.gz"
        write_person_vcf(out, self._person(), "GRCh38",
                         random.Random(0))
        plain = self.tmp / "person_0001.vcf"
        self.assertFalse(plain.exists(),
                         f"Stale {plain} found — writer should not "
                         f"leave a plain VCF behind")


if __name__ == "__main__":
    unittest.main(verbosity=2)
