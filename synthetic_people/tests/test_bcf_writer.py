"""Tests for syntheticgen/bcf_writer.py — multi-sample cohort BCF
streaming.

The Phase 5 plan moves cohort genotype state onto disk as BCF; this
module is the writer that lands them. Per-person VCFs are then derived
via ``bcftools view -s SAMPLE`` from the cohort BCF in the rest of the
pipeline.

Tests gate on bcftools / bgzip on PATH (matching the existing
``test_phase1_concurrency.py`` pattern) so CI hosts that already
install htslib for the rest of the suite exercise this path too. Hosts
without htslib skip cleanly.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.bcf_writer import (
    CohortBcfWriter,
    build_cohort_header,
)


_HAVE_BCFTOOLS = shutil.which("bcftools") is not None


def _site(pos: int, ref: str, alt: str, gts: list,
          chrom: str = "22", site_id: str = ".",
          **extra) -> dict:
    """Build a cohort-site dict matching simulate_cohort's output shape."""
    nalt = sum(int(t) for gt in gts for t in gt.split("|"))
    site = {
        "chrom": chrom,
        "pos": pos,
        "id": site_id,
        "ref": ref,
        "alts": [alt],
        "acs": [nalt],
        "gts": list(gts),
    }
    site.update(extra)
    return site


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
class CohortBcfWriterRoundTripTest(unittest.TestCase):
    """A site list written through the cohort writer should be readable
    by ``bcftools view`` with the same fields back."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="cohort_bcf_rt_")
        cls.samples = ["S1", "S2", "S3"]
        cls.sites = [
            _site(1000, "A", "G", ["0|0", "0|1", "1|1"], site_id="rs1"),
            _site(2000, "C", "T", ["0|0", "0|1", "0|0"]),
            _site(
                3000, "G", "<DEL>", ["0|0", "0|0", "0|1"], site_id="rs2",
                svtype="DEL", svlen=-500, end=3500, cipos=(-50, 50),
            ),
        ]
        cls.bcf = Path(cls.tmpdir) / "cohort.bcf"
        with CohortBcfWriter(cls.bcf, "GRCh38", cls.samples) as w:
            n = w.write_sites(cls.sites)
        cls.n_written = n

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_n_written_matches_input(self):
        self.assertEqual(self.n_written, len(self.sites))

    def test_bcf_and_index_exist(self):
        self.assertTrue(self.bcf.is_file())
        # CSI is the htslib-native index for BCF (TBI is text-only).
        self.assertTrue(Path(str(self.bcf) + ".csi").is_file())

    def test_bcftools_view_can_parse(self):
        proc = subprocess.run(
            ["bcftools", "view", "-h", str(self.bcf)],
            capture_output=True, text=True, check=True,
        )
        # All declared sample names should land in the #CHROM line.
        chrom_line = next(
            l for l in proc.stdout.splitlines() if l.startswith("#CHROM"))
        for s in self.samples:
            self.assertIn(s, chrom_line)

    def test_record_round_trip(self):
        proc = subprocess.run(
            ["bcftools", "view", "-H", str(self.bcf)],
            capture_output=True, text=True, check=True,
        )
        rows = [l.split("\t") for l in proc.stdout.splitlines() if l.strip()]
        self.assertEqual(len(rows), len(self.sites))
        # Spot-check the SV record — it carries the most metadata.
        sv_row = next(r for r in rows if r[4] == "<DEL>")
        info_kvs = dict(
            kv.split("=", 1) if "=" in kv else (kv, "")
            for kv in sv_row[7].split(";")
        )
        self.assertEqual(info_kvs["SVTYPE"], "DEL")
        self.assertEqual(info_kvs["SVLEN"], "-500")
        self.assertEqual(info_kvs["END"], "3500")
        self.assertEqual(info_kvs["CIPOS"], "-50,50")
        # Per-sample GT block: 3 cells after FORMAT.
        self.assertEqual(sv_row[9:], ["0|0", "0|0", "0|1"])

    def test_cohort_level_ac_an_af(self):
        proc = subprocess.run(
            ["bcftools", "view", "-H", str(self.bcf)],
            capture_output=True, text=True, check=True,
        )
        # First site has 1×0|0, 1×0|1, 1×1|1 → AC=3, AN=6, AF=0.5.
        first = proc.stdout.splitlines()[0].split("\t")
        info = dict(
            kv.split("=", 1) for kv in first[7].split(";") if "=" in kv)
        self.assertEqual(info["AC"], "3")
        self.assertEqual(info["AN"], "6")
        # bcftools strips trailing zeros on output — assert the float
        # value rather than the string form.
        self.assertAlmostEqual(float(info["AF"]), 0.5, places=5)


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
class CohortBcfPerSampleExtractTest(unittest.TestCase):
    """``bcftools view -s SAMPLE`` is the per-person derivation step the
    rest of Phase 5 uses; it must produce single-sample VCFs whose
    record set matches what the cohort writer originally received."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="cohort_bcf_extract_")
        samples = ["alice", "bob", "charlie"]
        cls.sites = [
            _site(1000, "A", "G", ["0|1", "0|0", "1|1"]),
            _site(2000, "C", "T", ["0|0", "0|1", "0|0"]),
        ]
        cls.bcf = Path(cls.tmpdir) / "cohort.bcf"
        with CohortBcfWriter(cls.bcf, "GRCh38", samples) as w:
            w.write_sites(cls.sites)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_per_sample_view_isolates_one_column(self):
        proc = subprocess.run(
            ["bcftools", "view", "-s", "bob", str(self.bcf)],
            capture_output=True, text=True, check=True,
        )
        chrom_line = next(
            l for l in proc.stdout.splitlines() if l.startswith("#CHROM"))
        self.assertTrue(chrom_line.endswith("\tbob"))
        rows = [l.split("\t") for l in proc.stdout.splitlines()
                if l.strip() and not l.startswith("#")]
        # Bob's GTs from the input fixtures.
        self.assertEqual([r[9] for r in rows], ["0|0", "0|1"])


class CohortBcfWriterArgValidationTest(unittest.TestCase):
    """Catches misshape mistakes early so a streaming caller doesn't
    silently emit a corrupt BCF for one chromosome and fail downstream."""

    @unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
    def test_wrong_gt_count_raises(self):
        tmpdir = tempfile.mkdtemp(prefix="cohort_bcf_bad_")
        try:
            bcf = Path(tmpdir) / "bad.bcf"
            samples = ["S1", "S2", "S3"]
            site_with_two_gts = _site(1000, "A", "G", ["0|1", "1|1"])
            with self.assertRaisesRegex(ValueError, "expected 3"):
                with CohortBcfWriter(bcf, "GRCh38", samples) as w:
                    w.write_site(site_with_two_gts)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_header_lists_all_samples(self):
        # Doesn't need bcftools — header is just text assembly.
        hdr = build_cohort_header("GRCh38", ["alpha", "beta", "gamma"])
        chrom_line = [l for l in hdr.splitlines()
                      if l.startswith("#CHROM")][0]
        self.assertIn("alpha", chrom_line)
        self.assertIn("beta", chrom_line)
        self.assertIn("gamma", chrom_line)
        # Tab-separated, in order, after the eight fixed columns.
        fields = chrom_line.split("\t")
        self.assertEqual(fields[-3:], ["alpha", "beta", "gamma"])

    def test_header_declares_format_and_info(self):
        hdr = build_cohort_header("GRCh38", ["S1"])
        # Sanity: every tag the writer / per-person derivation will use
        # is declared. If any of these gets dropped the bcftools-stats
        # path in the variant-scan pipeline would silently underreport.
        for tag in ("##FORMAT=<ID=GT", "##INFO=<ID=AC", "##INFO=<ID=AN",
                    "##INFO=<ID=AF", "##INFO=<ID=SVTYPE",
                    "##INFO=<ID=CLNSIG", "##ALT=<ID=DEL"):
            self.assertIn(tag, hdr)


if __name__ == "__main__":
    unittest.main()
