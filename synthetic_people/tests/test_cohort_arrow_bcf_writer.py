"""Parity tests for the Arrow-source parallel BCF write (Phase 5d.1).

The load-bearing claim of Phase 5d.1 is that
``write_cohort_bcf_parallel_from_arrow(arrow_path, ...)`` produces
output indistinguishable from ``write_cohort_bcf_parallel(...,
sites)`` given the same cohort. These tests pin that claim against
fixtures of varying shape (biallelic, multi-allelic, with/without
overlay INFO fields) and at multiple worker counts.

Skipped cleanly if pyarrow or bcftools is unavailable.
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

try:
    import pyarrow  # noqa: F401
    HAS_PYARROW = True
except ImportError:
    HAS_PYARROW = False

_HAVE_BCFTOOLS = shutil.which("bcftools") is not None

from syntheticgen.bcf_writer import (  # noqa: E402
    write_cohort_bcf_parallel,
)
from syntheticgen.cohort_sites import carriers_from_dense_gts  # noqa: E402

if HAS_PYARROW:
    from syntheticgen.bcf_writer import (  # noqa: E402
        write_cohort_bcf_parallel_from_arrow,
    )
    from syntheticgen.cohort_arrow import write_arrow_file  # noqa: E402


def _site(pos: int, gts: list, ref: str = "A", alt: str = "C",
          chrom: str = "22", site_id: str = ".", **extra) -> dict:
    """Build a cohort-site dict in the carriers shape — mirrors the
    helper in test_cohort_parallel_write.py."""
    nalt = sum(int(t) for gt in gts for t in gt.split("|"))
    site = {
        "chrom": chrom,
        "pos": pos,
        "id": site_id,
        "ref": ref,
        "alts": [alt],
        "afs": [nalt / (2 * len(gts)) if gts else 0],
        "acs": [nalt],
        "n_haplotypes": 2 * len(gts),
        "carriers": carriers_from_dense_gts(list(gts)),
    }
    site.update(extra)
    return site


def _bcf_data_md5(bcf_path: Path) -> str:
    """Content hash of a BCF's data section — same helper as
    test_cohort_parallel_write.py. Pulls fields into a fixed-order
    line format via ``bcftools query`` so the hash isn't sensitive
    to INFO-field ordering."""
    out = subprocess.check_output(
        ["bcftools", "query",
         "-f", "%CHROM\t%POS\t%REF\t%ALT\t"
               "%INFO/AC\t%INFO/AN\t%INFO/AF[\t%GT]\n",
         str(bcf_path)])
    return hashlib.md5(out).hexdigest()


def _bcf_sample_columns(bcf_path: Path) -> list[str]:
    out = subprocess.check_output(
        ["bcftools", "query", "-l", str(bcf_path)],
        text=True,
    )
    return [s for s in out.strip().splitlines() if s]


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
@unittest.skipUnless(HAS_PYARROW, "pyarrow not installed")
class ArrowBcfWriteParityTest(unittest.TestCase):
    """The core parity claim: same cohort -> same merged BCF, by
    content, regardless of source (sites-list vs Arrow IPC)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.sample_ids = [f"S{i:02d}" for i in range(8)]
        self.sites = [
            _site(100, ["0|0", "0|1", "1|1", "0|0",
                        "1|0", "0|0", "0|1", "1|1"]),
            _site(200, ["1|1", "0|0", "0|1", "0|1",
                        "1|0", "1|1", "0|0", "0|1"]),
            _site(300, ["0|0", "0|0", "0|0", "1|0",
                        "0|1", "1|1", "0|0", "0|0"]),
            _site(400, ["1|1", "1|1", "0|0", "0|0",
                        "0|1", "0|1", "1|0", "0|0"]),
            _site(500, ["0|0", "1|0", "0|1", "1|1",
                        "0|0", "0|0", "1|1", "0|1"]),
        ]

    def tearDown(self):
        self.tmp.cleanup()

    def _write_via_sites_list(self, workers: int) -> Path:
        out = self.dir / f"sites_list_w{workers}.bcf"
        write_cohort_bcf_parallel(
            out, "GRCh37", self.sample_ids, list(self.sites), workers,
        )
        return out

    def _write_via_arrow(self, workers: int) -> Path:
        arrow = self.dir / f"cohort_w{workers}.arrow"
        write_arrow_file(arrow, "22", len(self.sample_ids), iter(self.sites))
        out = self.dir / f"arrow_w{workers}.bcf"
        write_cohort_bcf_parallel_from_arrow(
            arrow, out, "GRCh37", self.sample_ids, workers,
        )
        return out

    def test_parity_workers_1(self):
        a = self._write_via_sites_list(1)
        b = self._write_via_arrow(1)
        self.assertEqual(_bcf_data_md5(a), _bcf_data_md5(b))
        self.assertEqual(_bcf_sample_columns(a), _bcf_sample_columns(b))

    def test_parity_workers_2(self):
        a = self._write_via_sites_list(2)
        b = self._write_via_arrow(2)
        self.assertEqual(_bcf_data_md5(a), _bcf_data_md5(b))
        self.assertEqual(_bcf_sample_columns(a), _bcf_sample_columns(b))

    def test_parity_workers_4(self):
        a = self._write_via_sites_list(4)
        b = self._write_via_arrow(4)
        self.assertEqual(_bcf_data_md5(a), _bcf_data_md5(b))
        self.assertEqual(_bcf_sample_columns(a), _bcf_sample_columns(b))

    def test_parity_workers_8(self):
        # workers == n_samples — every slice has one person.
        a = self._write_via_sites_list(8)
        b = self._write_via_arrow(8)
        self.assertEqual(_bcf_data_md5(a), _bcf_data_md5(b))
        self.assertEqual(_bcf_sample_columns(a), _bcf_sample_columns(b))

    def test_parity_workers_exceeds_n_samples(self):
        # workers > n_samples — both paths fall back to serial.
        a = self._write_via_sites_list(16)
        b = self._write_via_arrow(16)
        self.assertEqual(_bcf_data_md5(a), _bcf_data_md5(b))
        self.assertEqual(_bcf_sample_columns(a), _bcf_sample_columns(b))


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
@unittest.skipUnless(HAS_PYARROW, "pyarrow not installed")
class ArrowBcfMultiAllelicParityTest(unittest.TestCase):
    """Multi-allelic sites round-trip through Arrow + the BCF write
    correctly. Sparse carriers preserve allele indices > 1."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.sample_ids = [f"S{i:02d}" for i in range(6)]
        # Multi-allelic site: ref=A, alts=[C, G]; carriers cover both.
        self.sites = [{
            "chrom": "22",
            "pos": 100,
            "id": ".",
            "ref": "A",
            "alts": ["C", "G"],
            "acs": [3, 2],
            "afs": [3 / 12, 2 / 12],
            "n_haplotypes": 12,
            "carriers": [(0, 1), (2, 1), (4, 1), (5, 2), (7, 2)],
        }]

    def tearDown(self):
        self.tmp.cleanup()

    def test_multi_allelic_parity_workers_3(self):
        out_a = self.dir / "ml_sites.bcf"
        out_b = self.dir / "ml_arrow.bcf"
        arrow = self.dir / "ml.arrow"
        write_cohort_bcf_parallel(
            out_a, "GRCh37", self.sample_ids, list(self.sites), 3,
        )
        write_arrow_file(arrow, "22", len(self.sample_ids), iter(self.sites))
        write_cohort_bcf_parallel_from_arrow(
            arrow, out_b, "GRCh37", self.sample_ids, 3,
        )
        self.assertEqual(_bcf_data_md5(out_a), _bcf_data_md5(out_b))


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
@unittest.skipUnless(HAS_PYARROW, "pyarrow not installed")
class ArrowBcfOverlayParityTest(unittest.TestCase):
    """Overlay INFO fields (CLNSIG, CLNDN, COSMIC_ID, COSMIC_GENE,
    SVTYPE, SVLEN, END, CIPOS) survive the Arrow round-trip into
    the final BCF."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.sample_ids = [f"S{i:02d}" for i in range(4)]
        self.sites = [
            _site(100, ["0|0", "0|1", "1|1", "1|0"],
                  site_id="rs100",
                  clnsig="Pathogenic",
                  clndn="Some_Disease"),
            _site(200, ["1|1", "1|1", "0|0", "0|0"],
                  cosmic_id="COSV99",
                  cosmic_gene="TP53"),
            _site(300, ["0|0", "0|0", "1|0", "0|0"],
                  svtype="DEL", svlen=-200, end=500,
                  cipos=(-50, 50)),
            _site(400, ["0|1", "0|1", "0|1", "0|1"]),  # plain site, no overlays
        ]

    def tearDown(self):
        self.tmp.cleanup()

    def test_overlay_fields_present_after_arrow_round_trip(self):
        out_a = self.dir / "ov_sites.bcf"
        out_b = self.dir / "ov_arrow.bcf"
        arrow = self.dir / "ov.arrow"
        write_cohort_bcf_parallel(
            out_a, "GRCh37", self.sample_ids, list(self.sites), 2,
        )
        write_arrow_file(arrow, "22", len(self.sample_ids), iter(self.sites))
        write_cohort_bcf_parallel_from_arrow(
            arrow, out_b, "GRCh37", self.sample_ids, 2,
        )

        # Pull INFO fields explicitly to verify overlay survival.
        info_query = (
            "%CHROM\t%POS\t%INFO/CLNSIG\t%INFO/CLNDN\t"
            "%INFO/COSMIC_ID\t%INFO/COSMIC_GENE\t"
            "%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/END\t%INFO/CIPOS\n"
        )
        out_a_info = subprocess.check_output(
            ["bcftools", "query", "-f", info_query, str(out_a)],
            text=True,
        )
        out_b_info = subprocess.check_output(
            ["bcftools", "query", "-f", info_query, str(out_b)],
            text=True,
        )
        self.assertEqual(out_a_info, out_b_info)


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
@unittest.skipUnless(HAS_PYARROW, "pyarrow not installed")
class ArrowBcfWorkerErrorPropagationTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_missing_arrow_file_raises(self):
        # Missing Arrow file: workers fail when read_arrow_slice tries
        # to mmap. The writer must surface the failure rather than
        # silently produce a missing-partial.
        sample_ids = [f"S{i:02d}" for i in range(4)]
        out = self.dir / "out.bcf"
        bogus = self.dir / "does_not_exist.arrow"
        with self.assertRaises(Exception):
            write_cohort_bcf_parallel_from_arrow(
                bogus, out, "GRCh37", sample_ids, 2,
            )


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
@unittest.skipUnless(HAS_PYARROW, "pyarrow not installed")
class ArrowBcfSerialFallbackTest(unittest.TestCase):
    """workers <= 1 should use the in-process write path; output must
    match the parallel path on the same input (the parallel path's
    own determinism is covered by its existing tests)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.sample_ids = [f"S{i:02d}" for i in range(4)]
        self.sites = [
            _site(100, ["0|0", "0|1", "1|0", "1|1"]),
            _site(200, ["1|1", "0|0", "0|1", "0|0"]),
        ]

    def tearDown(self):
        self.tmp.cleanup()

    def test_workers_zero_falls_back_to_serial(self):
        # workers=0 should still produce a valid BCF (treated as
        # serial), matching the existing path's behaviour.
        arrow = self.dir / "cohort.arrow"
        out_par = self.dir / "par.bcf"
        out_ser = self.dir / "ser.bcf"
        write_arrow_file(
            arrow, "22", len(self.sample_ids), iter(self.sites)
        )
        write_cohort_bcf_parallel_from_arrow(
            arrow, out_par, "GRCh37", self.sample_ids, 4,
        )
        write_cohort_bcf_parallel_from_arrow(
            arrow, out_ser, "GRCh37", self.sample_ids, 0,
        )
        self.assertEqual(_bcf_data_md5(out_par), _bcf_data_md5(out_ser))


if __name__ == "__main__":
    unittest.main()
