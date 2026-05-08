"""Tests for Phase 5e Phase A — parallel sample-slice cohort BCF write.

The structural change replaces the per-chrom serial
``CohortBcfWriter`` block with ``write_cohort_bcf_parallel``, which
forks W workers (each writing a contiguous sample-slice partial BCF
for the same site set) and then runs ``bcftools merge`` to combine
them. The contract these tests pin:

- **Output equivalence:** at any ``workers``, the final cohort BCF
  has the same sites with the same per-sample GT columns as the
  serial-write reference. ``bcftools merge`` is asked to combine
  sample-disjoint partials, so the join collapses to a sample-column
  concatenation in the requested order.
- **Determinism across worker counts:** ``workers ∈ {1, 2, 4, 8}``
  on the same input list produces byte-identical decoded output
  (modulo ``bcftools view -H`` for the data section — bgzip block
  boundaries can drift).
- **Sample-slice helper parity:** ``dense_gts_from_carriers_slice``
  matches ``dense_gts_from_carriers`` on the corresponding range
  for any slice bounds.

bcftools is required for these tests (write + view + index + merge);
they skip cleanly on hosts without it.
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

from syntheticgen.bcf_writer import (
    CohortBcfWriter,
    _split_into_slices,
    write_cohort_bcf_parallel,
)
from syntheticgen.cohort_sites import (
    carriers_from_dense_gts,
    dense_gts_from_carriers,
    dense_gts_from_carriers_slice,
)


_HAVE_BCFTOOLS = shutil.which("bcftools") is not None


def _site(pos: int, ref: str, alt: str, gts: list,
          chrom: str = "22", site_id: str = ".",
          **extra) -> dict:
    """Build a cohort-site dict in the carriers shape that the
    streamed cohort writer expects (mirrors ``test_bcf_writer.py``)."""
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
    """Deterministic content hash of a BCF's data section.

    Pulls the data fields we care about into a fixed-order line
    format via ``bcftools query`` so the hash isn't sensitive to
    INFO-field ordering — ``bcftools merge`` re-sorts INFO fields
    relative to ``bcftools view``, but the values are identical.
    What we want to pin is: per-row chrom/pos/ref/alt + the AC/AN/AF
    values + the per-sample GT block.
    """
    out = subprocess.check_output(
        ["bcftools", "query",
         "-f", "%CHROM\t%POS\t%REF\t%ALT\t"
               "%INFO/AC\t%INFO/AN\t%INFO/AF[\t%GT]\n",
         str(bcf_path)])
    return hashlib.md5(out).hexdigest()


def _bcf_sample_columns(bcf_path: Path) -> list[str]:
    """Sample names in the order the BCF declares them."""
    out = subprocess.check_output(
        ["bcftools", "query", "-l", str(bcf_path)],
        text=True,
    )
    return [s for s in out.strip().splitlines() if s]


class DenseGtsFromCarriersSliceTest(unittest.TestCase):
    """Phase 5e: the slice helper has to agree with the full-expansion
    helper on the slice's range, exactly. This is the unit-level
    correctness gate for the parallel-write path."""

    def setUp(self):
        # Hand-built cohort with 6 people × 1 multi-allelic site so
        # the helper exercises both hom-ref and carrier slots.
        self.gts = ["0|0", "0|1", "1|1", "0|0", "1|0", "1|1"]
        self.carriers = carriers_from_dense_gts(self.gts)
        self.full = dense_gts_from_carriers(self.carriers, len(self.gts))

    def test_full_range_matches_full_helper(self):
        # Slicing [0, n) should produce the same list as the
        # full-expansion helper.
        n = len(self.gts)
        sliced = dense_gts_from_carriers_slice(self.carriers, 0, n)
        self.assertEqual(sliced, self.full)

    def test_prefix_slice(self):
        sliced = dense_gts_from_carriers_slice(self.carriers, 0, 3)
        self.assertEqual(sliced, self.full[0:3])

    def test_suffix_slice(self):
        sliced = dense_gts_from_carriers_slice(self.carriers, 3, 6)
        self.assertEqual(sliced, self.full[3:6])

    def test_middle_slice(self):
        sliced = dense_gts_from_carriers_slice(self.carriers, 2, 5)
        self.assertEqual(sliced, self.full[2:5])

    def test_singleton_slice(self):
        for i in range(len(self.gts)):
            sliced = dense_gts_from_carriers_slice(
                self.carriers, i, i + 1)
            self.assertEqual(sliced, [self.full[i]])

    def test_empty_slice(self):
        # An empty slice (lo == hi) should return an empty list, not
        # raise. Matters because _split_into_slices may produce
        # no entries when workers > n; we want the helper to be
        # defensive.
        self.assertEqual(
            dense_gts_from_carriers_slice(self.carriers, 3, 3), [])


class SplitIntoSlicesTest(unittest.TestCase):
    """Slice-splitting is the determinism source for the parallel
    write — same n + same workers must always produce the same
    boundary list."""

    def test_even_split(self):
        self.assertEqual(_split_into_slices(8, 4),
                         [(0, 2), (2, 4), (4, 6), (6, 8)])

    def test_uneven_split_distributes_remainder_to_first_slices(self):
        # n=10, W=3 → sizes [4, 3, 3] (the +1 lands on the first ``rem``
        # slices, deterministic). We pin the exact shape because the
        # parallel-write merge order depends on it.
        self.assertEqual(_split_into_slices(10, 3),
                         [(0, 4), (4, 7), (7, 10)])

    def test_workers_exceeds_n_drops_empty_slices(self):
        # n=3, W=8 → 3 single-person slices, 5 empty ones dropped.
        self.assertEqual(_split_into_slices(3, 8),
                         [(0, 1), (1, 2), (2, 3)])

    def test_one_worker_returns_full_range(self):
        self.assertEqual(_split_into_slices(5, 1), [(0, 5)])

    def test_zero_n_returns_empty(self):
        self.assertEqual(_split_into_slices(0, 4), [])

    def test_zero_workers_returns_empty(self):
        # Defensive — caller bug shouldn't blow up here.
        self.assertEqual(_split_into_slices(5, 0), [])


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
class ParallelWriteParityTest(unittest.TestCase):
    """End-to-end: a parallel-written cohort BCF must be content-
    equivalent to a serial-written one for the same input."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(
            prefix="cohort_parallel_parity_"))
        # 8 people across 12 sites with a mix of carrier shapes —
        # all-hom-ref, sparse, dense, multi-allelic on the boundary.
        cls.samples = [f"S{i:03d}" for i in range(8)]
        cls.sites = [
            _site(1000, "A", "G", ["0|0", "0|1", "1|1", "0|0",
                                   "0|0", "1|0", "0|1", "1|1"],
                  site_id="rs1"),
            _site(2000, "C", "T", ["0|0", "0|0", "0|1", "1|1",
                                   "0|1", "0|0", "1|1", "0|0"]),
            _site(3000, "G", "A", ["1|1", "1|1", "0|1", "0|0",
                                   "1|1", "0|1", "0|0", "1|0"]),
            _site(4000, "T", "C", ["0|0", "0|0", "0|0", "0|0",
                                   "0|0", "0|1", "0|0", "0|0"],
                  site_id="rs2"),
            _site(5000, "A", "T", ["0|1", "0|1", "0|1", "0|1",
                                   "0|1", "0|1", "0|1", "0|1"]),
        ]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def _serial_reference(self) -> Path:
        """Write a serial-reference BCF to compare parallel runs against."""
        out = self.tmpdir / "serial.bcf"
        with CohortBcfWriter(out, "GRCh38", self.samples) as w:
            w.write_sites(self.sites)
        return out

    def test_workers_1_falls_back_to_serial_path(self):
        # workers=1 should produce byte-identical output to the
        # original CohortBcfWriter — that's the no-parallelism
        # backstop.
        ref = self._serial_reference()
        out = self.tmpdir / "w1.bcf"
        write_cohort_bcf_parallel(
            out, "GRCh38", self.samples, self.sites, workers=1)
        self.assertEqual(_bcf_data_md5(out), _bcf_data_md5(ref))

    def test_workers_2_matches_serial(self):
        ref = self._serial_reference()
        out = self.tmpdir / "w2.bcf"
        write_cohort_bcf_parallel(
            out, "GRCh38", self.samples, self.sites, workers=2)
        self.assertEqual(_bcf_data_md5(out), _bcf_data_md5(ref))
        # Sample column order must match the input sample list — the
        # parallel-merge collapse must not re-shuffle samples.
        self.assertEqual(_bcf_sample_columns(out), self.samples)

    def test_workers_4_matches_serial(self):
        ref = self._serial_reference()
        out = self.tmpdir / "w4.bcf"
        write_cohort_bcf_parallel(
            out, "GRCh38", self.samples, self.sites, workers=4)
        self.assertEqual(_bcf_data_md5(out), _bcf_data_md5(ref))
        self.assertEqual(_bcf_sample_columns(out), self.samples)

    def test_workers_8_matches_serial(self):
        # workers == n_samples — each worker writes exactly one
        # sample's column.
        ref = self._serial_reference()
        out = self.tmpdir / "w8.bcf"
        write_cohort_bcf_parallel(
            out, "GRCh38", self.samples, self.sites, workers=8)
        self.assertEqual(_bcf_data_md5(out), _bcf_data_md5(ref))
        self.assertEqual(_bcf_sample_columns(out), self.samples)

    def test_workers_exceeds_n_falls_back_cleanly(self):
        # workers=16 with 8 samples — orchestrator should clamp
        # to len(slices) and still produce a valid BCF identical to
        # the serial reference.
        ref = self._serial_reference()
        out = self.tmpdir / "w16.bcf"
        write_cohort_bcf_parallel(
            out, "GRCh38", self.samples, self.sites, workers=16)
        self.assertEqual(_bcf_data_md5(out), _bcf_data_md5(ref))
        self.assertEqual(_bcf_sample_columns(out), self.samples)

    def test_partials_dir_is_cleaned_on_success(self):
        # Mid-run the writer creates ``cohort_dir/.partials/<stem>/``;
        # on success it must be removed so a resume doesn't see stale
        # partials. The .partials/ container also gets cleaned if it
        # ends up empty.
        out = self.tmpdir / "cleanup_check.bcf"
        write_cohort_bcf_parallel(
            out, "GRCh38", self.samples, self.sites, workers=4)
        self.assertFalse((self.tmpdir / ".partials").exists(),
                         ".partials directory should be cleaned up "
                         "after a successful parallel write")


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
class ParallelWriteCohortAFTest(unittest.TestCase):
    """The per-site INFO ``AN`` and ``AF`` must reflect the *full*
    cohort, not just the worker's sample slice. This is the regression
    boundary for the ``cohort_size`` arg threading through
    ``CohortBcfWriter``."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(
            prefix="cohort_parallel_af_"))
        # 4 people, 1 site. Site has AC=3 (one 1|1 + one 0|1 = 3 alt
        # haplotypes). Full cohort AN=8 → AF=0.375.
        cls.samples = ["S0", "S1", "S2", "S3"]
        cls.sites = [
            _site(1000, "A", "G", ["1|1", "0|1", "0|0", "0|0"],
                  site_id="rs1"),
        ]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_an_af_reflect_full_cohort_not_slice(self):
        # workers=2 → slices [(0, 2), (2, 4)]. If a worker wrote AN
        # using its own slice size (4 instead of 8), bcftools merge
        # would surface the mismatch. Pin AN=8 and AF≈0.375 in the
        # final merged BCF.
        out = self.tmpdir / "af_check.bcf"
        write_cohort_bcf_parallel(
            out, "GRCh38", self.samples, self.sites, workers=2)
        info = subprocess.check_output(
            ["bcftools", "query", "-f",
             "%INFO/AN\t%INFO/AF\n", str(out)],
            text=True,
        ).strip().splitlines()
        self.assertEqual(len(info), 1)
        an_str, af_str = info[0].split("\t")
        self.assertEqual(int(an_str), 8)
        self.assertAlmostEqual(float(af_str), 0.375, places=3)


if __name__ == "__main__":
    unittest.main()
