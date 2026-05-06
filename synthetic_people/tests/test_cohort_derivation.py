"""Tests for cohort_derivation — per-person records from cohort BCFs.

The Phase 5b2 derivation pipeline replaces the in-memory
``person_records_from_cohort`` for callers that source their cohort
sites from disk-backed BCFs (Phase 5b1's streamed output). This file
covers the contract:

- The per-person record dict shape is identical to what
  ``person_records_from_cohort`` returns (so ``write_person_vcf``
  consumes either source identically).
- Hom-ref records are dropped, missing (``./.``) records are kept —
  matching the in-memory function's drop-only-all-zero semantics.
- Overlay metadata (CLNSIG / COSMIC / SVTYPE / SVLEN / CIPOS) flows
  from the BCF's INFO field through to the per-person record dict
  with the right key shapes.

Tests gate on bcftools/bgzip on PATH (matching test_bcf_writer.py).
"""

from __future__ import annotations

import shutil
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.bcf_writer import CohortBcfWriter
from syntheticgen.cohort import person_records_from_cohort
from syntheticgen.cohort_derivation import derive_person_records
from syntheticgen.cohort_sites import carriers_from_dense_gts


_HAVE_BCFTOOLS = shutil.which("bcftools") is not None


def _site(pos: int, ref: str, alt: str, gts: list,
          chrom: str = "22", site_id: str = ".",
          **extra) -> dict:
    nalt = sum(int(t) for gt in gts for t in gt.split("|"))
    # Phase 5c: site stores sparse carriers — derived from the
    # human-readable dense GTs the test passed in. The BCF writer's
    # dense-gts fallback accepts either form, but
    # person_records_from_cohort now requires carriers, so we go
    # through the sparse shape here for both branches to share a
    # fixture.
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


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
class CohortDerivationParityTest(unittest.TestCase):
    """Round-trip parity vs in-memory ``person_records_from_cohort``.

    Build a cohort, write it to a BCF, derive each sample from the
    BCF, and compare the per-person record list to what the in-memory
    function returns from the same cohort. Per-record fields must
    match for chrom / pos / id / ref / alts / gt and for any overlay
    metadata (CLNSIG / SVTYPE etc.) that's present.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(prefix="cohort_deriv_parity_"))
        cls.samples = ["alice", "bob", "carol"]
        cls.sites = [
            _site(1000, "A", "G", ["0|0", "0|1", "1|1"], site_id="rs1"),
            _site(2000, "C", "T", ["0|0", "0|1", "0|0"]),
            _site(
                3000, "G", "<DEL>", ["0|0", "0|0", "0|1"], site_id="rs2",
                svtype="DEL", svlen=-500, end=3500, cipos=(-50, 50),
                clnsig="Pathogenic", clndn="foo",
            ),
        ]
        cls.bcf = cls.tmpdir / "cohort.chr22.bcf"
        with CohortBcfWriter(cls.bcf, "GRCh38", cls.samples) as w:
            w.write_sites(cls.sites)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def _ref_records(self, person_index: int) -> list:
        return person_records_from_cohort(self.sites, person_index)

    def test_alice_drops_all_hom_ref(self):
        # Alice is 0|0 at every site → in-memory derivation returns
        # an empty list; the BCF derivation must do the same.
        derived = derive_person_records([self.bcf], "alice")
        self.assertEqual(derived, [])
        self.assertEqual(self._ref_records(0), [])

    def test_bob_record_set_matches_in_memory(self):
        derived = derive_person_records([self.bcf], "bob")
        ref = self._ref_records(1)
        # Same record count.
        self.assertEqual(len(derived), len(ref))
        # Per-record: chrom/pos/id/ref/alts/gt agree.
        for d, r in zip(derived, ref):
            for key in ("chrom", "pos", "id", "ref", "alts", "gt"):
                self.assertEqual(
                    d[key], r[key],
                    msg=f"key={key} d={d} r={r}")

    def test_carol_overlay_metadata_round_trips(self):
        derived = derive_person_records([self.bcf], "carol")
        # Carol has the SV record at pos 3000 — it should carry every
        # overlay field through the BCF round-trip.
        sv_rec = next(r for r in derived if r["pos"] == 3000)
        self.assertEqual(sv_rec["svtype"], "DEL")
        self.assertEqual(sv_rec["svlen"], -500)
        self.assertEqual(sv_rec["end"], 3500)
        self.assertEqual(sv_rec["cipos"], (-50, 50))
        self.assertEqual(sv_rec["clnsig"], "Pathogenic")
        self.assertEqual(sv_rec["clndn"], "foo")


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
class CohortDerivationMultiChromTest(unittest.TestCase):
    """When deriving across multiple per-chrom BCFs, records emerge in
    BCF-iteration order (chr1 first, then chr2, etc. — same order the
    caller passes the path list)."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(prefix="cohort_deriv_multi_"))
        cls.samples = ["S1", "S2"]
        # One BCF per chromosome with a single alt-bearing record per
        # chrom for S1.
        cls.bcfs = []
        for i, chrom in enumerate(("21", "22"), start=1):
            sites = [_site(1000 * i, "A", "G", ["0|1", "0|0"],
                           chrom=chrom, site_id=f"rs{chrom}")]
            bcf = cls.tmpdir / f"cohort.chr{chrom}.bcf"
            with CohortBcfWriter(bcf, "GRCh38", cls.samples) as w:
                w.write_sites(sites)
            cls.bcfs.append(bcf)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_records_come_back_in_bcf_order(self):
        recs = derive_person_records(self.bcfs, "S1")
        self.assertEqual([r["chrom"] for r in recs], ["21", "22"])

    def test_unrelated_sample_still_returns_empty(self):
        recs = derive_person_records(self.bcfs, "S2")
        self.assertEqual(recs, [])


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
class CohortDerivationErrorPathTest(unittest.TestCase):
    """A missing or corrupt BCF should surface a clear error rather
    than silently producing an empty list."""

    def test_missing_bcf_raises_runtime_error(self):
        with self.assertRaisesRegex(RuntimeError, "bcftools"):
            derive_person_records(["/nonexistent/cohort.chr22.bcf"], "S1")


if __name__ == "__main__":
    unittest.main()
