"""Tests for cohort_dispatch — Phase 5g.4 single-pass dispatch.

The headline contract is **byte-equivalence with
``derive_persons_batch``**: dispatch+read must produce identical
per-person record lists for the same cohort BCFs and sample IDs.
That gates everything else in Phase 5g.4, so it gets the most test
weight here.

Additional standing hazards from the plan are pinned individually:

* M13.5 MT carve-out (every MT record reaches every sample
  regardless of original GT).
* ``afs=[None]`` shape (write_person_vcf's MT lineage-carrier
  fallback depends on it).
* Full INFO field carriage (CLNSIG / CLNDN / COSMIC_ID /
  COSMIC_GENE / SVTYPE / SVLEN / END / CIPOS).
* ``cipos`` tuple type preservation across the JSON round-trip.
* Record ordering matches ``derive_persons_batch`` across multi-
  chromosome cohort BCF lists.

Tests gate on bcftools/bgzip on PATH (matching the existing
``test_cohort_derivation.py`` convention).
"""

from __future__ import annotations

import shutil
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.bcf_writer import CohortBcfWriter
from syntheticgen.cohort_derivation import derive_persons_batch
from syntheticgen.cohort_dispatch import (
    dispatch_cohort_to_staging,
    read_person_staging,
)
from syntheticgen.cohort_sites import carriers_from_dense_gts


_HAVE_BCFTOOLS = shutil.which("bcftools") is not None


def _site(pos: int, ref: str, alt: str, gts: list,
          chrom: str = "22", site_id: str = ".",
          **extra) -> dict:
    """Build a cohort site dict the same way
    ``test_cohort_derivation`` does — keeps the two test suites
    comparable site-for-site.
    """
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


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
class DispatchByteEquivalenceTest(unittest.TestCase):
    """Load-bearing test: dispatch+read produces identical per-person
    record lists to ``derive_persons_batch`` for the same inputs.

    This is the gating property for Phase 5g.4. A regression here
    means the dispatch path silently produces different data — the
    whole point of the new path is byte-equivalence, so this test
    must stay green for the PR to merge.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(prefix="cohort_dispatch_be_"))
        cls.samples = ["alice", "bob", "carol", "dave"]
        # Mix of all-hom-ref / mixed / all-carry sites to exercise
        # every per-record dispatch path.
        cls.sites = [
            _site(1000, "A", "G", ["0|0", "0|1", "1|1", "0|0"],
                  site_id="rs1"),
            _site(2000, "C", "T", ["0|0", "0|0", "0|1", "1|0"]),
            _site(
                3000, "G", "<DEL>",
                ["0|0", "0|0", "0|1", "0|1"], site_id="rs2",
                svtype="DEL", svlen=-500, end=3500, cipos=(-50, 50),
                clnsig="Pathogenic", clndn="bar",
            ),
            _site(4000, "T", "A", ["0|0", "1|1", "0|0", "0|0"],
                  site_id="rs3"),
            _site(5000, "G", "C", ["0|1", "0|1", "0|1", "0|1"]),
        ]
        cls.bcf = cls.tmpdir / "cohort.chr22.bcf"
        with CohortBcfWriter(cls.bcf, "GRCh38", cls.samples) as w:
            w.write_sites(cls.sites)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_each_sample_matches_derive_persons_batch(self):
        # Byte-equivalence: per-sample record lists from the two
        # paths must be element-for-element equal.
        batched = derive_persons_batch([self.bcf], self.samples)
        staging_dir = self.tmpdir / "staging_be"
        paths = dispatch_cohort_to_staging(
            [self.bcf], self.samples, staging_dir,
        )
        for sid in self.samples:
            with self.subTest(sample=sid):
                dispatched = read_person_staging(paths[sid])
                self.assertEqual(
                    dispatched, batched[sid],
                    f"dispatch and batched records diverge for {sid}",
                )

    def test_returns_path_for_every_sample(self):
        staging_dir = self.tmpdir / "staging_paths"
        paths = dispatch_cohort_to_staging(
            [self.bcf], self.samples, staging_dir,
        )
        self.assertEqual(set(paths.keys()), set(self.samples))
        for sid in self.samples:
            self.assertTrue(
                paths[sid].exists(),
                f"staging file missing for {sid}",
            )


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
class DispatchMtCarveOutTest(unittest.TestCase):
    """M13.5 contract: every MT record must reach every sample's
    staging file regardless of original simulator GT, so the
    write-time lineage clonality override sees the same MT record
    set across same-lineage members.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(prefix="cohort_dispatch_mt_"))
        cls.samples = ["alice", "bob", "carol"]
        # Two MT records: one where everyone is hom-ref (the case
        # the carve-out exists for), one where carol is the only
        # carrier. Both must land in every sample's staging.
        cls.sites = [
            _site(100, "A", "G", ["0|0", "0|0", "0|0"],
                  chrom="MT", site_id="mt_homref"),
            _site(200, "T", "C", ["0|0", "0|0", "1|1"],
                  chrom="MT", site_id="mt_carol_alt"),
        ]
        cls.bcf = cls.tmpdir / "cohort.MT.bcf"
        with CohortBcfWriter(cls.bcf, "GRCh38", cls.samples) as w:
            w.write_sites(cls.sites)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_every_mt_record_lands_in_every_sample_staging(self):
        staging_dir = self.tmpdir / "staging_mt"
        paths = dispatch_cohort_to_staging(
            [self.bcf], self.samples, staging_dir,
        )
        for sid in self.samples:
            with self.subTest(sample=sid):
                recs = read_person_staging(paths[sid])
                positions = [(r["chrom"], r["pos"]) for r in recs]
                self.assertEqual(
                    positions,
                    [("MT", 100), ("MT", 200)],
                    f"MT record set mismatched for {sid}",
                )

    def test_mt_carve_out_matches_derive_persons_batch(self):
        # Cross-check: the in-memory derive_persons_batch must do the
        # same carve-out; dispatch's behaviour must match site-for-
        # site on MT.
        batched = derive_persons_batch([self.bcf], self.samples)
        staging_dir = self.tmpdir / "staging_mt_xcheck"
        paths = dispatch_cohort_to_staging(
            [self.bcf], self.samples, staging_dir,
        )
        for sid in self.samples:
            with self.subTest(sample=sid):
                self.assertEqual(
                    read_person_staging(paths[sid]),
                    batched[sid],
                    f"MT carve-out diverges between paths for {sid}",
                )


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
class DispatchInfoRoundTripTest(unittest.TestCase):
    """Full ``_INFO_FIELDS_TO_CARRY`` set must round-trip through
    staging with the right value types — string for CLNSIG / CLNDN /
    COSMIC_* / SVTYPE, int for SVLEN / END, ``tuple[int, int]`` for
    CIPOS.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(prefix="cohort_dispatch_info_"))
        cls.samples = ["carrier"]
        cls.sites = [
            _site(
                3000, "G", "<DEL>", ["0|1"], site_id="rs_sv",
                svtype="DEL", svlen=-1234, end=4234,
                cipos=(-25, 75),
                clnsig="Pathogenic", clndn="DiseaseName",
            ),
            _site(
                7000, "C", "T", ["1|1"], site_id="rs_cosmic",
                cosmic_id="COSV12345", cosmic_gene="TP53",
            ),
        ]
        cls.bcf = cls.tmpdir / "cohort.chr22.bcf"
        with CohortBcfWriter(cls.bcf, "GRCh38", cls.samples) as w:
            w.write_sites(cls.sites)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_sv_metadata_round_trips_with_correct_types(self):
        staging_dir = self.tmpdir / "staging_info_sv"
        paths = dispatch_cohort_to_staging(
            [self.bcf], self.samples, staging_dir,
        )
        recs = read_person_staging(paths["carrier"])
        sv = next(r for r in recs if r["pos"] == 3000)
        self.assertEqual(sv["svtype"], "DEL")
        self.assertEqual(sv["svlen"], -1234)
        self.assertEqual(sv["end"], 4234)
        # cipos must be a tuple, not a list — write_person_vcf may
        # rely on the type elsewhere; the in-memory derivation emits
        # tuples and the dispatch path must too.
        self.assertEqual(sv["cipos"], (-25, 75))
        self.assertIsInstance(sv["cipos"], tuple)
        self.assertEqual(sv["clnsig"], "Pathogenic")
        self.assertEqual(sv["clndn"], "DiseaseName")

    def test_cosmic_metadata_round_trips(self):
        staging_dir = self.tmpdir / "staging_info_cos"
        paths = dispatch_cohort_to_staging(
            [self.bcf], self.samples, staging_dir,
        )
        recs = read_person_staging(paths["carrier"])
        cos = next(r for r in recs if r["pos"] == 7000)
        self.assertEqual(cos["cosmic_id"], "COSV12345")
        self.assertEqual(cos["cosmic_gene"], "TP53")


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
class DispatchAfsShapeTest(unittest.TestCase):
    """``write_person_vcf``'s MT lineage-carrier fallback depends on
    the ``afs=[None] * len(alts)`` shape that
    ``derive_persons_batch`` emits. Dispatch must preserve it.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(prefix="cohort_dispatch_afs_"))
        cls.samples = ["s1"]
        cls.sites = [
            _site(1000, "A", "G", ["1|1"]),
            _site(2000, "C", "T", ["0|1"]),
        ]
        cls.bcf = cls.tmpdir / "cohort.chr22.bcf"
        with CohortBcfWriter(cls.bcf, "GRCh38", cls.samples) as w:
            w.write_sites(cls.sites)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_every_record_has_afs_none_per_alt(self):
        staging_dir = self.tmpdir / "staging_afs"
        paths = dispatch_cohort_to_staging(
            [self.bcf], self.samples, staging_dir,
        )
        recs = read_person_staging(paths["s1"])
        self.assertEqual(len(recs), 2)
        for r in recs:
            self.assertEqual(
                r["afs"], [None] * len(r["alts"]),
                f"afs shape diverged for record {r!r}",
            )


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
class DispatchMultiChromOrderingTest(unittest.TestCase):
    """Records emerge per-sample in the BCF list's order — same as
    ``derive_persons_batch``."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = Path(tempfile.mkdtemp(prefix="cohort_dispatch_mc_"))
        cls.samples = ["S1", "S2"]
        cls.bcfs = []
        # Build cohort BCFs for chroms 20, 21, 22 — both samples
        # carry an alt at one site per chrom so ordering is
        # observable.
        for chrom in ("20", "21", "22"):
            sites = [_site(1000, "A", "G", ["0|1", "0|1"],
                           chrom=chrom, site_id=f"rs_{chrom}")]
            bcf = cls.tmpdir / f"cohort.chr{chrom}.bcf"
            with CohortBcfWriter(bcf, "GRCh38", cls.samples) as w:
                w.write_sites(sites)
            cls.bcfs.append(bcf)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_records_appear_in_bcf_list_order(self):
        staging_dir = self.tmpdir / "staging_mc"
        paths = dispatch_cohort_to_staging(
            self.bcfs, self.samples, staging_dir,
        )
        for sid in self.samples:
            with self.subTest(sample=sid):
                recs = read_person_staging(paths[sid])
                self.assertEqual(
                    [r["chrom"] for r in recs],
                    ["20", "21", "22"],
                )

    def test_matches_derive_persons_batch_across_chroms(self):
        batched = derive_persons_batch(self.bcfs, self.samples)
        staging_dir = self.tmpdir / "staging_mc_xcheck"
        paths = dispatch_cohort_to_staging(
            self.bcfs, self.samples, staging_dir,
        )
        for sid in self.samples:
            with self.subTest(sample=sid):
                self.assertEqual(
                    read_person_staging(paths[sid]),
                    batched[sid],
                )


@unittest.skipUnless(_HAVE_BCFTOOLS, "bcftools not on PATH")
class DispatchEdgeCaseTest(unittest.TestCase):
    """Empty / unrelated / missing-input cases."""

    def test_empty_sample_list_returns_empty_dict(self):
        # No samples → no scan, no staging, empty dict. Matches
        # derive_persons_batch's empty-list short-circuit.
        with tempfile.TemporaryDirectory() as td:
            out = dispatch_cohort_to_staging([], [], Path(td))
            self.assertEqual(out, {})

    def test_unrelated_sample_returns_empty_staging(self):
        # A sample that's hom-ref at every site has an empty staging
        # file. The staging file still exists (the contract is
        # "every sample_id gets a path"), it's just zero bytes.
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            samples = ["carrier", "noncarrier"]
            sites = [_site(1000, "A", "G", ["0|1", "0|0"])]
            bcf = tmp / "cohort.chr22.bcf"
            with CohortBcfWriter(bcf, "GRCh38", samples) as w:
                w.write_sites(sites)
            paths = dispatch_cohort_to_staging(
                [bcf], samples, tmp / "staging",
            )
            self.assertEqual(
                read_person_staging(paths["noncarrier"]), [])
            self.assertTrue(paths["noncarrier"].exists())

    def test_missing_bcf_raises_runtime_error(self):
        # Mirror cohort_derivation's error-path contract: bcftools
        # failure surfaces as RuntimeError with the captured stderr.
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaisesRegex(RuntimeError, "bcftools"):
                dispatch_cohort_to_staging(
                    ["/nonexistent/cohort.chr22.bcf"],
                    ["S1"], Path(td) / "staging",
                )


if __name__ == "__main__":
    unittest.main()
