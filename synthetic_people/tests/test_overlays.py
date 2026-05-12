"""Tests for M7 overlays (ClinVar annotate/inject, dbSNP rsID inject,
COSMIC inject).

These tests build pure in-memory pools / site lists, so they don't
depend on bcftools or downloaded VCFs and run fast in any environment.
The bcftools-driven `load_*` functions are exercised by the CLI smoke
test in the exit check rather than here.
"""

from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.clinvar import annotate_clinvar, inject_clinvar
from syntheticgen.cosmic import inject_cosmic
from syntheticgen.dbsnp import _normalise_rsid, inject_rsids


def _site(chrom: str, pos: int, ref: str, alt: str,
          gts: list | None = None) -> dict:
    """Build a cohort site dict matching the post-Phase-5c shape.

    Tests still pass a human-readable list of GT strings;
    ``carriers_from_dense_gts`` converts to the sparse storage shape
    that production code emits, so overlay round-trips exercise the
    same data path.
    """
    from syntheticgen.cohort_sites import carriers_from_dense_gts
    if gts is None:
        gts = ["0|0", "0|1", "1|1"]
    nalt = sum(int(t) for gt in gts for t in gt.split("|"))
    n_haps = 2 * len(gts)
    return {
        "chrom": chrom, "pos": pos, "id": ".",
        "ref": ref, "alts": [alt],
        "afs": [nalt / n_haps], "acs": [nalt],
        "n_haplotypes": n_haps,
        "carriers": carriers_from_dense_gts(gts),
    }


def _clinvar_rec(chrom: str, pos: int, ref: str, alt: str,
                 vid: str = "VCV1", clnsig: str = "Pathogenic",
                 clndn: str = "Disease_X", rsid: str = "") -> dict:
    return {
        "chrom": chrom, "pos": pos, "id": vid,
        "ref": ref, "alt": alt,
        "clnsig": clnsig, "clndn": clndn, "rsid": rsid,
    }


class TestAnnotateClinvar(unittest.TestCase):
    def test_annotates_collisions(self):
        sites = [
            _site("22", 100, "A", "G"),
            _site("22", 200, "C", "T"),
            _site("22", 300, "G", "A"),
        ]
        records = [
            _clinvar_rec("22", 200, "C", "T", "VCV200",
                         clnsig="Pathogenic", clndn="Test"),
        ]
        n = annotate_clinvar(sites, records)
        self.assertEqual(n, 1)
        self.assertEqual(sites[1]["clnsig"], "Pathogenic")
        self.assertEqual(sites[1]["clndn"], "Test")
        self.assertEqual(sites[1]["id"], "VCV200")
        # Non-matching sites untouched
        self.assertNotIn("clnsig", sites[0])
        self.assertNotIn("clnsig", sites[2])

    def test_no_match_returns_zero(self):
        sites = [_site("22", 100, "A", "G")]
        records = [_clinvar_rec("22", 999, "A", "G")]
        self.assertEqual(annotate_clinvar(sites, records), 0)
        self.assertNotIn("clnsig", sites[0])

    def test_alt_mismatch_does_not_annotate(self):
        sites = [_site("22", 100, "A", "G")]
        records = [_clinvar_rec("22", 100, "A", "T")]
        self.assertEqual(annotate_clinvar(sites, records), 0)


class TestInjectClinvar(unittest.TestCase):
    def test_injects_density_records(self):
        sites = [_site("22", i * 100, "A", "G") for i in range(1, 11)]
        records = [_clinvar_rec("22", 1_000_000 + i, "C", "T",
                                vid=f"VCV{i}") for i in range(20)]
        rng = random.Random(0)
        n = inject_clinvar(sites, records, density=0.3, rng=rng)
        self.assertEqual(n, 3)
        # Three sites carry CLNSIG and have moved to ClinVar coordinates
        with_clnsig = [s for s in sites if "clnsig" in s]
        self.assertEqual(len(with_clnsig), 3)
        for s in with_clnsig:
            self.assertGreaterEqual(s["pos"], 1_000_000)
            self.assertEqual(s["ref"], "C")
            self.assertEqual(s["alts"], ["T"])
            self.assertTrue(s["id"].startswith("VCV"))

    def test_preserves_gt_block(self):
        from syntheticgen.cohort_sites import carriers_from_dense_gts
        gts = ["0|0", "0|1", "1|1", "0|0"]
        expected_carriers = carriers_from_dense_gts(gts)
        sites = [_site("22", i * 100, "A", "G", gts=gts)
                 for i in range(1, 6)]
        records = [_clinvar_rec("22", 1_000_000 + i, "C", "T")
                   for i in range(10)]
        rng = random.Random(1)
        inject_clinvar(sites, records, density=0.4, rng=rng)
        # Every site still carries the original GT block — sparse
        # carriers under Phase 5c, but representing the same per-
        # person genotypes.
        for s in sites:
            self.assertTrue(
                np.array_equal(s["carriers"], expected_carriers),
            )

    def test_keeps_sites_sorted(self):
        sites = [_site("22", i * 100, "A", "G") for i in range(1, 21)]
        records = [_clinvar_rec("22", 5_000_000 + i, "C", "T")
                   for i in range(50)]
        rng = random.Random(2)
        inject_clinvar(sites, records, density=0.5, rng=rng)
        positions = [s["pos"] for s in sites]
        self.assertEqual(positions, sorted(positions))

    def test_zero_density_noop(self):
        sites = [_site("22", 100, "A", "G")]
        records = [_clinvar_rec("22", 200, "C", "T")]
        rng = random.Random(0)
        self.assertEqual(inject_clinvar(sites, records, 0.0, rng), 0)
        self.assertEqual(sites[0]["pos"], 100)
        self.assertEqual(sites[0]["alts"], ["G"])

    def test_skips_records_off_chromosome(self):
        sites = [_site("22", i * 100, "A", "G") for i in range(1, 6)]
        records = [_clinvar_rec("21", 5_000, "C", "T") for _ in range(10)]
        rng = random.Random(3)
        n = inject_clinvar(sites, records, density=0.5, rng=rng)
        self.assertEqual(n, 0)
        # All sites remain unchanged.
        for s in sites:
            self.assertNotIn("clnsig", s)


class TestNormaliseRsid(unittest.TestCase):
    def test_id_column_already_prefixed(self):
        self.assertEqual(_normalise_rsid("rs12345", ""), "rs12345")

    def test_id_column_bare_digits(self):
        self.assertEqual(_normalise_rsid("12345", ""), "rs12345")

    def test_clinvar_info_rs_bare_digits(self):
        self.assertEqual(_normalise_rsid(".", "98765"), "rs98765")

    def test_id_preferred_over_info(self):
        self.assertEqual(_normalise_rsid("rs1", "rs2"), "rs1")

    def test_missing_returns_empty(self):
        self.assertEqual(_normalise_rsid(".", "."), "")
        self.assertEqual(_normalise_rsid("", ""), "")

    def test_first_of_semicolon_list(self):
        self.assertEqual(_normalise_rsid("rs1;rs2", ""), "rs1")

    def test_first_of_comma_info(self):
        self.assertEqual(_normalise_rsid(".", "111,222"), "rs111")


class TestInjectRsids(unittest.TestCase):
    def _pool(self, n: int = 10, chrom: str = "22") -> list:
        return [
            {"chrom": chrom, "pos": 2_000_000 + i,
             "ref": "A", "alt": "G", "rsid": f"rs{i+1000}"}
            for i in range(n)
        ]

    def test_injects_density_rsids(self):
        sites = [_site("22", i * 100, "A", "G") for i in range(1, 11)]
        rng = random.Random(0)
        n = inject_rsids(sites, self._pool(), density=0.3, rng=rng)
        self.assertEqual(n, 3)
        with_id = [s for s in sites if s["id"].startswith("rs")]
        self.assertEqual(len(with_id), 3)

    def test_preserves_gt_block(self):
        from syntheticgen.cohort_sites import carriers_from_dense_gts
        gts = ["0|1", "1|0", "0|0"]
        expected_carriers = carriers_from_dense_gts(gts)
        sites = [_site("22", i * 100, "A", "G", gts=gts)
                 for i in range(1, 6)]
        rng = random.Random(1)
        inject_rsids(sites, self._pool(), density=0.6, rng=rng)
        for s in sites:
            self.assertTrue(
                np.array_equal(s["carriers"], expected_carriers),
            )

    def test_reserve_indices_excluded(self):
        sites = [_site("22", i * 100, "A", "G") for i in range(1, 6)]
        rng = random.Random(2)
        # Reserve the first three indices — only two candidates remain
        # so the dbSNP injection should land on at most those two.
        n = inject_rsids(sites, self._pool(),
                         density=1.0, rng=rng,
                         reserve_indices={0, 1, 2})
        # After sort the original indices may have moved, but the count
        # of rsID-bearing rows should equal the injection count.
        self.assertLessEqual(n, 2)
        rsid_count = sum(1 for s in sites if s["id"].startswith("rs"))
        self.assertEqual(rsid_count, n)

    def test_keeps_sorted(self):
        sites = [_site("22", i * 100, "A", "G") for i in range(1, 21)]
        rng = random.Random(3)
        inject_rsids(sites, self._pool(20), density=0.5, rng=rng)
        positions = [s["pos"] for s in sites]
        self.assertEqual(positions, sorted(positions))

    def test_zero_density_noop(self):
        sites = [_site("22", 100, "A", "G")]
        rng = random.Random(0)
        n = inject_rsids(sites, self._pool(), density=0, rng=rng)
        self.assertEqual(n, 0)
        self.assertEqual(sites[0]["id"], ".")


class TestInjectCosmic(unittest.TestCase):
    def test_injects_with_gene_and_id(self):
        sites = [_site("22", i * 100, "A", "G") for i in range(1, 11)]
        records = [{
            "chrom": "22", "pos": 4_000_000 + i,
            "ref": "C", "alt": "T",
            "id": f"COSV{i}",
            "gene": "TP53", "cds": "", "aa": "",
        } for i in range(10)]
        rng = random.Random(0)
        n = inject_cosmic(sites, records, density=0.3, rng=rng)
        self.assertEqual(n, 3)
        with_cosmic = [s for s in sites if s.get("cosmic_id")]
        self.assertEqual(len(with_cosmic), 3)
        for s in with_cosmic:
            self.assertEqual(s["cosmic_gene"], "TP53")
            self.assertTrue(s["cosmic_id"].startswith("COSV"))
            self.assertEqual(s["ref"], "C")
            self.assertEqual(s["alts"], ["T"])

    def test_zero_density_noop(self):
        sites = [_site("22", 100, "A", "G")]
        records = [{"chrom": "22", "pos": 200, "ref": "C", "alt": "T",
                    "id": "COSV1", "gene": "X", "cds": "", "aa": ""}]
        rng = random.Random(0)
        self.assertEqual(inject_cosmic(sites, records, 0, rng), 0)
        self.assertNotIn("cosmic_id", sites[0])


class TestOverlayInteraction(unittest.TestCase):
    """Ensure ClinVar + rsID injection don't fight for the same rows."""

    def test_reserve_indices_keeps_clinvar_rows_intact(self):
        sites = [_site("22", i * 100, "A", "G") for i in range(1, 21)]
        clinvar_records = [_clinvar_rec("22", 3_000_000 + i, "C", "T",
                                        vid=f"VCV{i}")
                           for i in range(10)]
        rsid_pool = [
            {"chrom": "22", "pos": 5_000_000 + i,
             "ref": "T", "alt": "A", "rsid": f"rs{i+9000}"}
            for i in range(20)
        ]
        rng = random.Random(7)
        inject_clinvar(sites, clinvar_records, density=0.2, rng=rng)
        clinvar_indices = {i for i, s in enumerate(sites)
                           if s.get("clnsig")}
        inject_rsids(sites, rsid_pool, density=0.3, rng=rng,
                     reserve_indices=clinvar_indices)
        # Every CLNSIG-bearing row should still carry CLNSIG (i.e. the
        # rsID injection didn't overwrite it).
        for s in sites:
            if s.get("clnsig"):
                # ClinVar-injected rows shouldn't have an rsID-style ID
                self.assertFalse(s["id"].startswith("rs"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
