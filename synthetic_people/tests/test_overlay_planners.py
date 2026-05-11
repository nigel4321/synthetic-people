"""Tests for the planner halves of the M7 overlay injectors
(``plan_inject_clinvar`` / ``plan_inject_rsids`` /
``plan_inject_cosmic``).

These are the rng-consuming halves extracted from the in-place
``inject_*`` functions so they can run without a materialised sites
list — the streaming-cohort path uses them directly. The legacy
``inject_*`` functions are thin wrappers that call the planner and
then apply the result in place, so their existing tests in
``test_overlays.py`` already lock the in-place behaviour. These tests
verify two additional invariants on top:

  1. The planner produces a ``{site_index: overlay_record}`` mapping
     consistent with what the legacy in-place mutator would have
     written for the same inputs at the same seed (parity).
  2. The fields each overlay record carries match the per-overlay
     contract (ClinVar emits clnsig/clndn; rsID emits id only;
     COSMIC emits id/cosmic_id/cosmic_gene conditionally).
"""

from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.clinvar import (  # noqa: E402
    inject_clinvar, plan_inject_clinvar,
)
from syntheticgen.cosmic import (  # noqa: E402
    inject_cosmic, plan_inject_cosmic,
)
from syntheticgen.dbsnp import (  # noqa: E402
    inject_rsids, plan_inject_rsids,
)


def _site(chrom: str, pos: int, ref: str = "A", alt: str = "G") -> dict:
    """Minimal site dict — overlays only touch coordinate / id / overlay
    INFO fields, never carriers, so the carriers slot stays empty in
    these tests. Matches the data the planner cares about for its
    decisions (chrom + pos via sites_meta), plus the fields the
    in-place wrapper mutates."""
    return {
        "chrom": chrom, "pos": pos, "id": ".",
        "ref": ref, "alts": [alt], "acs": [1], "afs": [0.1],
        "n_haplotypes": 6, "carriers": [],
    }


def _clinvar_rec(chrom: str, pos: int, ref: str, alt: str,
                 vid: str = "VCV000001", clnsig: str = "Pathogenic",
                 clndn: str = "Some_Disease") -> dict:
    return {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "id": vid, "clnsig": clnsig, "clndn": clndn}


def _rsid_rec(chrom: str, pos: int, ref: str, alt: str,
              rsid: str = "rs1") -> dict:
    return {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "rsid": rsid}


def _cosmic_rec(chrom: str, pos: int, ref: str, alt: str,
                cid: str = "COSV99", gene: str = "TP53") -> dict:
    return {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "id": cid, "gene": gene}


# --- Planner output shape -------------------------------------------


class PlanInjectClinvarShapeTest(unittest.TestCase):
    """The planner returns a {index: overlay_record} dict; each
    overlay_record has exactly the fields the in-place wrapper mutates
    (pos / ref / alts / id / clnsig / clndn)."""

    def test_returns_dict_keyed_by_site_index(self):
        sites_meta = [("22", i * 100) for i in range(1, 11)]
        records = [_clinvar_rec("22", 1_000_000 + i, "C", "T",
                                vid=f"VCV{i}") for i in range(20)]
        plan = plan_inject_clinvar(
            sites_meta, records, density=0.3,
            rng=random.Random(0),
        )
        self.assertEqual(len(plan), 3)
        for idx, rec in plan.items():
            self.assertIsInstance(idx, int)
            self.assertGreaterEqual(idx, 0)
            self.assertLess(idx, len(sites_meta))
            self.assertEqual(
                set(rec.keys()),
                {"pos", "ref", "alts", "id", "clnsig", "clndn"},
            )
            self.assertGreaterEqual(rec["pos"], 1_000_000)
            self.assertEqual(rec["ref"], "C")
            self.assertEqual(rec["alts"], ["T"])
            self.assertTrue(rec["id"].startswith("VCV"))

    def test_zero_density_returns_empty(self):
        sites_meta = [("22", 100)]
        records = [_clinvar_rec("22", 200, "C", "T")]
        plan = plan_inject_clinvar(
            sites_meta, records, 0.0, rng=random.Random(0),
        )
        self.assertEqual(plan, {})

    def test_off_chromosome_records_skipped(self):
        sites_meta = [("22", i * 100) for i in range(1, 6)]
        records = [_clinvar_rec("21", 5_000, "C", "T") for _ in range(10)]
        plan = plan_inject_clinvar(
            sites_meta, records, density=0.5,
            rng=random.Random(3),
        )
        self.assertEqual(plan, {})


class PlanInjectRsidsShapeTest(unittest.TestCase):
    def test_returns_dict_with_id_field(self):
        sites_meta = [("22", i * 100) for i in range(1, 11)]
        pool = [_rsid_rec("22", 1_000_000 + i, "C", "T",
                          rsid=f"rs{i}") for i in range(20)]
        plan = plan_inject_rsids(
            sites_meta, pool, density=0.3,
            rng=random.Random(0),
        )
        self.assertEqual(len(plan), 3)
        for idx, rec in plan.items():
            self.assertEqual(
                set(rec.keys()), {"pos", "ref", "alts", "id"},
            )
            self.assertTrue(rec["id"].startswith("rs"))

    def test_reserve_indices_excluded(self):
        sites_meta = [("22", i * 100) for i in range(1, 11)]
        pool = [_rsid_rec("22", 1_000_000 + i, "C", "T",
                          rsid=f"rs{i}") for i in range(20)]
        reserved = {0, 1, 2}
        plan = plan_inject_rsids(
            sites_meta, pool, density=0.5,
            rng=random.Random(0), reserve_indices=reserved,
        )
        # None of the reserved indices appear in the plan keys.
        for idx in plan:
            self.assertNotIn(idx, reserved)


class PlanInjectCosmicShapeTest(unittest.TestCase):
    def test_conditional_id_gene_fields(self):
        sites_meta = [("22", i * 100) for i in range(1, 11)]
        # Mix of records: with id+gene, with id only, with gene only,
        # and one with id="." (should not set id / cosmic_id).
        pool = [
            _cosmic_rec("22", 1_000_001, "C", "T",
                        cid="COSV1", gene="TP53"),
            _cosmic_rec("22", 1_000_002, "C", "T",
                        cid="COSV2", gene=""),
            _cosmic_rec("22", 1_000_003, "C", "T",
                        cid=".", gene="BRCA1"),
            _cosmic_rec("22", 1_000_004, "C", "T",
                        cid=".", gene=""),
        ]
        plan = plan_inject_cosmic(
            sites_meta, pool, density=1.0,  # try to inject all
            rng=random.Random(0),
        )
        # Records are picked without replacement; up to len(pool).
        self.assertLessEqual(len(plan), len(pool))
        # Each overlay record's optional fields follow the contract.
        for rec in plan.values():
            self.assertIn("pos", rec)
            self.assertIn("ref", rec)
            self.assertIn("alts", rec)
            if "id" in rec:
                self.assertNotEqual(rec["id"], ".")
                self.assertIn("cosmic_id", rec)
                self.assertEqual(rec["id"], rec["cosmic_id"])
            if "cosmic_gene" in rec:
                self.assertNotEqual(rec["cosmic_gene"], "")


# --- Parity with the in-place wrappers ------------------------------


class WrapperParityTest(unittest.TestCase):
    """At a given seed, the legacy in-place ``inject_*`` mutator and a
    plan-then-apply-manually invocation produce identical site dicts.

    Why this matters: the wrapper now goes plan -> apply internally. We
    re-implement the apply path in the test to verify the planner
    output is the single source of truth — any future refactor that
    bypasses ``inject_*`` and uses ``plan_*`` directly (which is
    exactly what the streaming-cohort path will do) gets the same
    bytes."""

    def _baseline_inject_clinvar_then_capture(self, seed: int):
        sites = [_site("22", i * 100) for i in range(1, 11)]
        records = [_clinvar_rec("22", 1_000_000 + i, "C", "T",
                                vid=f"VCV{i}") for i in range(20)]
        inject_clinvar(sites, records, density=0.3,
                       rng=random.Random(seed))
        return sites

    def _planner_apply_manually_then_capture(self, seed: int):
        sites = [_site("22", i * 100) for i in range(1, 11)]
        records = [_clinvar_rec("22", 1_000_000 + i, "C", "T",
                                vid=f"VCV{i}") for i in range(20)]
        sites_meta = [(s["chrom"], s["pos"]) for s in sites]
        plan = plan_inject_clinvar(
            sites_meta, records, density=0.3,
            rng=random.Random(seed),
        )
        for idx, rec in plan.items():
            site = sites[idx]
            site["pos"] = rec["pos"]
            site["ref"] = rec["ref"]
            site["alts"] = rec["alts"]
            site["id"] = rec["id"]
            site["clnsig"] = rec["clnsig"]
            site["clndn"] = rec["clndn"]
        sites.sort(key=lambda s: (s["chrom"], s["pos"]))
        return sites

    def test_clinvar_wrapper_parity(self):
        for seed in (0, 1, 42, 100):
            a = self._baseline_inject_clinvar_then_capture(seed)
            b = self._planner_apply_manually_then_capture(seed)
            self.assertEqual(a, b, f"seed={seed}")

    def test_rsids_wrapper_parity(self):
        for seed in (0, 1, 42, 100):
            sites_a = [_site("22", i * 100) for i in range(1, 11)]
            sites_b = [_site("22", i * 100) for i in range(1, 11)]
            pool = [_rsid_rec("22", 1_000_000 + i, "C", "T",
                              rsid=f"rs{i}") for i in range(20)]

            inject_rsids(sites_a, pool, density=0.3,
                         rng=random.Random(seed))

            sites_meta = [(s["chrom"], s["pos"]) for s in sites_b]
            plan = plan_inject_rsids(
                sites_meta, pool, density=0.3,
                rng=random.Random(seed),
            )
            for idx, rec in plan.items():
                site = sites_b[idx]
                site["pos"] = rec["pos"]
                site["ref"] = rec["ref"]
                site["alts"] = rec["alts"]
                site["id"] = rec["id"]
            sites_b.sort(key=lambda s: (s["chrom"], s["pos"]))

            self.assertEqual(sites_a, sites_b, f"seed={seed}")

    def test_cosmic_wrapper_parity(self):
        for seed in (0, 1, 42, 100):
            sites_a = [_site("22", i * 100) for i in range(1, 11)]
            sites_b = [_site("22", i * 100) for i in range(1, 11)]
            pool = [_cosmic_rec("22", 1_000_000 + i, "C", "T",
                                cid=f"COSV{i}",
                                gene=f"GENE{i % 3}")
                    for i in range(20)]

            inject_cosmic(sites_a, pool, density=0.3,
                          rng=random.Random(seed))

            sites_meta = [(s["chrom"], s["pos"]) for s in sites_b]
            plan = plan_inject_cosmic(
                sites_meta, pool, density=0.3,
                rng=random.Random(seed),
            )
            for idx, rec in plan.items():
                site = sites_b[idx]
                site["pos"] = rec["pos"]
                site["ref"] = rec["ref"]
                site["alts"] = rec["alts"]
                if "id" in rec:
                    site["id"] = rec["id"]
                if "cosmic_id" in rec:
                    site["cosmic_id"] = rec["cosmic_id"]
                if "cosmic_gene" in rec:
                    site["cosmic_gene"] = rec["cosmic_gene"]
            sites_b.sort(key=lambda s: (s["chrom"], s["pos"]))

            self.assertEqual(sites_a, sites_b, f"seed={seed}")


# --- Reserve-indices contract for chained overlays ------------------


class ReserveIndicesChainTest(unittest.TestCase):
    """The cli runs ClinVar -> rsID -> COSMIC, passing each later
    overlay the set of indices the earlier ones already claimed. The
    planners must honour ``reserve_indices`` identically to the
    in-place mutators so the streaming caller can build the same
    chain."""

    def test_rsids_planner_skips_reserved_indices(self):
        sites_meta = [("22", i * 100) for i in range(1, 21)]
        pool = [_rsid_rec("22", 5_000_000 + i, "C", "T",
                          rsid=f"rs{i}") for i in range(40)]
        clinvar_reserved = {3, 7, 11}
        plan = plan_inject_rsids(
            sites_meta, pool, density=0.5,
            rng=random.Random(5),
            reserve_indices=clinvar_reserved,
        )
        for idx in plan:
            self.assertNotIn(idx, clinvar_reserved)

    def test_cosmic_planner_skips_reserved_indices(self):
        sites_meta = [("22", i * 100) for i in range(1, 21)]
        pool = [_cosmic_rec("22", 5_000_000 + i, "C", "T",
                            cid=f"COSV{i}", gene="TP53")
                for i in range(40)]
        all_reserved = {1, 5, 9, 13, 17}
        plan = plan_inject_cosmic(
            sites_meta, pool, density=0.5,
            rng=random.Random(5),
            reserve_indices=all_reserved,
        )
        for idx in plan:
            self.assertNotIn(idx, all_reserved)


if __name__ == "__main__":
    unittest.main()
