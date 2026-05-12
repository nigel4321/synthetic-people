"""Tests for the Arrow IPC streaming intermediate (Phase 5d.1).

Skipped cleanly if pyarrow isn't installed — the cohort-arrow path is
a conditional dep, not a hard one. Install with ``pip install pyarrow``
to exercise these tests.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

try:
    import pyarrow  # noqa: F401
    HAS_PYARROW = True
except ImportError:
    HAS_PYARROW = False

from syntheticgen.cohort_sites import carriers_from_dense_gts  # noqa: E402


def _site_from_gts(pos: int, gts, **extras) -> dict:
    """Build a sparse-carriers site dict from human-readable GT strings.

    Mirrors ``test_cohort._site_with_dense_gts`` in style — readable
    in the test ("0|0", "0|1", "1|1") while storing carriers as
    production code does.
    """
    n_people = len(gts)
    site = {
        "chrom": "20",
        "pos": pos,
        "id": ".",
        "ref": "A",
        "alts": ["C"],
        "acs": [sum(1 for g in gts for a in g.split("|") if a != "0")],
        "afs": [0.0],  # filled in below
        "n_haplotypes": 2 * n_people,
        "carriers": carriers_from_dense_gts(gts),
    }
    if site["n_haplotypes"]:
        site["afs"] = [site["acs"][0] / site["n_haplotypes"]]
    site.update(extras)
    return site


@unittest.skipUnless(HAS_PYARROW, "pyarrow not installed")
class TestCohortSchema(unittest.TestCase):
    def test_metadata_present(self):
        from syntheticgen.cohort_arrow import (
            cohort_schema,
            META_KEY_CHROM,
            META_KEY_N_SAMPLES,
            META_KEY_N_HAPLOTYPES,
            META_KEY_FORMAT_VERSION,
        )

        schema = cohort_schema(n_samples=100, chrom="22")
        meta = schema.metadata or {}
        self.assertEqual(meta.get(META_KEY_CHROM), b"22")
        self.assertEqual(meta.get(META_KEY_N_SAMPLES), b"100")
        self.assertEqual(meta.get(META_KEY_N_HAPLOTYPES), b"200")
        self.assertEqual(meta.get(META_KEY_FORMAT_VERSION), b"1")

    def test_genotypes_field_is_fixed_size_haplotype_list(self):
        from syntheticgen.cohort_arrow import cohort_schema

        schema = cohort_schema(n_samples=50, chrom="1")
        gt_field = schema.field("genotypes")
        self.assertEqual(str(gt_field.type), "fixed_size_list<item: int8>[100]")

    def test_overlay_fields_nullable(self):
        from syntheticgen.cohort_arrow import cohort_schema

        schema = cohort_schema(n_samples=10, chrom="1")
        for name in (
            "clnsig", "clndn", "cosmic_id", "cosmic_gene",
            "svtype", "svlen", "end", "cipos_lo", "cipos_hi",
        ):
            self.assertTrue(
                schema.field(name).nullable,
                f"{name} should be nullable",
            )


@unittest.skipUnless(HAS_PYARROW, "pyarrow not installed")
class TestRoundTrip(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "cohort.arrow"

    def tearDown(self):
        self.tmp.cleanup()

    def test_simple_biallelic_round_trip(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_slice,
        )

        sites = [
            _site_from_gts(100, ["0|0", "0|1", "1|1"]),
            _site_from_gts(200, ["1|0", "0|0", "0|1"]),
            _site_from_gts(300, ["0|0", "0|0", "1|1"]),
        ]
        n_written = write_arrow_file(self.path, "20", 3, iter(sites))
        self.assertEqual(n_written, 3)

        read_back = list(read_arrow_slice(self.path, 0, 3))
        self.assertEqual(len(read_back), 3)

        for orig, got in zip(sites, read_back):
            self.assertEqual(got["chrom"], "20")
            self.assertEqual(got["pos"], orig["pos"])
            self.assertEqual(got["ref"], orig["ref"])
            self.assertEqual(got["alts"], orig["alts"])
            self.assertEqual(got["acs"], orig["acs"])
            self.assertAlmostEqual(got["afs"][0], orig["afs"][0], places=5)
            # gts re-derived from the haplotype matrix should match
            # the original dense-string view of the carriers
            from syntheticgen.cohort_sites import dense_gts_from_carriers
            expected_gts = dense_gts_from_carriers(orig["carriers"], 3)
            self.assertEqual(got["gts"], expected_gts)

    def test_slice_correctness(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_slice,
        )

        sites = [
            _site_from_gts(100, ["0|0", "0|1", "1|1", "1|0", "0|0"]),
            _site_from_gts(200, ["1|1", "0|0", "0|1", "0|1", "1|0"]),
        ]
        write_arrow_file(self.path, "20", 5, iter(sites))

        front = list(read_arrow_slice(self.path, 0, 2))
        back = list(read_arrow_slice(self.path, 2, 5))

        self.assertEqual([s["gts"] for s in front], [
            ["0|0", "0|1"],
            ["1|1", "0|0"],
        ])
        self.assertEqual([s["gts"] for s in back], [
            ["1|1", "1|0", "0|0"],
            ["0|1", "0|1", "1|0"],
        ])

    def test_multi_allelic_round_trip(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_slice,
        )

        sites = [{
            "chrom": "20",
            "pos": 100,
            "id": ".",
            "ref": "A",
            "alts": ["C", "G"],
            "acs": [3, 2],
            "afs": [0.3, 0.2],
            "n_haplotypes": 10,
            "carriers": [(0, 1), (2, 1), (4, 1), (5, 2), (7, 2)],
        }]
        write_arrow_file(self.path, "20", 5, iter(sites))
        read_back = list(read_arrow_slice(self.path, 0, 5))
        self.assertEqual(len(read_back), 1)
        s = read_back[0]
        self.assertEqual(s["alts"], ["C", "G"])
        self.assertEqual(s["acs"], [3, 2])
        self.assertEqual([round(f, 3) for f in s["afs"]], [0.3, 0.2])
        self.assertEqual(s["gts"], ["1|0", "1|0", "1|2", "0|2", "0|0"])

    def test_overlay_fields_round_trip(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_slice,
        )

        site_with = _site_from_gts(
            100, ["0|0", "0|1"],
            id="rs1234",
            clnsig="Pathogenic",
            clndn="Some_Disease",
            cosmic_id="COSV99",
            cosmic_gene="TP53",
            svtype="DEL",
            svlen=-200,
            end=300,
            cipos=(-50, 50),
        )
        site_without = _site_from_gts(200, ["0|0", "1|1"])
        write_arrow_file(self.path, "20", 2, iter([site_with, site_without]))

        read_back = list(read_arrow_slice(self.path, 0, 2))
        s0, s1 = read_back

        self.assertEqual(s0["id"], "rs1234")
        self.assertEqual(s0["clnsig"], "Pathogenic")
        self.assertEqual(s0["clndn"], "Some_Disease")
        self.assertEqual(s0["cosmic_id"], "COSV99")
        self.assertEqual(s0["cosmic_gene"], "TP53")
        self.assertEqual(s0["svtype"], "DEL")
        self.assertEqual(s0["svlen"], -200)
        self.assertEqual(s0["end"], 300)
        self.assertEqual(s0["cipos"], (-50, 50))

        # Site without overlays must not have those keys present —
        # bcf_writer's _format_info uses .get() so absence == omitted
        # INFO fields, which is the contract we want.
        for k in (
            "clnsig", "clndn", "cosmic_id", "cosmic_gene",
            "svtype", "svlen", "end", "cipos",
        ):
            self.assertNotIn(k, s1, f"{k} should be absent on site_without")

    def test_empty_sites_iter(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_slice,
            read_arrow_metadata,
        )

        n_written = write_arrow_file(self.path, "1", 5, iter([]))
        self.assertEqual(n_written, 0)
        self.assertEqual(list(read_arrow_slice(self.path, 0, 5)), [])
        meta = read_arrow_metadata(self.path)
        self.assertEqual(meta["num_rows"], 0)
        self.assertEqual(meta["num_record_batches"], 0)


@unittest.skipUnless(HAS_PYARROW, "pyarrow not installed")
class TestBatchBoundaries(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "cohort.arrow"

    def tearDown(self):
        self.tmp.cleanup()

    def test_partial_final_batch(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_metadata,
            read_arrow_slice,
        )

        sites = [_site_from_gts(100 + i * 10, ["0|0", "0|1"]) for i in range(7)]
        write_arrow_file(self.path, "20", 2, iter(sites), batch_size=3)

        meta = read_arrow_metadata(self.path)
        self.assertEqual(meta["num_record_batches"], 3)  # 3+3+1
        self.assertEqual(meta["num_rows"], 7)

        read_back = list(read_arrow_slice(self.path, 0, 2))
        self.assertEqual(len(read_back), 7)
        self.assertEqual([s["pos"] for s in read_back],
                          [100, 110, 120, 130, 140, 150, 160])

    def test_exact_multiple_batches(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_metadata,
        )

        sites = [_site_from_gts(100 + i * 10, ["0|0", "0|1"]) for i in range(6)]
        write_arrow_file(self.path, "20", 2, iter(sites), batch_size=3)

        meta = read_arrow_metadata(self.path)
        self.assertEqual(meta["num_record_batches"], 2)  # 3+3
        self.assertEqual(meta["num_rows"], 6)

    def test_invalid_batch_size_raises(self):
        from syntheticgen.cohort_arrow import write_arrow_file

        with self.assertRaises(ValueError):
            write_arrow_file(
                self.path, "20", 2, iter([]), batch_size=0,
            )


@unittest.skipUnless(HAS_PYARROW, "pyarrow not installed")
class TestSliceBounds(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "cohort.arrow"

    def tearDown(self):
        self.tmp.cleanup()

    def test_empty_slice_yields_nothing(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_slice,
        )

        sites = [_site_from_gts(100, ["0|0", "0|1"])]
        write_arrow_file(self.path, "20", 2, iter(sites))
        self.assertEqual(list(read_arrow_slice(self.path, 1, 1)), [])

    def test_sample_hi_too_large_raises(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_slice,
        )

        sites = [_site_from_gts(100, ["0|0", "0|1"])]
        write_arrow_file(self.path, "20", 2, iter(sites))
        with self.assertRaises(ValueError):
            list(read_arrow_slice(self.path, 0, 5))

    def test_sample_lo_negative_raises(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_slice,
        )

        sites = [_site_from_gts(100, ["0|0", "0|1"])]
        write_arrow_file(self.path, "20", 2, iter(sites))
        with self.assertRaises(ValueError):
            list(read_arrow_slice(self.path, -1, 2))


@unittest.skipUnless(HAS_PYARROW, "pyarrow not installed")
class TestDiagnosticCarriers(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "cohort.arrow"

    def tearDown(self):
        self.tmp.cleanup()

    def test_round_trip_carriers(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_carriers,
        )

        # 3 people, 6 haplotypes. Set carriers explicitly.
        original_carriers = [(1, 1), (3, 1), (5, 1)]
        site = {
            "chrom": "20", "pos": 12345, "id": ".",
            "ref": "A", "alts": ["C"], "acs": [3], "afs": [0.5],
            "n_haplotypes": 6, "carriers": original_carriers,
        }
        write_arrow_file(self.path, "20", 3, iter([site]))

        recovered = read_arrow_carriers(self.path, pos=12345)
        # ``recovered`` is a 2D np.int32 array; sort both sides as
        # tuple-of-tuples so the comparison is order-independent and
        # hashable.
        self.assertEqual(
            sorted(tuple(row) for row in recovered.tolist()),
            sorted(original_carriers),
        )

    def test_unknown_pos_returns_empty(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_carriers,
        )

        site = _site_from_gts(100, ["0|0", "0|1"])
        write_arrow_file(self.path, "20", 2, iter([site]))
        # Missing position → empty packed array, shape (0, 2).
        result = read_arrow_carriers(self.path, pos=999)
        self.assertEqual(len(result), 0)
        self.assertEqual(result.shape, (0, 2))

    def test_multi_allelic_carriers(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_carriers,
        )

        # Multi-allelic with allele_idx = 2
        original_carriers = [(0, 1), (3, 2), (4, 1)]
        site = {
            "chrom": "20", "pos": 555, "id": ".",
            "ref": "A", "alts": ["C", "G"], "acs": [2, 1], "afs": [0.33, 0.17],
            "n_haplotypes": 6, "carriers": original_carriers,
        }
        write_arrow_file(self.path, "20", 3, iter([site]))
        recovered = read_arrow_carriers(self.path, pos=555)
        # Multi-allelic round-trip: the (3, 2) and (4, 1) carriers
        # must survive. ``recovered`` is a 2D np.int32 array.
        self.assertEqual(
            sorted(tuple(row) for row in recovered.tolist()),
            sorted(original_carriers),
        )


@unittest.skipUnless(HAS_PYARROW, "pyarrow not installed")
class TestMetadata(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "cohort.arrow"

    def tearDown(self):
        self.tmp.cleanup()

    def test_round_trip_metadata(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_metadata,
        )

        sites = [_site_from_gts(100 + i, ["0|0", "0|1", "1|1"]) for i in range(5)]
        write_arrow_file(self.path, "22", 3, iter(sites), batch_size=2)

        meta = read_arrow_metadata(self.path)
        self.assertEqual(meta["chrom"], "22")
        self.assertEqual(meta["n_samples"], 3)
        self.assertEqual(meta["n_haplotypes"], 6)
        self.assertEqual(meta["format_version"], "1")
        self.assertEqual(meta["num_rows"], 5)
        # 5 sites, batch_size=2 -> 3 batches (2+2+1)
        self.assertEqual(meta["num_record_batches"], 3)


@unittest.skipUnless(HAS_PYARROW, "pyarrow not installed")
class TestGtsFallback(unittest.TestCase):
    """Sites can be written via the legacy ``gts`` field instead of
    ``carriers`` — mirrors ``CohortBcfWriter.write_site`` dual support.
    """
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "cohort.arrow"

    def tearDown(self):
        self.tmp.cleanup()

    def test_writes_from_gts_field(self):
        from syntheticgen.cohort_arrow import (
            write_arrow_file,
            read_arrow_slice,
        )

        site = {
            "chrom": "20", "pos": 100, "id": ".",
            "ref": "A", "alts": ["C"], "acs": [2], "afs": [0.33],
            "n_haplotypes": 6,
            "gts": ["0|0", "0|1", "1|1"],
        }
        write_arrow_file(self.path, "20", 3, iter([site]))
        read_back = list(read_arrow_slice(self.path, 0, 3))
        self.assertEqual(read_back[0]["gts"], ["0|0", "0|1", "1|1"])


if __name__ == "__main__":
    unittest.main()
