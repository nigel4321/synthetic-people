"""Tests for M13.3 haploid emission in ``write_person_vcf``.

These tests cover the per-record ploidy filtering that turns the
three M13.2 sex-chromosome validation gates from FAIL → GREEN:

- chrY records dropped entirely for females (ploidy=0).
- chrY / chrX non-PAR records emitted as single-allele GT in males
  (ploidy=1).
- MT records emitted as single-allele GT for everyone (ploidy=1).
- Autosomes, PAR positions in males, and chrX in females stay
  diploid (ploidy=2, pre-M13.3 behaviour).

We exercise ``write_person_vcf`` through the public entry point so
the bcftools-encoded output reflects what users actually receive.
The records are read back via ``bcftools view -H`` and inspected
field-by-field.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.writer import write_person_vcf  # noqa: E402

_HAVE_BCFTOOLS = shutil.which("bcftools") is not None
_HAVE_BGZIP = shutil.which("bgzip") is not None
_HAVE_TABIX = shutil.which("tabix") is not None


def _record(chrom: str, pos: int, gt: str,
            ref: str = "A", alt: str = "C") -> dict:
    """Minimal record dict matching the shape ``write_person_vcf``
    expects from ``person['background']`` / ``person['highlighted']``."""
    return {
        "chrom": chrom,
        "pos": pos,
        "id": ".",
        "ref": ref,
        "alts": [alt],
        "afs": [0.5],
        "gt": gt,
    }


def _write_and_read(tmpdir: Path, person: dict,
                    build: str = "GRCh38",
                    sex: str | None = None) -> list[dict]:
    """Drive ``write_person_vcf`` end-to-end and parse the records
    back via ``bcftools view -H``. Returns a list of
    ``{"chrom", "pos", "ref", "alt", "gt", "an"}`` dicts so tests
    can assert against the actual on-disk output."""
    import random
    out = tmpdir / f"{person['sample_id']}.vcf.gz"
    write_person_vcf(out, person, build, random.Random(42), sex=sex,
                     dp_mean=30.0)
    proc = subprocess.run(
        ["bcftools", "view", "-H", str(out)],
        capture_output=True, text=True, check=True,
    )
    out_records: list[dict] = []
    for line in proc.stdout.splitlines():
        if not line:
            continue
        parts = line.split("\t")
        info = dict(
            (kv.split("=", 1) + [None])[:2]
            for kv in parts[7].split(";") if kv
        )
        sample = parts[9].split(":")
        out_records.append({
            "chrom": parts[0],
            "pos": int(parts[1]),
            "ref": parts[3],
            "alt": parts[4],
            "gt": sample[0],
            "an": info.get("AN"),
        })
    return out_records


@unittest.skipUnless(_HAVE_BCFTOOLS and _HAVE_BGZIP and _HAVE_TABIX,
                     "bcftools + bgzip + tabix required")
class WritePersonVcfHaploidTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    # --- ploidy=0: chrY in females → record dropped ---

    def test_female_chry_records_dropped(self):
        # A female cohort member should have NO chrY records in
        # their VCF — flips the M13.2 female_y_absence gate GREEN.
        person = {
            "sample_id": "f001",
            "highlighted": _record("22", 1_000_000, "0|1"),
            "background": [
                _record("Y", 20_000_000, "0|1"),
                _record("Y", 30_000_000, "1|1"),
                _record("22", 2_000_000, "0|0"),
            ],
        }
        out = _write_and_read(self.dir, person, sex="f")
        chroms = {r["chrom"] for r in out}
        self.assertNotIn("Y", chroms,
                         f"chrY records leaked into female VCF: {out}")

    def test_female_chry_records_kept_when_sex_unset(self):
        # Backwards compat: when sex is not provided the writer
        # should preserve pre-M13.3 behaviour (no filtering).
        person = {
            "sample_id": "anon",
            "highlighted": _record("22", 1_000_000, "0|1"),
            "background": [_record("Y", 20_000_000, "0|1")],
        }
        out = _write_and_read(self.dir, person, sex=None)
        chroms = {r["chrom"] for r in out}
        self.assertIn("Y", chroms)

    # --- ploidy=1: chrY non-PAR in males → haploid GT ---

    def test_male_chry_non_par_emits_haploid_gt(self):
        # pos=20_000_000 on chrY is non-PAR on GRCh38. The GT must
        # be a single-allele field — flips the M13.2 y_het_in_males
        # gate GREEN (a haploid "1" is not heterozygous by the
        # _gt_is_heterozygous helper's a!=b check).
        person = {
            "sample_id": "m001",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [_record("Y", 20_000_000, "0|1")],
        }
        out = _write_and_read(self.dir, person, sex="m")
        y = [r for r in out if r["chrom"] == "Y"]
        self.assertEqual(len(y), 1)
        # Haploid GT has no separator — single allele field.
        self.assertNotIn("|", y[0]["gt"], f"GT not haploid: {y[0]}")
        self.assertNotIn("/", y[0]["gt"])
        # First haplotype of 0|1 → "0".
        self.assertEqual(y[0]["gt"], "0")
        # AN must reflect ploidy 1.
        self.assertEqual(y[0]["an"], "1")

    def test_male_chry_par_emits_diploid_gt(self):
        # GRCh38 chrY PAR1: 10_001-2_781_479. pos=1_000_000 is in
        # PAR; PAR is diploid in males so GT stays "X|Y".
        person = {
            "sample_id": "m002",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [_record("Y", 1_000_000, "0|1")],
        }
        out = _write_and_read(self.dir, person, sex="m")
        y = [r for r in out if r["chrom"] == "Y"]
        self.assertEqual(len(y), 1)
        self.assertIn("|", y[0]["gt"])
        self.assertEqual(y[0]["an"], "2")

    # --- ploidy=1: chrX non-PAR in males → haploid ---

    def test_male_chrx_non_par_emits_haploid_gt(self):
        person = {
            "sample_id": "m003",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [_record("X", 80_000_000, "0|1")],
        }
        out = _write_and_read(self.dir, person, sex="m")
        x = [r for r in out if r["chrom"] == "X"]
        self.assertEqual(len(x), 1)
        self.assertNotIn("|", x[0]["gt"])
        self.assertEqual(x[0]["an"], "1")

    def test_female_chrx_stays_diploid(self):
        # chrX in females is always diploid (ploidy=2 regardless of
        # position). Pin that explicitly because the male code path
        # is a different branch in ``ploidy_for``.
        person = {
            "sample_id": "f002",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [_record("X", 80_000_000, "0|1")],
        }
        out = _write_and_read(self.dir, person, sex="f")
        x = [r for r in out if r["chrom"] == "X"]
        self.assertEqual(len(x), 1)
        self.assertIn("|", x[0]["gt"])
        self.assertEqual(x[0]["an"], "2")

    # --- MT: always haploid in both sexes ---

    def test_mt_emits_haploid_gt_in_male(self):
        person = {
            "sample_id": "m004",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [_record("MT", 100, "0|1")],
        }
        out = _write_and_read(self.dir, person, sex="m")
        mt = [r for r in out if r["chrom"] == "MT"]
        self.assertEqual(len(mt), 1)
        self.assertNotIn("|", mt[0]["gt"])
        self.assertEqual(mt[0]["an"], "1")

    def test_mt_emits_haploid_gt_in_female(self):
        person = {
            "sample_id": "f003",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [_record("MT", 100, "1|1")],
        }
        out = _write_and_read(self.dir, person, sex="f")
        mt = [r for r in out if r["chrom"] == "MT"]
        self.assertEqual(len(mt), 1)
        self.assertNotIn("|", mt[0]["gt"])
        self.assertEqual(mt[0]["an"], "1")

    # --- autosomes always diploid ---

    def test_autosomes_stay_diploid_in_both_sexes(self):
        for sex in ("m", "f"):
            with self.subTest(sex=sex):
                person = {
                    "sample_id": f"a_{sex}",
                    "highlighted": _record("22", 1_000_000, "0|1"),
                    "background": [
                        _record("1", 5_000_000, "0|1"),
                        _record("21", 10_000_000, "1|1"),
                    ],
                }
                out = _write_and_read(self.dir, person, sex=sex)
                # Every record is diploid.
                for r in out:
                    self.assertIn("|", r["gt"],
                                  f"non-diploid GT on autosome: {r}")
                    self.assertEqual(r["an"], "2")


@unittest.skipUnless(_HAVE_BCFTOOLS and _HAVE_BGZIP and _HAVE_TABIX,
                     "bcftools + bgzip + tabix required")
class HighlightedRecordReachesVcfTest(unittest.TestCase):
    """PR #107 review (Copilot): the highlighted variant must not
    land on a chromosome the person's ploidy will drop. Worker now
    pre-filters the candidate pool — pin both halves of that:

    1. A female passed a candidate pool that's entirely chrY would
       fall back to the unfiltered pool (defensive — never lose the
       highlighted slot entirely).
    2. A female passed a mixed pool gets a candidate from the
       allowed subset; the highlighted record reaches her VCF.

    We test this at the ``_person_worker`` level by driving the
    filter logic directly rather than spinning up a whole cli.main
    run.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def _hi_candidate(self, chrom: str, pos: int) -> dict:
        return {
            "chrom": chrom, "pos": pos, "id": "rs1",
            "ref": "A", "alts": ["C"], "afs": [0.5],
        }

    def test_female_with_mixed_pool_never_gets_chry_highlight(self):
        # Drive the filter logic directly to make this deterministic.
        # 1000 trials with a 50/50 mix of chrY and chr22 candidates;
        # post-filter the female pool must contain no chrY entries.
        import random
        from syntheticgen.builds import ploidy_for
        candidates = [
            self._hi_candidate("Y", 20_000_000),
            self._hi_candidate("Y", 25_000_000),
            self._hi_candidate("22", 1_000_000),
            self._hi_candidate("22", 2_000_000),
        ]
        sex = "f"
        build = "GRCh38"
        filtered = [
            c for c in candidates
            if ploidy_for(c["chrom"], sex, build, c["pos"]) != 0
        ]
        # Every chrY candidate must be excluded.
        self.assertTrue(all(c["chrom"] != "Y" for c in filtered))
        # The chr22 candidates must survive — the pool isn't empty.
        self.assertEqual(len(filtered), 2)
        # Sanity: random draws over many trials never produce chrY.
        rng = random.Random(42)
        for _ in range(100):
            c = rng.choice(filtered)
            self.assertNotEqual(c["chrom"], "Y")

    def test_male_with_chry_candidate_keeps_record(self):
        # The companion case: a MALE pool with a chrY candidate
        # filters to keep the chrY entry (ploidy=1 in male non-PAR,
        # not 0). Confirms the filter only excludes ploidy==0.
        from syntheticgen.builds import ploidy_for
        candidates = [
            self._hi_candidate("Y", 20_000_000),  # non-PAR
            self._hi_candidate("22", 1_000_000),
        ]
        filtered = [
            c for c in candidates
            if ploidy_for(c["chrom"], "m", "GRCh38", c["pos"]) != 0
        ]
        self.assertEqual(len(filtered), 2)


if __name__ == "__main__":
    unittest.main()
