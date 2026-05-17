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

    def test_male_chry_par_input_records_are_dropped(self):
        # M13.4 changes the semantics: chrY PAR records from the
        # simulator are now dropped (they'd diverge from the
        # biologically-identical chrX PAR variants). Without a
        # chrX PAR input to mirror from, the male VCF has no chrY
        # PAR records at all.
        #
        # PAR diploid emission is now exercised via chrX PAR
        # records being mirrored onto chrY — see
        # ``test_male_chrx_par1_mirrored_onto_chry`` in
        # ``WritePersonVcfParCopyTest`` below.
        person = {
            "sample_id": "m002",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [_record("Y", 1_000_000, "0|1")],
        }
        out = _write_and_read(self.dir, person, sex="m")
        y = [r for r in out if r["chrom"] == "Y"]
        self.assertEqual(
            len(y), 0,
            "chrY PAR record from simulator must be dropped (M13.4)",
        )

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
class WritePersonVcfParCopyTest(unittest.TestCase):
    """M13.4: PAR1 / PAR2 copy mechanism. chrX PAR variants are
    mirrored onto chrY at the build-correct translated position
    for males; chrY PAR records from the simulator are dropped;
    females are unchanged (chrY entirely absent).

    These tests exercise the post-record-assembly pre-emission
    transformation in ``write_person_vcf`` end-to-end (write →
    ``bcftools view -H`` → assert chrom/pos pairs)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_male_chrx_par1_mirrored_onto_chry(self):
        # GRCh38 PAR1: chrX 10_001-2_781_479 ↔ chrY 10_001-2_781_479
        # (identical bp on this build). A chrX PAR1 variant must
        # appear on chrY at the same position in males.
        person = {
            "sample_id": "m_par1",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [_record("X", 1_000_000, "0|1")],
        }
        out = _write_and_read(self.dir, person, sex="m")
        x = [r for r in out if r["chrom"] == "X"]
        y = [r for r in out if r["chrom"] == "Y"]
        self.assertEqual(len(x), 1, "chrX original must survive")
        self.assertEqual(len(y), 1, "chrY mirror must be emitted")
        # Same position (PAR1 maps 1:1 on GRCh38).
        self.assertEqual(y[0]["pos"], 1_000_000)
        # Same GT, diploid, AN=2 (PAR in male is diploid).
        self.assertEqual(y[0]["gt"], x[0]["gt"])
        self.assertEqual(y[0]["an"], "2")

    def test_male_chrx_par2_mirrored_with_offset(self):
        # GRCh38 PAR2: chrX 155_701_383-156_030_895 ↔ chrY
        # 56_887_903-57_217_415. The lengths match but the start
        # positions don't, so the mirror needs an offset translation.
        # pos=155_800_000 on chrX → expected chrY pos:
        #   56_887_903 + (155_800_000 - 155_701_383) = 56_986_520
        person = {
            "sample_id": "m_par2",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [_record("X", 155_800_000, "0|1")],
        }
        out = _write_and_read(self.dir, person, sex="m")
        y = [r for r in out if r["chrom"] == "Y"]
        self.assertEqual(len(y), 1)
        self.assertEqual(y[0]["pos"], 56_986_520)

    def test_male_chrx_non_par_NOT_mirrored_onto_chry(self):
        # chrX non-PAR has no counterpart on chrY — it's a male-
        # specific X region. No chrY mirror should be emitted.
        person = {
            "sample_id": "m_nonpar",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [_record("X", 80_000_000, "0|1")],
        }
        out = _write_and_read(self.dir, person, sex="m")
        y = [r for r in out if r["chrom"] == "Y"]
        self.assertEqual(len(y), 0, f"chrY mirror leaked from non-PAR: {y}")

    def test_female_chrx_par_NOT_mirrored_onto_chry(self):
        # Females have no chrY at all (M13.3 drops it). The PAR
        # mirror must NOT reintroduce chrY for them.
        person = {
            "sample_id": "f_par",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [_record("X", 1_000_000, "0|1")],
        }
        out = _write_and_read(self.dir, person, sex="f")
        chroms = {r["chrom"] for r in out}
        self.assertNotIn("Y", chroms,
                         f"chrY mirror leaked into female VCF: {out}")

    def test_male_chry_par_from_simulator_dropped(self):
        # A chrY PAR record from independent simulator output must
        # be dropped, because PAR is supposed to be IDENTICAL to
        # chrX PAR — the independent simulator value would diverge.
        # Without any chrX PAR record to mirror from, the male VCF
        # has no chrY PAR records.
        person = {
            "sample_id": "m_y_par_only",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [_record("Y", 1_000_000, "0|1")],
        }
        out = _write_and_read(self.dir, person, sex="m")
        y = [r for r in out if r["chrom"] == "Y"]
        self.assertEqual(len(y), 0)

    def test_male_chry_par_replaced_by_chrx_par_mirror(self):
        # Both inputs present: chrX PAR + chrY PAR at the same
        # position. The chrY input is dropped; the chrX value is
        # mirrored. Net effect: chrY records carry the chrX GT,
        # not the original chrY GT.
        person = {
            "sample_id": "m_both_par",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [
                _record("X", 1_000_000, "0|1"),  # X PAR1 variant
                _record("Y", 1_000_000, "1|1"),  # different GT
            ],
        }
        out = _write_and_read(self.dir, person, sex="m")
        y = [r for r in out if r["chrom"] == "Y"]
        self.assertEqual(len(y), 1)
        # GT comes from the chrX record (0|1), not the dropped chrY
        # input (which was 1|1).
        self.assertEqual(y[0]["gt"], "0|1")

    def test_male_chry_non_par_preserved(self):
        # Non-PAR chrY records (which only exist on chrY, no chrX
        # counterpart) are NOT dropped by M13.4 — M13.3's haploid
        # emission handles them. Pin that here so a future refactor
        # doesn't over-aggressively drop non-PAR Y too.
        person = {
            "sample_id": "m_y_nonpar",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [_record("Y", 20_000_000, "0|1")],  # non-PAR
        }
        out = _write_and_read(self.dir, person, sex="m")
        y = [r for r in out if r["chrom"] == "Y"]
        self.assertEqual(len(y), 1)
        # Haploid (M13.3) — not "|"-separated.
        self.assertNotIn("|", y[0]["gt"])

    def test_mirror_does_not_carry_highlight_info(self):
        # PR #110 review (Copilot): the chrY PAR mirror must NEVER
        # carry the HIGHLIGHT INFO flag — only the chrX original
        # does — otherwise a single biological variant emits TWO
        # highlighted rows. Pin that with a dedicated INFO check.
        person = {
            "sample_id": "m_hi_par",
            # Highlighted = chrX PAR1 variant. The chrX row should
            # carry HIGHLIGHT; the mirrored chrY row must NOT.
            "highlighted": _record("X", 1_500_000, "0|1"),
            "background": [],
        }
        # We need the raw VCF text, not the parsed dict, to inspect
        # INFO contents — re-run bcftools view -H with the full
        # INFO column visible.
        import random
        out_path = self.dir / f"{person['sample_id']}.vcf.gz"
        write_person_vcf(out_path, person, "GRCh38",
                         random.Random(42), sex="m", dp_mean=30.0)
        proc = subprocess.run(
            ["bcftools", "view", "-H", str(out_path)],
            capture_output=True, text=True, check=True,
        )
        x_lines = [
            ln for ln in proc.stdout.splitlines()
            if ln.split("\t")[0] == "X"
        ]
        y_lines = [
            ln for ln in proc.stdout.splitlines()
            if ln.split("\t")[0] == "Y"
        ]
        self.assertEqual(len(x_lines), 1)
        self.assertEqual(len(y_lines), 1)
        self.assertIn("HIGHLIGHT", x_lines[0].split("\t")[7])
        self.assertNotIn("HIGHLIGHT", y_lines[0].split("\t")[7])

    def test_sv_records_not_mirrored(self):
        # PR #110 review (suppressed): SVs carry coordinate-bearing
        # INFO (END, SVLEN) that don't survive a chrom/pos swap
        # without translation, and the SV span can extend past the
        # PAR boundary. Skip mirroring for SVs to avoid emitting a
        # nonsensical chrY record with chrX-coordinate END.
        sv_variant = _record("X", 1_000_000, "0|1")
        sv_variant.update({
            "svtype": "DEL",
            "svlen": -500,
            "end": 1_000_500,
            "cipos": (-50, 50),
            "alts": ["<DEL>"],
        })
        person = {
            "sample_id": "m_sv",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [sv_variant],
        }
        out = _write_and_read(self.dir, person, sex="m")
        y = [r for r in out if r["chrom"] == "Y"]
        self.assertEqual(
            len(y), 0,
            f"SV record was mirrored onto chrY (END would be at "
            f"chrX coordinates): {y}",
        )

    def test_no_par_copy_when_sex_unset(self):
        # Back-compat: sex=None should preserve pre-M13.4 behaviour
        # — chrY records pass through, no chrX-to-chrY mirroring.
        person = {
            "sample_id": "anon",
            "highlighted": _record("22", 1_000_000, "0|0"),
            "background": [
                _record("X", 1_000_000, "0|1"),  # chrX PAR — would mirror in male
                _record("Y", 1_000_000, "1|1"),  # chrY PAR — would be dropped in any-sex
            ],
        }
        out = _write_and_read(self.dir, person, sex=None)
        y = [r for r in out if r["chrom"] == "Y"]
        # The original chrY record survives — no PAR drop.
        self.assertEqual(len(y), 1)
        self.assertEqual(y[0]["gt"], "1|1")


class HighlightedRecordReachesVcfTest(unittest.TestCase):
    """PR #107 + PR #108 review (Copilot): the highlighted variant
    must not land on a chromosome the person's ploidy will drop.
    Worker pre-filters the candidate pool by
    ``ploidy_for(...) != 0``; if the filter empties the pool, the
    worker fails fast with a clear error (PR #108 review — the
    earlier silent fallback to the unfiltered pool reintroduced
    the very bug the filter was added to prevent).

    Three behaviours pinned below:

    1. Mixed pool + female → chrY candidates excluded, autosomes
       survive, draws never land on chrY.
    2. Mixed pool + male → every candidate kept (chrY non-PAR is
       ploidy=1 for males).
    3. All-chrY pool + female → ``_person_worker`` raises
       ``RuntimeError`` with a clear error message.

    No bcftools subprocesses involved — these tests exercise only
    the filter logic + the fail-fast path. No ``skipUnless``
    decorator needed (PR #108 review: don't unnecessarily skip
    tests that don't shell out).
    """

    def _hi_candidate(self, chrom: str, pos: int) -> dict:
        return {
            "chrom": chrom, "pos": pos, "id": "rs1",
            "ref": "A", "alts": ["C"], "afs": [0.5],
        }

    def test_female_with_mixed_pool_never_gets_chry_highlight(self):
        # Drive the filter logic directly to make this deterministic.
        # 100 random draws over a 50/50 chrY + chr22 candidate set;
        # post-filter the female pool must contain no chrY entries
        # and the draws never land on chrY.
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
        self.assertTrue(all(c["chrom"] != "Y" for c in filtered))
        self.assertEqual(len(filtered), 2)
        rng = random.Random(42)
        for _ in range(100):
            c = rng.choice(filtered)
            self.assertNotEqual(c["chrom"], "Y")

    def test_male_with_chry_candidate_keeps_record(self):
        # The companion case: a MALE pool with a chrY NON-PAR
        # candidate filters to keep the chrY entry (ploidy=1 in
        # male non-PAR, not 0). Confirms the filter only excludes
        # ploidy==0 — and M13.4's additional chrY-PAR exclusion
        # (next test) doesn't fire for non-PAR positions.
        from syntheticgen.builds import is_in_par, ploidy_for
        candidates = [
            self._hi_candidate("Y", 20_000_000),  # non-PAR
            self._hi_candidate("22", 1_000_000),
        ]
        sex = "m"
        build = "GRCh38"
        filtered = [
            c for c in candidates
            if ploidy_for(c["chrom"], sex, build, c["pos"]) != 0
            and not (c["chrom"] == "Y"
                     and is_in_par("Y", c["pos"], build))
        ]
        self.assertEqual(len(filtered), 2)

    def test_male_chry_par_candidate_filtered_out(self):
        # PR #110 review (Copilot): a male could draw a chrY PAR
        # ClinVar candidate that has ploidy=2 in males (PAR is
        # diploid), so the M13.3 ploidy filter alone doesn't catch
        # it. But M13.4 drops chrY PAR records — so the highlighted
        # variant disappears silently. The M13.4 filter excludes
        # chrY PAR candidates from the highlighted draw entirely.
        from syntheticgen.builds import is_in_par, ploidy_for
        candidates = [
            # Three chrY PAR candidates (all in PAR1 on GRCh38)
            self._hi_candidate("Y", 100_000),
            self._hi_candidate("Y", 1_500_000),
            self._hi_candidate("Y", 2_500_000),
            # One chr22 candidate as the only valid option
            self._hi_candidate("22", 1_000_000),
        ]
        sex = "m"
        build = "GRCh38"
        filtered = [
            c for c in candidates
            if ploidy_for(c["chrom"], sex, build, c["pos"]) != 0
            and not (c["chrom"] == "Y"
                     and is_in_par("Y", c["pos"], build))
        ]
        # All three chrY PAR candidates must be excluded; only
        # the chr22 candidate survives.
        self.assertEqual(len(filtered), 1)
        self.assertEqual(filtered[0]["chrom"], "22")

    def test_empty_pool_for_female_raises_in_worker(self):
        # PR #108 review (Copilot): the empty-pool fallback used to
        # silently revert to the unfiltered pool, reintroducing the
        # exact bug the filter was added to prevent. Now it raises
        # ``RuntimeError`` with a clear actionable message.
        #
        # Drive ``_person_worker`` directly with a candidate pool
        # that's entirely chrY and a female sex assignment.
        from syntheticgen import cli as cli_module
        candidates = [
            self._hi_candidate("Y", 20_000_000),
            self._hi_candidate("Y", 25_000_000),
        ]
        # Minimum viable worker state: only the fields the candidate-
        # pool-filter path touches are populated. The worker will
        # raise before any of the unset fields get consulted.
        cli_module._PERSON_WORKER_STATE.clear()
        cli_module._PERSON_WORKER_STATE.update({
            "candidates": candidates,
            "cohort_sites": None,
            "cohort_bcfs": None,
            "build": "GRCh38",
            "output_dir": Path("/tmp"),
            "truth_dir": Path("/tmp"),
            "contig_order": {"Y": 23, "22": 21},
            "svs_per_person": 0,
            "sv_length_min": 50,
            "sv_length_max": 1000,
            "sv_chrom_span": 1000,
            "sv_chromosomes": ["Y"],
            "error_rate": 0.0,
            "dropout_rate": 0.0,
            "person_ancestry": [],
            "person_sexes": ["f"],
        })
        try:
            with self.assertRaises(RuntimeError) as ctx:
                cli_module._person_worker(0, "F001", seed=42)
            msg = str(ctx.exception)
            self.assertIn("highlighted-candidate pool empty", msg)
            self.assertIn("sex='f'", msg)
            self.assertIn("--clinvar-sig", msg)
        finally:
            cli_module._PERSON_WORKER_STATE.clear()


if __name__ == "__main__":
    unittest.main()
