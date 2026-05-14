"""Tests for the M12 reference-FASTA helpers.

The module's contract:
- Resolve the FASTA-side chromosome name regardless of cli/UCSC
  chr-prefix convention (``"22"`` vs ``"chr22"``).
- Return the uppercase REF base at a 1-based position.
- Fall back to ``"N"`` on missing chrom / out-of-range pos rather
  than raising per-variant (which would tank a multi-hour run).
- Pre-flight ``validate_fasta`` surfaces the typical "wrong FASTA"
  case at startup with a clear actionable message.

Tests gate on pysam being importable so the rest of the suite still
runs in stripped environments.
"""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


_HAVE_PYSAM = importlib.util.find_spec("pysam") is not None


def _write_fasta(dir_: Path, records: dict) -> Path:
    """Build a tiny FASTA + .fai index from ``{name: sequence}``."""
    import pysam
    dir_.mkdir(parents=True, exist_ok=True)
    path = dir_ / "tiny.fa"
    path.write_text(
        "\n".join(
            f">{name}\n{seq}" for name, seq in records.items()
        ) + "\n",
    )
    pysam.faidx(str(path))
    return path


@unittest.skipUnless(_HAVE_PYSAM, "pysam not installed")
class ReferenceFastaTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_have_pysam_true_when_installed(self):
        from syntheticgen.reference import have_pysam
        self.assertTrue(have_pysam())

    def test_load_fasta_round_trip(self):
        from syntheticgen.reference import load_fasta, fetch_ref_base
        path = _write_fasta(self.dir, {"22": "ACGTACGT"})
        fa = load_fasta(path)
        self.assertEqual(fetch_ref_base(fa, "22", 1), "A")
        self.assertEqual(fetch_ref_base(fa, "22", 2), "C")
        self.assertEqual(fetch_ref_base(fa, "22", 8), "T")

    def test_load_fasta_raises_on_missing_file(self):
        from syntheticgen.reference import load_fasta
        with self.assertRaises(FileNotFoundError):
            load_fasta(Path("/nonexistent.fa"))

    def test_fetch_ref_handles_chr_prefix_either_way(self):
        # FASTA uses unprefixed names; caller asks for "chr22".
        from syntheticgen.reference import load_fasta, fetch_ref_base
        path = _write_fasta(self.dir, {"22": "ACGT"})
        fa = load_fasta(path)
        self.assertEqual(fetch_ref_base(fa, "chr22", 1), "A")
        # And the reverse: FASTA uses "chr22", caller asks for "22".
        path2 = _write_fasta(
            self.dir / "p2", {"chr22": "TTTT"},
        )
        fa2 = load_fasta(path2)
        self.assertEqual(fetch_ref_base(fa2, "22", 1), "T")

    def test_fetch_ref_returns_N_on_missing_chrom(self):
        # Pre-flight validate_fasta should catch this at startup,
        # but the per-variant path needs to degrade gracefully too.
        from syntheticgen.reference import load_fasta, fetch_ref_base
        path = _write_fasta(self.dir, {"22": "ACGT"})
        fa = load_fasta(path)
        self.assertEqual(fetch_ref_base(fa, "99", 1), "N")

    def test_fetch_ref_returns_N_on_out_of_range_pos(self):
        from syntheticgen.reference import load_fasta, fetch_ref_base
        path = _write_fasta(self.dir, {"22": "ACGT"})
        fa = load_fasta(path)
        self.assertEqual(fetch_ref_base(fa, "22", 999), "N")

    def test_fetch_ref_normalises_lowercase_bases(self):
        # Soft-masked FASTAs use lowercase. Genome-wide consumers
        # always want uppercase.
        from syntheticgen.reference import load_fasta, fetch_ref_base
        path = _write_fasta(self.dir, {"22": "acgt"})
        fa = load_fasta(path)
        self.assertEqual(fetch_ref_base(fa, "22", 1), "A")

    def test_fetch_ref_none_fasta_returns_N(self):
        from syntheticgen.reference import fetch_ref_base
        self.assertEqual(fetch_ref_base(None, "22", 1), "N")

    def test_resolve_chrom_name_round_trip(self):
        from syntheticgen.reference import load_fasta, resolve_chrom_name
        path = _write_fasta(self.dir, {"22": "ACGT", "chrX": "TTTT"})
        fa = load_fasta(path)
        self.assertEqual(resolve_chrom_name(fa, "22"), "22")
        self.assertEqual(resolve_chrom_name(fa, "chr22"), "22")
        self.assertEqual(resolve_chrom_name(fa, "X"), "chrX")
        self.assertEqual(resolve_chrom_name(fa, "chrX"), "chrX")
        self.assertIsNone(resolve_chrom_name(fa, "99"))

    def test_validate_fasta_passes_when_chroms_present(self):
        from syntheticgen.reference import load_fasta, validate_fasta
        path = _write_fasta(self.dir, {"22": "A" * 1000})
        fa = load_fasta(path)
        # Should not raise.
        validate_fasta(fa, ["22"], chr_length_mb=0.0)

    def test_validate_fasta_raises_on_missing_chrom(self):
        from syntheticgen.reference import load_fasta, validate_fasta
        path = _write_fasta(self.dir, {"22": "A" * 100})
        fa = load_fasta(path)
        with self.assertRaises(ValueError) as ctx:
            validate_fasta(fa, ["22", "99"], chr_length_mb=0.0)
        self.assertIn("99", str(ctx.exception))
        self.assertIn("missing", str(ctx.exception).lower())

    def test_validate_fasta_raises_when_chrom_too_short(self):
        from syntheticgen.reference import load_fasta, validate_fasta
        # 1 kb FASTA but the caller wants 5 Mb of chr22.
        path = _write_fasta(self.dir, {"22": "A" * 1000})
        fa = load_fasta(path)
        with self.assertRaises(ValueError) as ctx:
            validate_fasta(fa, ["22"], chr_length_mb=5.0)
        self.assertIn("shorter than", str(ctx.exception))

    def test_validate_fasta_skips_length_check_when_chr_length_zero(self):
        # ``--chr-length-mb 0`` means "use full contig"; we don't
        # have a target length to enforce so length validation
        # is skipped (presence validation still applies).
        from syntheticgen.reference import load_fasta, validate_fasta
        path = _write_fasta(self.dir, {"22": "A" * 100})
        fa = load_fasta(path)
        validate_fasta(fa, ["22"], chr_length_mb=0.0)

    def test_validate_fasta_with_none_is_noop(self):
        # No FASTA passed → no validation needed.
        from syntheticgen.reference import validate_fasta
        validate_fasta(None, ["22"], chr_length_mb=5.0)


class ReferenceFastaWithoutPysamTest(unittest.TestCase):
    """When pysam isn't installed, ``load_fasta`` must raise
    ``ImportError`` with an actionable cli hint. The actual import
    is harder to fake here than to mock — patch ``have_pysam``
    and verify the error path."""

    def test_load_fasta_raises_actionable_when_pysam_missing(self):
        from unittest.mock import patch
        with patch(
            "syntheticgen.reference.have_pysam",
            return_value=False,
        ):
            from syntheticgen.reference import load_fasta
            with self.assertRaises(ImportError) as ctx:
                load_fasta(Path("/nonexistent.fa"))
            self.assertIn("pysam", str(ctx.exception))
            self.assertIn("--reference-fasta", str(ctx.exception))


@unittest.skipUnless(_HAVE_PYSAM, "pysam not installed")
class PickRefTest(unittest.TestCase):
    """``_pick_ref`` in coalescent.py is the single bottleneck where
    the FASTA path and the rng-fabricated path coexist. The key
    invariant: rng state advances identically with or without FASTA,
    so downstream rng consumers (overlay sampling, error model) see
    the same state regardless of REF source."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_rng_state_identical_with_and_without_fasta(self):
        # Build a tiny FASTA with a known base at pos 1.
        from syntheticgen.reference import load_fasta
        from syntheticgen.coalescent import _pick_ref
        import random
        path = _write_fasta(self.dir, {"22": "GACGT"})
        fa = load_fasta(path)

        # Same seed, two parallel rngs. Both _pick_ref calls
        # should consume exactly one rng draw each, leaving the
        # rngs in identical states despite returning different
        # REF values.
        rng_with = random.Random(42)
        rng_without = random.Random(42)
        ref_with = _pick_ref(rng_with, fa, "22", 1)
        ref_without = _pick_ref(rng_without, None, "22", 1)
        # FASTA returns the real REF; rng path returns whatever
        # rng.choice("ACGT") picked.
        self.assertEqual(ref_with, "G")
        # rng-only path got the rng draw — could be any base.
        self.assertIn(ref_without, ("A", "C", "G", "T"))
        # CRITICAL: rng states match — both rngs consumed one
        # ``choice("ACGT")`` call. Any downstream rng draw will
        # produce the same value from both.
        self.assertEqual(rng_with.random(), rng_without.random())

    def test_fasta_N_falls_back_to_rng_draw(self):
        # When FASTA returns N (missing chrom / oob), fall back
        # to the rng-drawn base. The rng was already consumed
        # for the fallback draw, so this is symmetric with the
        # no-FASTA path.
        from syntheticgen.reference import load_fasta
        from syntheticgen.coalescent import _pick_ref
        import random
        path = _write_fasta(self.dir, {"22": "ACGT"})
        fa = load_fasta(path)
        rng = random.Random(42)
        ref = _pick_ref(rng, fa, "99", 1)  # chrom doesn't exist
        # Must be one of the canonical bases (from the fallback rng
        # draw), not "N".
        self.assertIn(ref, ("A", "C", "G", "T"))


_HAVE_BCFTOOLS = importlib.util.find_spec("subprocess") is not None
import shutil as _shutil
_BCFTOOLS_BIN = _shutil.which("bcftools")
_HAVE_MSPRIME = importlib.util.find_spec("msprime") is not None
_HAVE_STDPOPSIM = importlib.util.find_spec("stdpopsim") is not None


@unittest.skipUnless(
    _HAVE_PYSAM and _BCFTOOLS_BIN and _HAVE_MSPRIME and _HAVE_STDPOPSIM,
    "needs pysam + bcftools + msprime + stdpopsim",
)
class ReferenceEndToEndTest(unittest.TestCase):
    """End-to-end: run the cli with ``--reference-fasta`` and
    verify the emitted BCF's REF column matches the FASTA at the
    corresponding POS.

    This is the empirical proof point Tier 1's REF-check gate was
    designed for — when this test passes, M12's wiring is
    confirmed end-to-end and the gate goes from "fails on every
    record by design" to "passes" on real runs.
    """

    def _build_fasta(self, tmp: Path) -> Path:
        # 1 Mb of deterministic chr22 — repeating "ACGT" so every
        # POS has a known base. Real FASTAs are GRCh38; here we
        # need just enough that --chr-length-mb 0.5 is satisfied.
        # Use ACGT repeats; the cli's variants will land at random
        # positions and each emitted REF must match the local FASTA
        # lookup.
        import pysam
        path = tmp / "synthetic_chr22.fa"
        seq = ("ACGT" * 250_000)  # 1 Mb
        path.write_text(f">22\n{seq}\n")
        pysam.faidx(str(path))
        return path

    def _common_args(self, out_dir: Path, fasta: Path) -> list:
        return [
            "--no-config",
            "--n", "3",
            "--seed", "42",
            "--build", "GRCh38",
            "--chromosomes", "22",
            "--chr-length-mb", "0.5",
            "--demo-model", "none",
            "--rsid-density", "0",
            "--clinvar-inject-density", "0",
            "--svs-per-person", "0",
            "--error-rate", "0",
            "--dropout-rate", "0",
            "--workers", "1",
            "--output-dir", str(out_dir),
            "--cache-dir", str(out_dir / "cache"),
            "--mode", "cohort",
            "--cohort-mode", "sites_list",
            "--reference-fasta", str(fasta),
        ]

    def test_emitted_ref_matches_fasta(self):
        import subprocess
        import pysam
        from syntheticgen import cli as cli_module
        with tempfile.TemporaryDirectory() as tmp_str:
            tmp = Path(tmp_str)
            fasta_path = self._build_fasta(tmp)
            out_dir = tmp / "out"
            rc = cli_module.main(self._common_args(out_dir, fasta_path))
            self.assertEqual(rc, 0)
            cohort_bcf = out_dir / "cohort" / "cohort.chr22.bcf"
            self.assertTrue(cohort_bcf.is_file(),
                            f"cohort BCF not written: {cohort_bcf}")
            # Pull (POS, REF) tuples from the emitted BCF.
            out = subprocess.check_output(
                ["bcftools", "query", "-f", "%POS\t%REF\n",
                 str(cohort_bcf)],
                text=True,
            )
            records = [
                tuple(line.split("\t")) for line in out.strip().splitlines()
                if line
            ]
            self.assertGreater(
                len(records), 0,
                "no variants emitted at the canary scale",
            )
            # Verify every emitted REF matches the FASTA base at POS.
            fa = pysam.FastaFile(str(fasta_path))
            mismatches = []
            for pos_s, ref in records:
                pos = int(pos_s)
                fa_base = fa.fetch("22", pos - 1, pos).upper()
                if ref != fa_base:
                    mismatches.append((pos, ref, fa_base))
            self.assertEqual(
                mismatches, [],
                f"REF mismatch in {len(mismatches)} of {len(records)} "
                f"records (first 5: {mismatches[:5]})",
            )


if __name__ == "__main__":
    unittest.main()
