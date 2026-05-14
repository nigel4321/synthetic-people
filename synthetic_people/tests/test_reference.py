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

from tests._shared_cache import SHARED_TEST_CACHE_DIR  # noqa: E402

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

    def test_fetch_ref_normalises_iupac_ambiguity_to_N(self):
        # PR #83 review: IUPAC ambiguity codes (R, Y, W, S, K, M,
        # B, D, H, V) appear in some assemblies. Before the fix
        # ``fetch_ref_base`` returned them verbatim, then
        # ``choose_alt`` returned ``None`` for non-ACGT input and
        # the producer's ``assert alt is not None`` tripped. The
        # callable contract now: only ``A``/``C``/``G``/``T``/``N``
        # come out, so ``_pick_ref``'s ``N`` fallback funnels
        # every non-canonical base into the rng draw cleanly.
        from syntheticgen.reference import load_fasta, fetch_ref_base
        # Mix canonical, IUPAC, and explicit N at known positions.
        path = _write_fasta(self.dir, {"22": "ARYWSKMBDHVN"})
        fa = load_fasta(path)
        self.assertEqual(fetch_ref_base(fa, "22", 1), "A")
        for pos in range(2, 12):  # R Y W S K M B D H V
            self.assertEqual(
                fetch_ref_base(fa, "22", pos), "N",
                f"pos {pos} expected N (IUPAC), got "
                f"{fetch_ref_base(fa, '22', pos)!r}",
            )
        self.assertEqual(fetch_ref_base(fa, "22", 12), "N")

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
        validate_fasta(fa, ["22"], chr_length_mb=0.0, build="GRCh38")

    def test_validate_fasta_raises_on_missing_chrom(self):
        from syntheticgen.reference import load_fasta, validate_fasta
        path = _write_fasta(self.dir, {"22": "A" * 100})
        fa = load_fasta(path)
        with self.assertRaises(ValueError) as ctx:
            validate_fasta(
                fa, ["22", "99"], chr_length_mb=0.0, build="GRCh38",
            )
        self.assertIn("99", str(ctx.exception))
        self.assertIn("missing", str(ctx.exception).lower())

    def test_validate_fasta_raises_when_chrom_too_short(self):
        from syntheticgen.reference import load_fasta, validate_fasta
        # 1 kb FASTA but the caller wants 5 Mb of chr22; chr22's
        # natural length on GRCh38 is ~50.8 Mb, so the effective
        # required length is min(50.8M, 5M) = 5M — still raises.
        path = _write_fasta(self.dir, {"22": "A" * 1000})
        fa = load_fasta(path)
        with self.assertRaises(ValueError) as ctx:
            validate_fasta(fa, ["22"], chr_length_mb=5.0, build="GRCh38")
        self.assertIn("shorter than", str(ctx.exception))

    def test_validate_fasta_skips_length_check_when_chr_length_zero(self):
        # ``--chr-length-mb 0`` means "use full contig"; we don't
        # have a target length to enforce so length validation
        # is skipped (presence validation still applies).
        from syntheticgen.reference import load_fasta, validate_fasta
        path = _write_fasta(self.dir, {"22": "A" * 100})
        fa = load_fasta(path)
        validate_fasta(fa, ["22"], chr_length_mb=0.0, build="GRCh38")

    def test_validate_fasta_with_none_is_noop(self):
        # No FASTA passed → no validation needed.
        from syntheticgen.reference import validate_fasta
        validate_fasta(None, ["22"], chr_length_mb=5.0, build="GRCh38")

    def test_validate_fasta_passes_when_chr_length_mb_exceeds_natural(self):
        # Regression: ``--chr-length-mb`` is a CAP, not a minimum.
        # When the cap exceeds a chromosome's natural length the
        # simulator uses min(natural, cap) — so the FASTA only
        # needs to cover the natural length. Pre-fix, validate_fasta
        # rejected the run with "shorter than --chr-length-mb" for
        # every chrom <70 Mb at ``--chr-length-mb 70`` on GRCh38
        # (i.e. chr19/20/21/22 + every sex chrom).
        from unittest.mock import patch
        from syntheticgen.reference import load_fasta, validate_fasta
        path = _write_fasta(self.dir, {"22": "A" * 1000})
        fa = load_fasta(path)
        # Stub the BUILDS lookup so we can declare chr22's "natural"
        # length as exactly 1000 bp without having to ship a
        # 50-MB-of-A FASTA in test fixtures.
        synthetic_builds = {"GRCh38": {"contigs": {"22": 1000}}}
        with patch("syntheticgen.reference.BUILDS", synthetic_builds):
            # Cap = 5 Mb but natural is only 1 kb. Effective sim
            # length = min(1000, 5_000_000) = 1000. FASTA covers it.
            # Should NOT raise.
            validate_fasta(
                fa, ["22"], chr_length_mb=5.0, build="GRCh38",
            )

    def test_validate_fasta_raises_when_below_min_of_natural_and_cap(self):
        # Companion to the regression test: when the FASTA is
        # genuinely shorter than min(natural, cap), the validation
        # still fires. Catches the "wrong FASTA / truncated FASTA"
        # case the validation was originally written for.
        from unittest.mock import patch
        from syntheticgen.reference import load_fasta, validate_fasta
        path = _write_fasta(self.dir, {"22": "A" * 100})
        fa = load_fasta(path)
        # Natural = 5 kb, cap = 5 Mb → required = min(5 kb, 5 Mb) =
        # 5 kb. FASTA only has 100 bp → too short.
        synthetic_builds = {"GRCh38": {"contigs": {"22": 5000}}}
        with patch("syntheticgen.reference.BUILDS", synthetic_builds):
            with self.assertRaises(ValueError) as ctx:
                validate_fasta(
                    fa, ["22"], chr_length_mb=5.0, build="GRCh38",
                )
        self.assertIn("shorter than required", str(ctx.exception))


@unittest.skipUnless(_HAVE_PYSAM, "pysam not installed")
class FetchReferenceFastaTest(unittest.TestCase):
    """``fetch_reference_fasta`` is the M12 auto-download surface.

    The real URL is a 3 GB Ensembl FASTA — far too large to pull in
    CI on every job. The cli pairs auto-fetch with ``actions/cache``
    so it lands once per CI image and is reused thereafter. Tests
    monkey-patch ``urllib.request.urlopen`` to a local ``BytesIO``
    backed by a hand-built tiny gzipped FASTA, which exercises every
    branch in the cache contract without hitting the network.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cache = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def _install_fake_urlopen(self, payload: bytes):
        """Replace ``urllib.request.urlopen`` with a ``BytesIO`` stub
        for the lifetime of one test. Restores on tearDown.

        The stub mirrors ``HTTPResponse`` enough that
        ``reference._download`` can call ``getheader`` + ``read`` and
        use the response as a context manager.
        """
        import io
        from syntheticgen import reference as ref_mod

        class _Resp(io.BytesIO):
            def __enter__(self_inner):
                return self_inner

            def __exit__(self_inner, *exc):
                self_inner.close()
                return False

            def getheader(self_inner, name):
                if name == "Content-Length":
                    return str(len(payload))
                return None

        original = ref_mod.urllib.request.urlopen
        ref_mod.urllib.request.urlopen = lambda url: _Resp(payload)
        self.addCleanup(
            setattr, ref_mod.urllib.request, "urlopen", original,
        )

    def _tiny_gz_fasta(self) -> bytes:
        """Return a bgzip-friendly gzipped FASTA whose ungzipped form
        is a short valid FASTA with one record.

        ``gzip.GzipFile`` output decompresses cleanly with
        ``gzip.open(...)`` which is what ``fetch_reference_fasta``
        uses; we don't need a real BGZF block structure here because
        the production decompress path is plain gzip.
        """
        import gzip
        import io
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
            gz.write(b">22\nACGTACGTACGT\n")
        return buf.getvalue()

    def _install_failing_urlopen(self):
        """Replace ``urlopen`` with one that asserts it's never called.

        PR #86 review (Copilot): the prior tests relied on the absence
        of a stub to prove the download branch wasn't entered, which
        means a future refactor that reorders the cache check would
        silently hit the real Ensembl URL on CI runners with network
        access. Install an explicit fail-if-called stub instead so
        the test fails loudly on any regression.
        """
        from syntheticgen import reference as ref_mod
        original = ref_mod.urllib.request.urlopen

        def _fail(url):
            self.fail(
                f"download branch fired (urlopen called with {url!r}) — "
                "this test should be a pure cache hit",
            )

        ref_mod.urllib.request.urlopen = _fail
        self.addCleanup(
            setattr, ref_mod.urllib.request, "urlopen", original,
        )

    def test_idempotent_when_fa_and_fai_already_present(self):
        # Both .fa and .fai already on disk → return path without
        # touching urlopen.
        from syntheticgen.reference import fetch_reference_fasta
        import pysam
        ref_dir = self.cache / "reference"
        ref_dir.mkdir(parents=True)
        fa = ref_dir / "GRCh38.fa"
        fa.write_text(">22\nACGT\n")
        pysam.faidx(str(fa))
        # Sanity: .fai now exists.
        self.assertTrue(fa.with_suffix(".fa.fai").is_file())
        self._install_failing_urlopen()
        out = fetch_reference_fasta(self.cache, "GRCh38")
        self.assertEqual(out, fa)

    def test_reindexes_when_fai_missing_but_fa_present(self):
        # Cache contract: a present-but-unindexed FASTA is re-indexed
        # in place instead of redownloaded.
        from syntheticgen.reference import fetch_reference_fasta
        ref_dir = self.cache / "reference"
        ref_dir.mkdir(parents=True)
        fa = ref_dir / "GRCh38.fa"
        fa.write_text(">22\nACGT\n")
        self.assertFalse(fa.with_suffix(".fa.fai").is_file())
        self._install_failing_urlopen()
        out = fetch_reference_fasta(self.cache, "GRCh38")
        self.assertEqual(out, fa)
        self.assertTrue(fa.with_suffix(".fa.fai").is_file())

    def test_overwrites_stale_gz_from_prior_partial_run(self):
        # PR #86 review (Copilot): if a previous run died after the
        # .gz landed but before decompress finished, the next run
        # used to ``rename`` the new download over the stale .gz —
        # which fails on Windows. ``Path.replace`` + pre-emptive
        # unlink of stale sidecars handles this.
        from syntheticgen.reference import fetch_reference_fasta
        ref_dir = self.cache / "reference"
        ref_dir.mkdir(parents=True)
        # Seed a stale .gz that a prior run would have left behind.
        (ref_dir / "GRCh38.fa.gz").write_bytes(b"junk-from-prior-run")
        self._install_fake_urlopen(self._tiny_gz_fasta())
        out = fetch_reference_fasta(self.cache, "GRCh38")
        self.assertTrue(out.is_file())
        # Stale .gz must have been overwritten then cleaned up.
        self.assertFalse((ref_dir / "GRCh38.fa.gz").exists())

    def test_full_download_decompress_and_index(self):
        from syntheticgen.reference import fetch_reference_fasta
        self._install_fake_urlopen(self._tiny_gz_fasta())
        out = fetch_reference_fasta(self.cache, "GRCh38")
        self.assertEqual(out, self.cache / "reference" / "GRCh38.fa")
        self.assertTrue(out.is_file())
        self.assertTrue(out.with_suffix(".fa.fai").is_file())
        # The .gz must be cleaned up post-decompress (saves ~900 MB
        # on the real download).
        self.assertFalse(
            (self.cache / "reference" / "GRCh38.fa.gz").exists(),
        )
        # And the .part sentinels must not linger.
        self.assertFalse(
            (self.cache / "reference" / "GRCh38.fa.gz.part").exists(),
        )
        self.assertFalse(
            (self.cache / "reference" / "GRCh38.fa.part").exists(),
        )
        # Smoke-test that the resulting FASTA actually contains what
        # we shipped through the stub.
        import pysam
        fa = pysam.FastaFile(str(out))
        try:
            self.assertEqual(fa.fetch("22", 0, 4), "ACGT")
        finally:
            fa.close()

    def test_rejects_unknown_build(self):
        from syntheticgen.reference import fetch_reference_fasta
        with self.assertRaises(ValueError):
            fetch_reference_fasta(self.cache, "GRCh99")


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
            # Share the ClinVar download across tests; see
            # tests/_shared_cache.py.
            "--cache-dir", str(SHARED_TEST_CACHE_DIR),
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


_HAVE_DEMES = importlib.util.find_spec("demes") is not None
_HAVE_TSKIT = importlib.util.find_spec("tskit") is not None


@unittest.skipUnless(
    _HAVE_PYSAM and _HAVE_MSPRIME and _HAVE_DEMES and _HAVE_TSKIT,
    # PR #84 review: the admixture path's deps are msprime + demes +
    # tskit (see ``syntheticgen/admixture.py:_require_deps``).
    # ``stdpopsim`` is NOT required — the original guard here was
    # too strict and skipped the test in environments where the
    # admixture stack is present but stdpopsim is not.
    "needs pysam + msprime + demes + tskit (admixture deps)",
)
class AdmixtureReferenceTest(unittest.TestCase):
    """PR #83 review #2: the admixture path has its own
    ``simulate_chromosome`` / ``_simulate_chromosome_from_seed`` /
    ``simulate_cohort`` chain that wasn't threaded with the FASTA
    in the original M12 PR. This test exercises the threading.

    The admixture chain uses ``ProcessPoolExecutor`` for parallel
    chromosomes (when ``workers > 1`` and ``len(chromosomes) > 1``),
    so the FASTA must be passed as a path (string) rather than a
    ``FastaFile`` handle — handles don't pickle. Each worker
    re-opens from the path; kernel mmap is shared.
    """

    def _build_fasta(self, tmp: Path, chroms: tuple = ("22",)) -> Path:
        # 1 Mb of CG repeats per requested chrom — distinct bytes
        # per pos so the test catches an "all-rng-fallback" bug
        # (which would produce a uniform distribution, not all C/G).
        import pysam
        path = tmp / "admix.fa"
        body = "".join(
            f">{c}\n{'CGCG' * 250_000}\n" for c in chroms
        )
        path.write_text(body)
        pysam.faidx(str(path))
        return path

    def _assert_all_refs_cg(self, sites: list, label: str) -> None:
        non_cg = [s["ref"] for s in sites if s["ref"] not in "CG"]
        self.assertEqual(
            non_cg, [],
            f"{label}: {len(non_cg)}/{len(sites)} REFs were not C/G "
            f"— FASTA path not threaded; rng fallback fired uniformly",
        )

    def test_admixture_serial_chain_uses_fasta(self):
        # workers=1 path → serial, no ProcessPoolExecutor.
        # Verifies the path-threading works end-to-end through
        # simulate_cohort → _simulate_chromosome_from_seed →
        # simulate_chromosome → _tree_sequence_to_sites.
        import random
        from syntheticgen.admixture import simulate_cohort
        with tempfile.TemporaryDirectory() as tmp_str:
            tmp = Path(tmp_str)
            fasta_path = self._build_fasta(tmp)
            rng = random.Random(42)
            sites, ancestry = simulate_cohort(
                chromosomes=["22"], build="GRCh38",
                n_people=3, length_mb=0.5,
                proportions=(0.6, 0.25, 0.15),
                rec_rate=1e-8, mu=1.29e-8,
                rng=rng, verbose=False, workers=1,
                fasta_path=str(fasta_path),
            )
            self.assertGreater(len(sites), 0)
            self._assert_all_refs_cg(sites, "serial")

    def test_admixture_parallel_chain_pickles_fasta_path(self):
        # PR #84 review #2: workers > 1 AND len(chromosomes) > 1
        # routes through ProcessPoolExecutor. The fasta_path
        # (string) pickles cleanly; an open FastaFile would not.
        # This test exercises that branch — if a future change
        # accidentally passes a handle through the executor, it'll
        # fail to pickle and this test catches it before WGS users
        # do.
        import random
        from syntheticgen.admixture import simulate_cohort
        with tempfile.TemporaryDirectory() as tmp_str:
            tmp = Path(tmp_str)
            fasta_path = self._build_fasta(tmp, chroms=("21", "22"))
            rng = random.Random(42)
            sites, ancestry = simulate_cohort(
                chromosomes=["21", "22"], build="GRCh38",
                n_people=2, length_mb=0.3,
                proportions=(0.6, 0.25, 0.15),
                rec_rate=1e-8, mu=1.29e-8,
                rng=rng, verbose=False, workers=2,
                fasta_path=str(fasta_path),
            )
            self.assertGreater(len(sites), 0)
            self._assert_all_refs_cg(sites, "parallel (workers=2)")


@unittest.skipUnless(
    _HAVE_PYSAM and _HAVE_MSPRIME and _HAVE_DEMES and _HAVE_TSKIT
    and _BCFTOOLS_BIN,
    "needs pysam + admixture deps + bcftools",
)
class AdmixtureCliFastaTest(unittest.TestCase):
    """PR #84 review #2: end-to-end cli invocation with
    ``--reference-fasta --admixture``. The previous regression
    test in ``AdmixtureReferenceTest`` calls ``simulate_cohort``
    directly, bypassing the cli's arg-routing + pre-flight
    validation. This test runs the cli's ``main`` so a future
    regression in the cli plumbing (e.g. forgetting to pass
    ``args.reference_fasta`` through) trips here.
    """

    def _build_fasta(self, tmp: Path) -> Path:
        import pysam
        path = tmp / "cli_admix.fa"
        path.write_text(f">22\n{'CGCG' * 250_000}\n")
        pysam.faidx(str(path))
        return path

    def test_cli_admixture_emits_real_refs(self):
        import subprocess
        from syntheticgen import cli as cli_module
        with tempfile.TemporaryDirectory() as tmp_str:
            tmp = Path(tmp_str)
            fasta_path = self._build_fasta(tmp)
            out_dir = tmp / "out"
            # Minimal cli invocation: 2 people, chr22 only, 0.5 Mb,
            # admixture proportions sum to 1.0, --reference-fasta
            # threaded all the way through.
            argv = [
                "--no-config",
                "--n", "2",
                "--seed", "42",
                "--build", "GRCh38",
                "--chromosomes", "22",
                "--chr-length-mb", "0.5",
                "--admixture",
                "--eur-frac", "0.6",
                "--sas-frac", "0.25",
                "--afr-frac", "0.15",
                "--rsid-density", "0",
                "--clinvar-inject-density", "0",
                "--svs-per-person", "0",
                "--error-rate", "0",
                "--dropout-rate", "0",
                "--workers", "1",
                "--output-dir", str(out_dir),
                # Share the ClinVar download across tests; see
                # tests/_shared_cache.py.
                "--cache-dir", str(SHARED_TEST_CACHE_DIR),
                "--reference-fasta", str(fasta_path),
            ]
            rc = cli_module.main(argv)
            self.assertEqual(rc, 0)
            # Pull every REF out of the first person's VCF;
            # they must all be C or G.
            vcf = out_dir / "person_0001.vcf.gz"
            self.assertTrue(vcf.is_file(), f"person VCF missing: {vcf}")
            out = subprocess.check_output(
                ["bcftools", "query", "-f", "%REF\n", str(vcf)],
                text=True,
            )
            refs = [r for r in out.strip().splitlines() if r]
            self.assertGreater(len(refs), 0)
            non_cg = [r for r in refs if r not in ("C", "G")]
            self.assertEqual(
                non_cg, [],
                f"cli admixture: {len(non_cg)}/{len(refs)} non-C/G "
                f"REFs — the --reference-fasta flag was not routed "
                f"into the admixture path",
            )

    def test_cli_admixture_rejects_missing_fasta_early(self):
        # PR #84 review #3: a missing FASTA path must fail at
        # startup, not after msprime starts working inside the
        # admixture worker. The streamed coalescent path
        # validates at startup; this test ensures the admixture
        # path now does too (PR #84 added the validate-then-discard
        # block to the admixture branch in cli.main).
        #
        # PR #86 review (Copilot): the cli now surfaces the
        # ``FileNotFoundError`` as a ``SystemExit`` with an
        # actionable hint mentioning the ``--no-reference-fasta``
        # opt-out, rather than letting the bare stack trace
        # bubble up.
        from syntheticgen import cli as cli_module
        with tempfile.TemporaryDirectory() as tmp_str:
            tmp = Path(tmp_str)
            out_dir = tmp / "out"
            argv = [
                "--no-config",
                "--n", "2",
                "--seed", "42",
                "--build", "GRCh38",
                "--chromosomes", "22",
                "--chr-length-mb", "0.5",
                "--admixture",
                "--eur-frac", "0.6",
                "--sas-frac", "0.25",
                "--afr-frac", "0.15",
                "--rsid-density", "0",
                "--clinvar-inject-density", "0",
                "--svs-per-person", "0",
                "--workers", "1",
                "--output-dir", str(out_dir),
                # Share the ClinVar download across tests; see
                # tests/_shared_cache.py.
                "--cache-dir", str(SHARED_TEST_CACHE_DIR),
                "--reference-fasta", "/nonexistent/path.fa",
            ]
            with self.assertRaises(SystemExit) as ctx:
                cli_module.main(argv)
            self.assertIn("--no-reference-fasta", str(ctx.exception))


class ResolveReferenceFastaFlagConflictTest(unittest.TestCase):
    """PR #86 review (Copilot): ``--no-reference-fasta`` and an
    explicit ``--reference-fasta <path>`` are contradictory. The
    resolver must reject the combination loudly so users don't
    silently get the legacy fabricated-REF path while wondering
    why their FASTA path was ignored.
    """

    def test_both_flags_set_exits_with_clear_message(self):
        from syntheticgen.cli import _resolve_reference_fasta
        import argparse
        args = argparse.Namespace(
            no_reference_fasta=True,
            reference_fasta=Path("/tmp/whatever.fa"),
            cache_dir=Path("/tmp/cache"),
            build="GRCh38",
            chr_length_mb=0.0,
        )
        with self.assertRaises(SystemExit) as ctx:
            _resolve_reference_fasta(args, ["22"])
        msg = str(ctx.exception)
        # Message must name both flags so the user knows which one
        # to drop.
        self.assertIn("--no-reference-fasta", msg)
        self.assertIn("--reference-fasta", msg)
        self.assertIn("mutually exclusive", msg)

    def test_only_no_reference_fasta_returns_none_pair(self):
        # Sanity: with just --no-reference-fasta (and no
        # --reference-fasta), the resolver returns the (None, None)
        # opt-out tuple — i.e. the flag-conflict check doesn't
        # mis-fire on the common opt-out path.
        from syntheticgen.cli import _resolve_reference_fasta
        import argparse
        args = argparse.Namespace(
            no_reference_fasta=True,
            reference_fasta=None,
            cache_dir=Path("/tmp/cache"),
            build="GRCh38",
            chr_length_mb=0.0,
        )
        fa, path = _resolve_reference_fasta(args, ["22"])
        self.assertIsNone(fa)
        self.assertIsNone(path)


@unittest.skipUnless(_HAVE_PYSAM, "pysam not installed")
class ResolveReferenceFastaActionableErrorsTest(unittest.TestCase):
    """PR #86 review (Copilot + claude review): when ``load_fasta``
    or ``validate_fasta`` raises, the cli must surface an actionable
    message (mentioning ``--no-reference-fasta``) and not leak the
    open ``pysam.FastaFile`` handle through ``sys.exit``.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_missing_fasta_path_exits_with_install_hint(self):
        from syntheticgen.cli import _resolve_reference_fasta
        import argparse
        args = argparse.Namespace(
            no_reference_fasta=False,
            reference_fasta=Path("/definitely/does/not/exist.fa"),
            cache_dir=self.dir,
            build="GRCh38",
            chr_length_mb=0.0,
        )
        with self.assertRaises(SystemExit) as ctx:
            _resolve_reference_fasta(args, ["22"])
        msg = str(ctx.exception)
        self.assertIn("--no-reference-fasta", msg)

    def test_validate_failure_closes_handle_before_exit(self):
        # The FASTA below is a valid FASTA but only contains chr22
        # at 8 bp; ask for a chromosome that isn't there so
        # validate_fasta raises. The just-opened pysam handle must
        # be closed before sys.exit() rather than relying on GC.
        from syntheticgen.cli import _resolve_reference_fasta
        import argparse
        import pysam
        fa_path = self.dir / "tiny.fa"
        fa_path.write_text(">22\nACGTACGT\n")
        pysam.faidx(str(fa_path))
        args = argparse.Namespace(
            no_reference_fasta=False,
            reference_fasta=fa_path,
            cache_dir=self.dir,
            build="GRCh38",
            chr_length_mb=0.0,
        )
        # Capture the FastaFile instance via the reference module's
        # ``load_fasta`` so we can interrogate it after sys.exit.
        # ``cli._resolve_reference_fasta`` does ``from .reference
        # import load_fasta`` inside the function body, so it
        # re-fetches the module attribute at call time — patching
        # the reference module is sufficient.
        from unittest.mock import patch
        captured: list = []
        from syntheticgen import reference as ref_mod
        original_load = ref_mod.load_fasta

        def _wrapped(path):
            handle = original_load(path)
            captured.append(handle)
            return handle

        with patch.object(ref_mod, "load_fasta", _wrapped):
            with self.assertRaises(SystemExit):
                _resolve_reference_fasta(args, ["99"])
        # validate raised → handle must have been closed before exit.
        self.assertEqual(len(captured), 1)
        handle = captured[0]
        # pysam.FastaFile exposes ``.is_open`` (False after close).
        self.assertFalse(handle.is_open())

    def test_auto_fetch_import_error_exits_with_install_hint(self):
        # PR #87 review (Copilot): ``fetch_reference_fasta`` runs
        # BEFORE ``load_fasta`` on the auto-fetch path and can raise
        # its own ImportError (pysam needed for ``pysam.faidx``).
        # That must surface the same actionable hint, not a bare
        # traceback.
        from syntheticgen.cli import _resolve_reference_fasta
        from syntheticgen import reference as ref_mod
        from unittest.mock import patch
        import argparse
        args = argparse.Namespace(
            no_reference_fasta=False,
            reference_fasta=None,  # forces auto-fetch path
            cache_dir=self.dir,
            build="GRCh38",
            chr_length_mb=0.0,
        )

        def _raise_import_error(cache_dir, build):
            raise ImportError("pysam is required to index the downloaded FASTA")

        with patch.object(ref_mod, "fetch_reference_fasta", _raise_import_error):
            with self.assertRaises(SystemExit) as ctx:
                _resolve_reference_fasta(args, ["22"])
        msg = str(ctx.exception)
        self.assertIn("--no-reference-fasta", msg)
        self.assertIn("pip install pysam", msg)

    def test_load_fasta_value_error_exits_with_hint(self):
        # PR #87 review (Copilot): ``load_fasta`` raises ValueError
        # when the FASTA is unreadable or unindexed (e.g. .fai
        # missing on a read-only filesystem). That path must also
        # surface the actionable hint and not escape as a traceback.
        from syntheticgen.cli import _resolve_reference_fasta
        from syntheticgen import reference as ref_mod
        from unittest.mock import patch
        import argparse
        fa_path = self.dir / "broken.fa"
        fa_path.write_text(">22\nACGT\n")
        args = argparse.Namespace(
            no_reference_fasta=False,
            reference_fasta=fa_path,
            cache_dir=self.dir,
            build="GRCh38",
            chr_length_mb=0.0,
        )

        def _raise_value_error(path):
            raise ValueError(
                f"--reference-fasta: could not open {path} — index missing."
            )

        with patch.object(ref_mod, "load_fasta", _raise_value_error):
            with self.assertRaises(SystemExit) as ctx:
                _resolve_reference_fasta(args, ["22"])
        msg = str(ctx.exception)
        self.assertIn("--no-reference-fasta", msg)


if __name__ == "__main__":
    unittest.main()
