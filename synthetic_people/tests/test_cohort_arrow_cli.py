"""Tests for the Phase 5d.1 cli.py orchestration: --cohort-mode
{sites_list,arrow,auto}, the pre-flight disk-space check, and the
full-pipeline parity claim.

Pure unit tests run anywhere. The full-pipeline parity test gates
on bcftools / tabix / bgzip / msprime / stdpopsim / pyarrow.
"""

from __future__ import annotations

import collections
import hashlib
import importlib.util
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen import cli as cli_module  # noqa: E402
from syntheticgen.cli import (  # noqa: E402
    _ARROW_AUTO_THRESHOLD,
    _estimate_arrow_chrom_scratch_bytes,
    _preflight_arrow_disk_check,
    _resolve_cohort_mode,
)

_HAVE_BCFTOOLS = shutil.which("bcftools") is not None
_HAVE_TABIX = shutil.which("tabix") is not None
_HAVE_BGZIP = shutil.which("bgzip") is not None
_HAVE_MSPRIME = importlib.util.find_spec("msprime") is not None
_HAVE_STDPOPSIM = importlib.util.find_spec("stdpopsim") is not None
_HAVE_PYARROW = importlib.util.find_spec("pyarrow") is not None


class ResolveCohortModeTest(unittest.TestCase):
    """Mode resolution is the only place ``--cohort-mode auto`` becomes
    a concrete mode. The threshold is the only thing the rest of the
    codebase depends on for picking; pin it explicitly."""

    def test_explicit_sites_list_overrides_n(self):
        self.assertEqual(_resolve_cohort_mode("sites_list", 1), "sites_list")
        self.assertEqual(
            _resolve_cohort_mode("sites_list", 10_000_000), "sites_list",
        )

    def test_explicit_arrow_overrides_n(self):
        self.assertEqual(_resolve_cohort_mode("arrow", 1), "arrow")
        self.assertEqual(
            _resolve_cohort_mode("arrow", 10_000_000), "arrow",
        )

    def test_auto_below_threshold_picks_sites_list(self):
        self.assertEqual(_resolve_cohort_mode("auto", 1), "sites_list")
        self.assertEqual(_resolve_cohort_mode("auto", 99_999), "sites_list")
        self.assertEqual(
            _resolve_cohort_mode("auto", _ARROW_AUTO_THRESHOLD - 1),
            "sites_list",
        )

    def test_auto_at_threshold_picks_arrow(self):
        self.assertEqual(
            _resolve_cohort_mode("auto", _ARROW_AUTO_THRESHOLD), "arrow",
        )
        self.assertEqual(_resolve_cohort_mode("auto", 1_000_000), "arrow")

    def test_threshold_is_100k(self):
        # The threshold value itself is part of the public contract — pin it
        # so any future change is a deliberate decision, not an accident.
        self.assertEqual(_ARROW_AUTO_THRESHOLD, 100_000)


class EstimateArrowScratchTest(unittest.TestCase):
    def test_scales_with_n(self):
        a = _estimate_arrow_chrom_scratch_bytes(1_000, 1.0)
        b = _estimate_arrow_chrom_scratch_bytes(2_000, 1.0)
        self.assertEqual(b, 2 * a)

    def test_scales_with_chr_length(self):
        a = _estimate_arrow_chrom_scratch_bytes(1_000, 1.0)
        b = _estimate_arrow_chrom_scratch_bytes(1_000, 10.0)
        self.assertEqual(b, 10 * a)

    def test_n1m_70mb_in_expected_range(self):
        # Plan §5d says ~250-300 GB after Arrow encoding for n=1M ×
        # 70Mb. The estimator over-estimates raw uncompressed (which
        # is the safe direction for a pre-flight check).
        b = _estimate_arrow_chrom_scratch_bytes(1_000_000, 70.0)
        self.assertGreater(b, 100 * 1e9)   # > 100 GB
        self.assertLess(b, 1500 * 1e9)     # < 1.5 TB upper sanity

    def test_arrow_mode_excludes_sidecar_term(self):
        # PR #77 review: only ``arrow-streaming`` allocates a
        # carriers sidecar. ``arrow`` mode keeps carriers in RAM
        # and must not be charged for it — would falsely fail
        # tight-disk runs.
        arrow_only = _estimate_arrow_chrom_scratch_bytes(
            1_000, 1.0, cohort_mode="arrow",
        )
        streaming = _estimate_arrow_chrom_scratch_bytes(
            1_000, 1.0, cohort_mode="arrow-streaming",
        )
        # Streaming budget is roughly 2× arrow (one Arrow file +
        # one sidecar of comparable size).
        self.assertAlmostEqual(streaming, 2 * arrow_only)

    def test_default_cohort_mode_is_arrow(self):
        # Defensive: callers passing no cohort_mode should get the
        # arrow-only (no-sidecar) estimate, not the wider streaming
        # one. The sidecar charge has to be opt-in for the runtime
        # not to surprise legacy callers.
        defaulted = _estimate_arrow_chrom_scratch_bytes(1_000, 1.0)
        arrow_only = _estimate_arrow_chrom_scratch_bytes(
            1_000, 1.0, cohort_mode="arrow",
        )
        self.assertEqual(defaulted, arrow_only)


class PreflightArrowDiskCheckTest(unittest.TestCase):
    """The pre-flight check has three observable behaviours: pass
    silently when free disk >= 2x estimate, warn when it's between
    1x and 2x, fail-with-SystemExit when below 1x."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cohort_dir = Path(self.tmp.name) / "cohort"

    def tearDown(self):
        self.tmp.cleanup()

    def _patch_free_bytes(self, free_bytes: int):
        DiskUsage = collections.namedtuple(
            "DiskUsage", "total used free",
        )
        return mock.patch.object(
            cli_module.shutil, "disk_usage",
            return_value=DiskUsage(total=10**12, used=0, free=free_bytes),
        )

    def test_plenty_of_free_disk_passes_silently(self):
        with self._patch_free_bytes(10**12):  # 1 TB free
            with mock.patch.object(cli_module.sys, "stderr") as err:
                _preflight_arrow_disk_check(self.cohort_dir, 100, 0.05)
                err.write.assert_not_called()

    def test_tight_free_disk_warns(self):
        # Estimate at n=100, chr=0.05Mb is 100 × 0.05 × 5000 × 2 × 1.05
        # = ~52 KB. Set free to 1.5x of that so it warns but doesn't fail.
        per_chrom = _estimate_arrow_chrom_scratch_bytes(100, 0.05)
        with self._patch_free_bytes(int(per_chrom * 1.5)):
            with mock.patch.object(cli_module.sys, "stderr") as err:
                _preflight_arrow_disk_check(self.cohort_dir, 100, 0.05)
                # At least one stderr write happened with the warning.
                self.assertTrue(err.write.called)
                printed = "".join(
                    a.args[0] if a.args else ""
                    for a in err.write.call_args_list
                )
                self.assertIn("WARNING", printed)
                self.assertIn("arrow", printed)

    def test_insufficient_free_disk_exits(self):
        per_chrom = _estimate_arrow_chrom_scratch_bytes(100, 0.05)
        with self._patch_free_bytes(per_chrom // 2):  # below 1x
            with self.assertRaises(SystemExit) as ctx:
                _preflight_arrow_disk_check(self.cohort_dir, 100, 0.05)
            msg = str(ctx.exception)
            self.assertIn("--cohort-mode arrow needs", msg)
            self.assertIn("Free up disk", msg)

    def test_streaming_mode_message_mentions_sidecar(self):
        # PR #78 review #1: the failure message used to hardcode
        # ``--cohort-mode arrow`` even when ``arrow-streaming``
        # triggered the check. Now ``cohort_mode`` is interpolated
        # AND the streaming variant mentions the sidecar so
        # operators understand why the budget is ~2× the
        # materialised-arrow budget.
        per_chrom = _estimate_arrow_chrom_scratch_bytes(
            100, 0.05, cohort_mode="arrow-streaming",
        )
        with self._patch_free_bytes(per_chrom // 2):
            with self.assertRaises(SystemExit) as ctx:
                _preflight_arrow_disk_check(
                    self.cohort_dir, 100, 0.05,
                    cohort_mode="arrow-streaming",
                )
            msg = str(ctx.exception)
            self.assertIn("--cohort-mode arrow-streaming needs", msg)
            # Sidecar is the user-visible reason the streaming
            # estimate is wider than the materialised one.
            self.assertIn("sidecar", msg)
            self.assertIn("sites_list", msg)


def _common_args(
    out_dir: Path,
    cohort_mode: str,
    n: int = 3,
    chroms: str = "22",
) -> list:
    """Mirrors test_cohort_streaming._common_args; adds --cohort-mode."""
    return [
        "--n", str(n),
        "--seed", "42",
        "--build", "GRCh38",
        "--chromosomes", chroms,
        "--chr-length-mb", "0.05",
        "--demo-model", "none",
        "--rsid-density", "0",
        "--clinvar-inject-density", "0",
        "--svs-per-person", "0",
        "--error-rate", "0",
        "--dropout-rate", "0",
        "--workers", "2",
        "--output-dir", str(out_dir),
        "--cache-dir", str(out_dir / "cache"),
        "--mode", "cohort",
        "--cohort-mode", cohort_mode,
    ]


def _bcf_data_md5(bcf_path: Path) -> str:
    """Same content-hash helper as test_cohort_arrow_bcf_writer /
    test_cohort_parallel_write — fixed-order line format via
    bcftools query so the hash is insensitive to BCF metadata trivia."""
    out = subprocess.check_output(
        ["bcftools", "query",
         "-f", "%CHROM\t%POS\t%REF\t%ALT\t"
               "%INFO/AC\t%INFO/AN\t%INFO/AF[\t%GT]\n",
         str(bcf_path)])
    return hashlib.md5(out).hexdigest()


@unittest.skipUnless(
    _HAVE_BCFTOOLS and _HAVE_TABIX and _HAVE_BGZIP
    and _HAVE_MSPRIME and _HAVE_STDPOPSIM and _HAVE_PYARROW,
    "needs bcftools/tabix/bgzip + msprime + stdpopsim + pyarrow",
)
class CohortModeArrowParityTest(unittest.TestCase):
    """The integration parity claim: end-to-end cli runs at the same
    seed produce the same cohort BCF content under
    ``--cohort-mode sites_list`` and ``--cohort-mode arrow``.

    Intentionally tiny scale (n=3, 0.05 Mb, chr22) so the two runs
    finish in seconds; the byte-identical assertion is what carries
    the load-bearing claim, not coverage of large-n behaviour.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir_a = Path(
            tempfile.mkdtemp(prefix="cohort_mode_sites_list_"),
        )
        cls.tmpdir_b = Path(tempfile.mkdtemp(prefix="cohort_mode_arrow_"))
        rc_a = cli_module.main(
            _common_args(cls.tmpdir_a, "sites_list"),
        )
        if rc_a != 0:
            raise RuntimeError(f"cli.main(sites_list) exited {rc_a}")
        rc_b = cli_module.main(_common_args(cls.tmpdir_b, "arrow"))
        if rc_b != 0:
            raise RuntimeError(f"cli.main(arrow) exited {rc_b}")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir_a, ignore_errors=True)
        shutil.rmtree(cls.tmpdir_b, ignore_errors=True)

    def test_byte_identical_cohort_bcf(self):
        bcf_a = self.tmpdir_a / "cohort" / "cohort.chr22.bcf"
        bcf_b = self.tmpdir_b / "cohort" / "cohort.chr22.bcf"
        self.assertTrue(bcf_a.exists(), f"missing {bcf_a}")
        self.assertTrue(bcf_b.exists(), f"missing {bcf_b}")
        self.assertEqual(_bcf_data_md5(bcf_a), _bcf_data_md5(bcf_b))

    def test_arrow_scratch_cleaned_up_on_success(self):
        # The Arrow path creates cohort/.arrow/cohort.chr<N>.arrow
        # transiently; on success both the file and the empty .arrow
        # directory should be gone.
        arrow_dir = self.tmpdir_b / "cohort" / ".arrow"
        self.assertFalse(
            arrow_dir.exists(),
            f"Arrow scratch dir was not cleaned up: {arrow_dir}",
        )

    def test_sites_list_run_does_not_create_arrow_scratch(self):
        # sites_list mode must never touch the .arrow scratch.
        arrow_dir = self.tmpdir_a / "cohort" / ".arrow"
        self.assertFalse(arrow_dir.exists())


@unittest.skipUnless(
    _HAVE_BCFTOOLS and _HAVE_TABIX and _HAVE_BGZIP
    and _HAVE_MSPRIME and _HAVE_STDPOPSIM and _HAVE_PYARROW,
    "needs bcftools/tabix/bgzip + msprime + stdpopsim + pyarrow",
)
class CohortModeArrowStreamingParityTest(unittest.TestCase):
    """Phase 5d.1 PR 3: ``--cohort-mode arrow-streaming`` end-to-end.

    Same byte-identical contract as :class:`CohortModeArrowParityTest`:
    a cli run with ``--cohort-mode arrow-streaming`` at the same seed
    produces the exact same cohort BCF as ``--cohort-mode arrow`` (the
    materialised path). This locks the streaming refactor as a pure
    performance change with zero observable output difference.

    Intentionally tiny scale (n=3, 0.05 Mb, chr22) — the byte-
    identical assertion is what matters; coverage of WGS-scale
    behaviour lives in the unit-level streaming-cohort tests + the
    memprof-based smoke runs the operator triggers manually.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir_arrow = Path(
            tempfile.mkdtemp(prefix="cohort_mode_arrow_"),
        )
        cls.tmpdir_stream = Path(
            tempfile.mkdtemp(prefix="cohort_mode_arrow_streaming_"),
        )
        rc_arrow = cli_module.main(_common_args(cls.tmpdir_arrow, "arrow"))
        if rc_arrow != 0:
            raise RuntimeError(f"cli.main(arrow) exited {rc_arrow}")
        rc_stream = cli_module.main(
            _common_args(cls.tmpdir_stream, "arrow-streaming"),
        )
        if rc_stream != 0:
            raise RuntimeError(
                f"cli.main(arrow-streaming) exited {rc_stream}",
            )

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir_arrow, ignore_errors=True)
        shutil.rmtree(cls.tmpdir_stream, ignore_errors=True)

    def test_byte_identical_cohort_bcf(self):
        bcf_arrow = self.tmpdir_arrow / "cohort" / "cohort.chr22.bcf"
        bcf_stream = self.tmpdir_stream / "cohort" / "cohort.chr22.bcf"
        self.assertTrue(bcf_arrow.exists(), f"missing {bcf_arrow}")
        self.assertTrue(bcf_stream.exists(), f"missing {bcf_stream}")
        self.assertEqual(_bcf_data_md5(bcf_arrow),
                         _bcf_data_md5(bcf_stream))

    def test_streaming_run_cleans_up_arrow_scratch(self):
        arrow_dir = self.tmpdir_stream / "cohort" / ".arrow"
        self.assertFalse(
            arrow_dir.exists(),
            f"streaming Arrow scratch dir was not cleaned up: {arrow_dir}",
        )


class ResolveCohortModeAutoPickStreamingTest(unittest.TestCase):
    """Pin the auto-pick threshold logic: when the predicted
    materialised parent peak exceeds 50% of host RAM, ``--cohort-mode
    auto`` resolves to ``arrow-streaming``. Below that, the existing
    n>=100k → ``arrow`` / else → ``sites_list`` cascade applies."""

    def test_explicit_arrow_streaming_passes_through(self):
        from syntheticgen.cli import _resolve_cohort_mode
        self.assertEqual(
            _resolve_cohort_mode("arrow-streaming", 1), "arrow-streaming",
        )

    def test_auto_picks_arrow_streaming_when_predicted_peak_exceeds_half_ram(self):
        from syntheticgen.cli import _resolve_cohort_mode
        # n=3000 × 70 Mb predicts ~18 GB parent peak. With a 32 GB
        # host (half = 16 GB), auto must pick arrow-streaming.
        self.assertEqual(
            _resolve_cohort_mode(
                "auto", 3000, chr_length_mb=70,
                host_ram_bytes=32 * 1024**3,
            ),
            "arrow-streaming",
        )

    def test_auto_picks_arrow_for_large_n_when_peak_below_half_ram(self):
        from syntheticgen.cli import _resolve_cohort_mode
        # n=100k × 1 Mb predicts ~9 GB. On a 64 GB host (half = 32 GB),
        # we stay on the materialised+arrow path.
        self.assertEqual(
            _resolve_cohort_mode(
                "auto", 100_000, chr_length_mb=1,
                host_ram_bytes=64 * 1024**3,
            ),
            "arrow",
        )

    def test_auto_picks_sites_list_at_small_n_low_chr_length(self):
        from syntheticgen.cli import _resolve_cohort_mode
        # Small everything on a big host — stays on the in-RAM path.
        self.assertEqual(
            _resolve_cohort_mode(
                "auto", 1_000, chr_length_mb=5,
                host_ram_bytes=128 * 1024**3,
            ),
            "sites_list",
        )


class ResolveCohortModeFullLengthAutoPickTest(unittest.TestCase):
    """Regression test for PR #52 Copilot review comments 1 & 4.

    Pre-fix bug: ``--chr-length-mb 0`` (full contig length) skipped
    the predicted-peak check entirely, so the auto-pick fell through
    to the n>=100k cascade. For WGS-n=3000 (the exact scenario the
    streaming refactor was built for) this auto-picked ``sites_list``
    — the OOM-prone path. Now ``_resolve_cohort_mode`` derives the
    effective per-chrom length from the build's contig table when
    ``chr_length_mb <= 0``."""

    def test_full_length_wgs_n3000_picks_streaming_on_32gb_host(self):
        from syntheticgen.cli import _resolve_cohort_mode
        # chr1 on GRCh38 is ~249 Mb. n=3000 × 249 Mb × ~90 KB/sample/Mb
        # ≈ 67 GB predicted parent peak. On a 32 GB host this is well
        # above the 50% threshold (16 GB), so auto must pick streaming.
        self.assertEqual(
            _resolve_cohort_mode(
                "auto", 3_000, chr_length_mb=0,
                chromosomes=["1"], build="GRCh38",
                host_ram_bytes=32 * 1024**3,
            ),
            "arrow-streaming",
        )

    def test_full_length_chr22_only_small_n_stays_on_sites_list(self):
        from syntheticgen.cli import _resolve_cohort_mode
        # chr22 on GRCh38 is ~51 Mb. n=1000 × 51 × 90 KB ≈ 4.6 GB
        # — well under the 50% threshold on any host. Falls through
        # to the n>=100k cascade; n=1000 picks sites_list.
        self.assertEqual(
            _resolve_cohort_mode(
                "auto", 1_000, chr_length_mb=0,
                chromosomes=["22"], build="GRCh38",
                host_ram_bytes=32 * 1024**3,
            ),
            "sites_list",
        )

    def test_full_length_multi_chrom_uses_largest_as_binding(self):
        from syntheticgen.cli import _resolve_cohort_mode
        # User passes --chromosomes 1,22. We process serially, so the
        # per-chrom peak is bound by chr1 (longest). The effective
        # length should be ~249 Mb, not the average or chr22.
        self.assertEqual(
            _resolve_cohort_mode(
                "auto", 3_000, chr_length_mb=0,
                chromosomes=["1", "22"], build="GRCh38",
                host_ram_bytes=32 * 1024**3,
            ),
            "arrow-streaming",
        )

    def test_explicit_length_still_used_when_provided(self):
        from syntheticgen.cli import _resolve_cohort_mode
        # When --chr-length-mb is set, it overrides the contig lookup
        # — this is the original (pre-fix) behaviour, preserved.
        self.assertEqual(
            _resolve_cohort_mode(
                "auto", 3_000, chr_length_mb=70,
                chromosomes=["1"], build="GRCh38",
                host_ram_bytes=32 * 1024**3,
            ),
            "arrow-streaming",
        )


class ResolveCohortModeChunkingInteractionTest(unittest.TestCase):
    """Regression test for PR #52 Copilot review comment 2.

    Pre-fix bug: ``--cohort-mode arrow-streaming`` (or auto picking
    it) combined with ``--chr-chunk-mb > 0`` that would split the
    simulation would crash inside ``simulate_chromosome_ts`` with
    ``NotImplementedError`` mid-run. Now the auto-pick avoids
    streaming when chunking would split, and the cli pre-flight
    catches the explicit-explicit combination early."""

    def test_auto_skips_streaming_when_chunking_would_split(self):
        from syntheticgen.cli import _resolve_cohort_mode
        # WGS-n3000 scenario that would normally pick streaming, but
        # the user requested --chr-chunk-mb 5 (which would split a
        # 249 Mb chromosome 50 ways). Auto must fall back to arrow.
        self.assertEqual(
            _resolve_cohort_mode(
                "auto", 3_000, chr_length_mb=0,
                chromosomes=["1"], build="GRCh38",
                chunk_size_mb=5,
                host_ram_bytes=32 * 1024**3,
            ),
            "arrow",
        )

    def test_auto_picks_streaming_when_chunk_equal_or_above_length(self):
        from syntheticgen.cli import _resolve_cohort_mode
        # --chr-chunk-mb >= chrom length is a no-op (no actual split).
        # Auto-pick should not avoid streaming in this case.
        self.assertEqual(
            _resolve_cohort_mode(
                "auto", 3_000, chr_length_mb=70,
                chromosomes=["22"], build="GRCh38",
                chunk_size_mb=70,  # equal to chr_length
                host_ram_bytes=32 * 1024**3,
            ),
            "arrow-streaming",
        )

    def test_explicit_streaming_passes_through_regardless_of_chunking(self):
        from syntheticgen.cli import _resolve_cohort_mode
        # _resolve_cohort_mode itself doesn't enforce the explicit-
        # combination check; the cli main does that pre-flight.
        # The function just passes through.
        self.assertEqual(
            _resolve_cohort_mode(
                "arrow-streaming", 3_000, chr_length_mb=0,
                chromosomes=["1"], build="GRCh38",
                chunk_size_mb=5,
            ),
            "arrow-streaming",
        )


class EffectiveChrLengthTest(unittest.TestCase):
    def test_explicit_chr_length_passes_through(self):
        from syntheticgen.cli import _effective_chr_length_mb
        self.assertEqual(
            _effective_chr_length_mb(70.0, ["1"], "GRCh38"), 70.0,
        )

    def test_zero_chr_length_resolves_to_max_contig(self):
        from syntheticgen.cli import _effective_chr_length_mb
        # GRCh38 chr1 ≈ 248.96 Mb
        result = _effective_chr_length_mb(0, ["1"], "GRCh38")
        self.assertGreater(result, 240.0)
        self.assertLess(result, 260.0)

    def test_zero_chr_length_with_multi_chrom_uses_max(self):
        from syntheticgen.cli import _effective_chr_length_mb
        # max(chr22 ~51, chr1 ~249) → chr1's length
        result = _effective_chr_length_mb(0, ["1", "22"], "GRCh38")
        self.assertGreater(result, 240.0)

    def test_zero_chr_length_no_chromosomes_returns_zero(self):
        from syntheticgen.cli import _effective_chr_length_mb
        # Defensive fall-through: returns 0 so callers' heuristics
        # skip chr_length-driven checks rather than misfire.
        self.assertEqual(
            _effective_chr_length_mb(0, None, None), 0.0,
        )
        self.assertEqual(
            _effective_chr_length_mb(0, [], "GRCh38"), 0.0,
        )


class PreflightDiskCheckUsesEffectiveLengthTest(unittest.TestCase):
    """Regression test for PR #55 Copilot review comment 1.

    Pre-fix bug: cli.main called ``_preflight_arrow_disk_check`` with
    raw ``args.chr_length_mb``. When the user passed ``--chr-length-mb
    0`` (full contig) the scratch estimate collapsed to ~0 bytes and
    the pre-flight became a silent no-op for exactly the WGS-scale
    runs that needed it most. Post-fix wraps the length through
    ``_effective_chr_length_mb`` so the check sees the real per-chrom
    length (~249 Mb chr1 on GRCh38)."""

    def test_raw_zero_length_estimate_is_trivial(self):
        # Documents the broken path: at chr_length=0 the scratch
        # estimator clamps to ``max(1, variants)`` and the resulting
        # byte count is in the low-KB range, which defeats the disk
        # check on any host with non-empty filesystem.
        broken = _estimate_arrow_chrom_scratch_bytes(3000, 0)
        self.assertLess(broken, 100_000)  # < 100 KB

    def test_effective_length_estimate_is_multi_gigabyte(self):
        # Post-fix: the effective length (chr1 ≈ 249 Mb) gives a
        # scratch estimate in the multi-GB range for WGS-scale n,
        # which actually catches scratch-disk pressure.
        from syntheticgen.cli import _effective_chr_length_mb
        eff_len = _effective_chr_length_mb(0, ["1"], "GRCh38")
        self.assertGreater(eff_len, 240.0)
        fixed = _estimate_arrow_chrom_scratch_bytes(3000, eff_len)
        self.assertGreater(fixed, 5 * 1024**3)  # > 5 GB

    def test_preflight_with_resolved_length_trips_disk_check(self):
        # End-to-end the contract: when we pass the resolved length
        # into _preflight_arrow_disk_check and there's no free disk,
        # SystemExit fires. Without the fix (passing 0) it would
        # silently pass.
        from syntheticgen.cli import _effective_chr_length_mb
        eff_len = _effective_chr_length_mb(0, ["1"], "GRCh38")
        with tempfile.TemporaryDirectory() as tmp:
            cohort_dir = Path(tmp) / "cohort"
            DiskUsage = collections.namedtuple(
                "DiskUsage", "total used free",
            )
            with mock.patch.object(
                cli_module.shutil, "disk_usage",
                return_value=DiskUsage(
                    total=10**12, used=0, free=1024,  # 1 KB free
                ),
            ):
                with self.assertRaises(SystemExit):
                    _preflight_arrow_disk_check(
                        cohort_dir, 3000, eff_len,
                    )


class CheckCohortModeChunkingCompatTest(unittest.TestCase):
    """Helper extracted from cli main so the streaming-+-chunking
    guard rails are unit-testable. PR #55 Copilot review comment 2
    flagged that the explicit-mode hard-error path and the auto-mode
    demote path were untested."""

    def test_resolved_mode_not_streaming_is_noop(self):
        from syntheticgen.cli import _check_cohort_mode_chunking_compat
        # When the resolved mode isn't arrow-streaming the helper
        # short-circuits to a no-op regardless of chunk size.
        final, msg = _check_cohort_mode_chunking_compat(
            cli_cohort_mode="arrow",
            resolved_cohort_mode="arrow",
            chunk_size_mb=5,
            chromosomes=["1"],
            build="GRCh38",
            chr_length_mb=0,
        )
        self.assertEqual(final, "arrow")
        self.assertIsNone(msg)

    def test_streaming_with_chunk_zero_is_noop(self):
        from syntheticgen.cli import _check_cohort_mode_chunking_compat
        # At the pre-flight call site (before chunk auto-pick) the
        # cli passes ``args.chr_chunk_mb`` directly — and the
        # user-facing flag treats 0 as the *auto-pick sentinel*, not
        # as "no chunking". The helper must treat 0 as "chunk not
        # yet resolved, defer": only an explicit user-supplied
        # ``--chr-chunk-mb`` should be able to fire ERROR at
        # pre-flight; the auto-picked-chunk case is caught later by
        # the second-pass call site after chunk_size_mb is set.
        final, msg = _check_cohort_mode_chunking_compat(
            cli_cohort_mode="arrow-streaming",
            resolved_cohort_mode="arrow-streaming",
            chunk_size_mb=0,
            chromosomes=["1"],
            build="GRCh38",
            chr_length_mb=0,
        )
        self.assertEqual(final, "arrow-streaming")
        self.assertIsNone(msg)

    def test_streaming_with_full_chunk_is_noop(self):
        from syntheticgen.cli import _check_cohort_mode_chunking_compat
        # Chunk >= chrom length is a no-op (no actual split).
        final, msg = _check_cohort_mode_chunking_compat(
            cli_cohort_mode="arrow-streaming",
            resolved_cohort_mode="arrow-streaming",
            chunk_size_mb=70,
            chromosomes=["22"],
            build="GRCh38",
            chr_length_mb=70,
        )
        self.assertEqual(final, "arrow-streaming")
        self.assertIsNone(msg)

    def test_explicit_streaming_with_splitting_chunk_returns_error(self):
        from syntheticgen.cli import _check_cohort_mode_chunking_compat
        # User explicitly asked for --cohort-mode arrow-streaming
        # AND --chr-chunk-mb 5 against a 70 Mb chrom: would split.
        # Helper returns the original mode + ERROR message; the
        # caller is responsible for sys.exiting the body.
        final, msg = _check_cohort_mode_chunking_compat(
            cli_cohort_mode="arrow-streaming",
            resolved_cohort_mode="arrow-streaming",
            chunk_size_mb=5,
            chromosomes=["22"],
            build="GRCh38",
            chr_length_mb=70,
        )
        self.assertEqual(final, "arrow-streaming")
        self.assertIsNotNone(msg)
        self.assertTrue(
            msg.startswith("ERROR: "), f"got: {msg!r}",
        )
        self.assertIn(
            "arrow-streaming does not yet support", msg,
        )
        self.assertIn("chunk_size_mb=5.00", msg)
        self.assertIn("70.0 Mb", msg)
        # PR #58 review: the message must not recommend
        # ``--chr-chunk-mb 0`` since that's the auto-pick sentinel,
        # not a no-chunking flag. It must instead point the user at
        # an explicit chunk >= eff_len_mb (70 here, rounded to 70)
        # or a mode switch.
        self.assertNotIn("--chr-chunk-mb 0", msg)
        self.assertIn("--chr-chunk-mb 70", msg)
        self.assertIn("--cohort-mode arrow", msg)

    def test_auto_streaming_with_splitting_chunk_demotes_to_arrow(self):
        from syntheticgen.cli import _check_cohort_mode_chunking_compat
        # User picked --cohort-mode auto and the resolver settled on
        # arrow-streaming. The auto-picked chunk would now split, so
        # demote to arrow (which preserves mmap-share via the
        # chunked materialised path) and emit INFO so the user can
        # see what happened.
        final, msg = _check_cohort_mode_chunking_compat(
            cli_cohort_mode="auto",
            resolved_cohort_mode="arrow-streaming",
            chunk_size_mb=5,
            chromosomes=["22"],
            build="GRCh38",
            chr_length_mb=70,
        )
        self.assertEqual(final, "arrow")
        self.assertIsNotNone(msg)
        self.assertTrue(
            msg.startswith("INFO: "), f"got: {msg!r}",
        )
        self.assertIn("demoted arrow-streaming → arrow", msg)
        self.assertIn("chunk_size_mb=5.00", msg)

    def test_explicit_streaming_full_contig_uses_effective_length(self):
        # When the user passes --chr-length-mb 0 (full contig) the
        # helper must resolve the length via chr1's contig size
        # (~249 Mb on GRCh38) — otherwise small chunks would look
        # like "no split" and the ERROR wouldn't fire. Regression
        # test for the Comment 1 wiring at the actual check call.
        from syntheticgen.cli import _check_cohort_mode_chunking_compat
        final, msg = _check_cohort_mode_chunking_compat(
            cli_cohort_mode="arrow-streaming",
            resolved_cohort_mode="arrow-streaming",
            chunk_size_mb=50,
            chromosomes=["1"],
            build="GRCh38",
            chr_length_mb=0,
        )
        self.assertEqual(final, "arrow-streaming")
        self.assertTrue(msg.startswith("ERROR: "))
        # The message should report a real per-chrom length, not 0.
        self.assertNotIn("0.0 Mb", msg)

    def test_suggested_chunk_is_ceil_not_round_for_non_integer_length(self):
        # PR #60 review: the suggested chunk in the ERROR message was
        # ``{eff_len_mb:.0f}``, which rounds (Python's banker's rule
        # at .5) and can produce a value BELOW eff_len — so the
        # advised remedy would still split the chrom. Use a pinned
        # explicit chr_length_mb=70.1 to make the bug visible: the
        # suggestion must be 71 (ceil), never 70.
        from syntheticgen.cli import _check_cohort_mode_chunking_compat
        final, msg = _check_cohort_mode_chunking_compat(
            cli_cohort_mode="arrow-streaming",
            resolved_cohort_mode="arrow-streaming",
            chunk_size_mb=5,
            chromosomes=["22"],
            build="GRCh38",
            chr_length_mb=70.1,
        )
        self.assertEqual(final, "arrow-streaming")
        self.assertTrue(msg.startswith("ERROR: "))
        self.assertIn("70.1 Mb", msg)
        self.assertIn("--chr-chunk-mb 71", msg)
        # And critically: must NOT recommend 70 (would still split).
        self.assertNotIn("--chr-chunk-mb 70 ", msg)

    def test_suggested_chunk_never_zero_for_sub_mb_length(self):
        # The original ``.0f`` formatter would have rounded a 0.4 Mb
        # effective length down to 0 and reintroduced the forbidden
        # ``--chr-chunk-mb 0`` recommendation that PR #58 already
        # removed. Pin a sub-Mb explicit length so the helper has to
        # ceil to 1 here.
        from syntheticgen.cli import _check_cohort_mode_chunking_compat
        final, msg = _check_cohort_mode_chunking_compat(
            cli_cohort_mode="arrow-streaming",
            resolved_cohort_mode="arrow-streaming",
            chunk_size_mb=0.1,
            chromosomes=["22"],
            build="GRCh38",
            chr_length_mb=0.4,
        )
        self.assertEqual(final, "arrow-streaming")
        self.assertTrue(msg.startswith("ERROR: "))
        self.assertIn("--chr-chunk-mb 1", msg)
        self.assertNotIn("--chr-chunk-mb 0", msg)

    def test_no_chrom_info_defers_silently(self):
        # Defensive: when chromosomes/build are missing and the user
        # didn't pin chr_length_mb, the effective length resolves to
        # 0 and the helper defers (silent no-op) rather than firing
        # spurious messages off stale data.
        from syntheticgen.cli import _check_cohort_mode_chunking_compat
        final, msg = _check_cohort_mode_chunking_compat(
            cli_cohort_mode="auto",
            resolved_cohort_mode="arrow-streaming",
            chunk_size_mb=5,
            chromosomes=None,
            build=None,
            chr_length_mb=0,
        )
        self.assertEqual(final, "arrow-streaming")
        self.assertIsNone(msg)


if __name__ == "__main__":
    unittest.main()
