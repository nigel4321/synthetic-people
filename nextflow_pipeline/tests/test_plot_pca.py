"""Tests for bin/plot_pca.py — PCA happy path, skip paths, error path.

Heavy-dep tests (numpy + sklearn + matplotlib) self-skip when those are
not importable. Skip-path and missing-input tests run on any host with
bcftools / tabix / bgzip on PATH.
"""

import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fixtures import (
    bin_script,
    default_cohort,
    require_tools,
    standard_filename,
)
from synthetic_vcf import Variant, write_vcf


def _heavy_deps_available() -> bool:
    return all(
        importlib.util.find_spec(m) is not None
        for m in ("numpy", "sklearn", "matplotlib")
    )


def _run_plot_pca(vcf: str, name: str, out_dir: str,
                  min_samples: int = 3, min_variants: int = 10):
    out_png = os.path.join(out_dir, f"{name}.pca.png")
    out_json = os.path.join(out_dir, f"{name}.pca.json")
    cmd = [
        sys.executable, bin_script("plot_pca.py"),
        "--vcf", vcf, "--name", name,
        "--out-png", out_png, "--out-json", out_json,
        "--min-samples", str(min_samples),
        "--min-variants", str(min_variants),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    summary = None
    if os.path.isfile(out_json):
        with open(out_json) as fh:
            summary = json.load(fh)
    png_size = os.path.getsize(out_png) if os.path.isfile(out_png) else None
    return proc, summary, png_size


def _structured_variants(n_variants: int = 60) -> list[Variant]:
    """Build a list of variants where the first half of the cohort sits at
    one allele frequency and the second half sits at another.

    Drives a clear PC1 signal so the happy-path test can assert variance
    explained is non-trivial without depending on a specific seed.
    """
    out: list[Variant] = []
    samples = default_cohort().samples
    half = len(samples) // 2
    for i in range(n_variants):
        # Deterministic per-variant genotypes — population-divergent at
        # every third site, neutral elsewhere. The divergent sites carry
        # the PC1 signal; the neutral ones broaden the spectrum so PCA
        # has full-rank input.
        gts: list[str] = []
        for s_idx in range(len(samples)):
            in_first_half = s_idx < half
            if i % 3 == 0:
                gts.append("1|1" if in_first_half else "0|0")
            elif i % 3 == 1:
                gts.append("0|0" if in_first_half else "1|1")
            else:
                gts.append("0|1")
        out.append(Variant(
            pos=28_000_000 + i * 1000, ref="A", alt="G", variant_id=".",
            genotypes=gts,
        ))
    return out


@unittest.skipUnless(_heavy_deps_available(),
                     "numpy / sklearn / matplotlib not available")
class PlotPcaHappyPathTest(unittest.TestCase):
    """A cohort with a built-in two-group split should yield a real PCA."""

    @classmethod
    def setUpClass(cls):
        require_tools("bcftools", "tabix", "bgzip")
        cls.tmpdir = tempfile.mkdtemp(prefix="plot_pca_happy_")
        cls.vcf = write_vcf(
            standard_filename(cls.tmpdir, "15"), "15",
            _structured_variants(), default_cohort(),
        )

    def test_writes_non_empty_png(self):
        proc, _, png_size = _run_plot_pca(self.vcf, "happy", self.tmpdir)
        self.assertEqual(proc.returncode, 0, msg=proc.stderr)
        # PNG file format is at least ~150 bytes for a 1×1 image. Our 6×5
        # plot will be 1 KB+. Anything smaller indicates the plotting path
        # never ran.
        self.assertIsNotNone(png_size)
        self.assertGreater(png_size, 1000)

    def test_summary_shape(self):
        _, summary, _ = _run_plot_pca(self.vcf, "happy", self.tmpdir)
        self.assertNotIn("skipped", summary)
        self.assertEqual(summary["name"], "happy")
        self.assertEqual(summary["n_samples"], 20)
        self.assertGreater(summary["n_variants_used"], 10)
        self.assertEqual(len(summary["explained_variance_pct"]), 2)
        self.assertEqual(len(summary["samples"]), 20)
        for entry in summary["samples"]:
            self.assertIn("sample", entry)
            self.assertIn("pc1", entry)
            self.assertIn("pc2", entry)

    def test_pc1_separates_the_two_groups(self):
        _, summary, _ = _run_plot_pca(self.vcf, "happy", self.tmpdir)
        # First 10 samples are group A, last 10 are group B by
        # construction. PC1 should put them on opposite sides of zero.
        first_half_pc1 = [s["pc1"] for s in summary["samples"][:10]]
        second_half_pc1 = [s["pc1"] for s in summary["samples"][10:]]
        # All-positive-vs-all-negative is the strict assertion; but PCA
        # can flip sign so check separation rather than direction.
        a_mean = sum(first_half_pc1) / len(first_half_pc1)
        b_mean = sum(second_half_pc1) / len(second_half_pc1)
        self.assertGreater(
            abs(a_mean - b_mean), 1.0,
            msg="PC1 should separate the two structured groups; "
                f"got group means a={a_mean:.2f}, b={b_mean:.2f}",
        )


class PlotPcaSkipPathsTest(unittest.TestCase):
    """Too-small inputs should skip cleanly (exit 0, JSON marks reason)."""

    @classmethod
    def setUpClass(cls):
        require_tools("bcftools", "tabix", "bgzip")
        cls.tmpdir = tempfile.mkdtemp(prefix="plot_pca_skip_")

    def test_too_few_variants(self):
        # Build a VCF with only 2 variants (well below default min=10).
        vcf = write_vcf(
            standard_filename(self.tmpdir, "22") + ".small.vcf.gz", "22",
            [
                Variant(pos=16050075, ref="A", alt="G",
                        af_by_pop={"ALL": 0.1}),
                Variant(pos=16060075, ref="C", alt="T",
                        af_by_pop={"ALL": 0.2}),
            ],
            default_cohort(),
        )
        proc, summary, _ = _run_plot_pca(vcf, "small_v", self.tmpdir)
        self.assertEqual(proc.returncode, 0, msg=proc.stderr)
        self.assertIn("skipped", summary)
        self.assertIn("too few variants", summary["skipped"])

    def test_too_few_samples(self):
        # A 2-sample cohort can't carry a meaningful 2-component PCA.
        from synthetic_vcf import SyntheticCohort
        small_cohort = SyntheticCohort(samples=["S1", "S2"])
        vcf = write_vcf(
            standard_filename(self.tmpdir, "22") + ".tiny.vcf.gz", "22",
            [
                Variant(pos=p, ref="A", alt="G", variant_id=".",
                        genotypes=["0|0", "1|1"])
                for p in range(16_000_000, 16_000_000 + 30 * 1000, 1000)
            ],
            small_cohort,
        )
        proc, summary, _ = _run_plot_pca(vcf, "tiny", self.tmpdir,
                                         min_samples=3)
        self.assertEqual(proc.returncode, 0, msg=proc.stderr)
        self.assertIn("skipped", summary)
        self.assertIn("too few samples", summary["skipped"])


class PlotPcaErrorPathTest(unittest.TestCase):
    """A genuinely unreadable input should fail fast (non-zero exit)."""

    @classmethod
    def setUpClass(cls):
        require_tools("bcftools", "tabix", "bgzip")
        cls.tmpdir = tempfile.mkdtemp(prefix="plot_pca_err_")

    def test_missing_input_file_exits_nonzero(self):
        missing = os.path.join(self.tmpdir, "does_not_exist.vcf.gz")
        proc, _, _ = _run_plot_pca(missing, "ghost", self.tmpdir)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("missing input", proc.stderr)

    def test_corrupt_vcf_surfaces_bcftools_stderr(self):
        # Build a valid VCF, then truncate it so bcftools fails on read.
        # The error message should carry bcftools's own diagnostic so a
        # user can tell what's wrong without consulting the work dir.
        vcf = write_vcf(
            standard_filename(self.tmpdir, "15") + ".corrupt.vcf.gz",
            "15", _structured_variants(n_variants=20)[:20], default_cohort(),
        )
        with open(vcf, "r+b") as fh:
            fh.truncate(100)  # bgzip footer chopped off → unreadable
        proc, _, _ = _run_plot_pca(vcf, "broken", self.tmpdir)
        self.assertNotEqual(proc.returncode, 0)
        # Single-line plot_pca prefix + bcftools's own stderr should both
        # appear; we don't care exactly which BGZF error bcftools raises
        # so long as it's a real diagnostic and not a phantom failure.
        self.assertIn("[plot_pca]", proc.stderr)
        self.assertIn("bcftools query", proc.stderr)


@unittest.skipUnless(_heavy_deps_available(),
                     "numpy / sklearn / matplotlib not available")
class PlotPcaMaxVariantsTest(unittest.TestCase):
    """``--max-variants`` should cap the matrix without failing the
    process — regression test for the SIGTERM-on-cap bug where bcftools
    was being killed by ``proc.terminate()`` and the resulting -15
    returncode was being raised as a phantom failure.
    """

    @classmethod
    def setUpClass(cls):
        require_tools("bcftools", "tabix", "bgzip")
        cls.tmpdir = tempfile.mkdtemp(prefix="plot_pca_cap_")
        cls.vcf = write_vcf(
            standard_filename(cls.tmpdir, "15"), "15",
            _structured_variants(n_variants=80), default_cohort(),
        )

    def test_cap_below_available_succeeds(self):
        # Cap at 25 even though 80 are available. Pre-fix this exited 1
        # with `bcftools query failed (exit -15)`.
        out_png = os.path.join(self.tmpdir, "capped.pca.png")
        out_json = os.path.join(self.tmpdir, "capped.pca.json")
        proc = subprocess.run(
            [sys.executable, bin_script("plot_pca.py"),
             "--vcf", self.vcf, "--name", "capped",
             "--out-png", out_png, "--out-json", out_json,
             "--max-variants", "25"],
            capture_output=True, text=True,
        )
        self.assertEqual(proc.returncode, 0, msg=proc.stderr)
        with open(out_json) as fh:
            summary = json.load(fh)
        self.assertNotIn("skipped", summary)
        # Used variants should be at most the cap (some may be pruned
        # for zero variance after mean imputation).
        self.assertLessEqual(summary["n_variants_used"], 25)
        self.assertGreater(summary["n_variants_used"], 5)


if __name__ == "__main__":
    unittest.main()
