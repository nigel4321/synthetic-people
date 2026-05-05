"""End-to-end: invoke `nextflow run main.nf` against synthetic VCFs.

Generates a small chr15 + chr22 cohort, runs the pipeline, and verifies the
three published reports + carriers.tsv look right. Skipped automatically if
nextflow is not on PATH.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fixtures import (
    PIPELINE_DIR,
    default_cohort,
    in_range_variants,
    require_tools,
    standard_filename,
)
from synthetic_vcf import Variant, write_vcf


class PipelineE2ETest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        require_tools("nextflow", "bcftools", "tabix", "bgzip")
        cls.tmpdir = tempfile.mkdtemp(prefix="pipeline_e2e_")
        cls.input_dir = os.path.join(cls.tmpdir, "vcfs")
        os.makedirs(cls.input_dir)

        # chr15: target variant in range.
        write_vcf(
            standard_filename(cls.input_dir, "15"),
            "15", in_range_variants(), default_cohort(),
        )
        # chr22: does not contain chr15 → should classify as not_applicable.
        write_vcf(
            standard_filename(cls.input_dir, "22"),
            "22",
            [Variant(pos=16050075, ref="A", alt="G",
                     af_by_pop={"ALL": 0.05})],
            default_cohort(),
        )

        cls.outdir = os.path.join(cls.tmpdir, "results")
        cls.workdir = os.path.join(cls.tmpdir, "work")

    @classmethod
    def tearDownClass(cls):
        # Leave the work dir behind if NF_KEEP is set, so a failing run can
        # be inspected. Otherwise clean up — it is tens of MB per run.
        if not os.environ.get("NF_KEEP"):
            shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_pipeline_runs_and_produces_reports(self):
        glob = os.path.join(self.input_dir, "*.vcf.gz")
        cmd = [
            "nextflow", "run", "main.nf",
            "--input", glob,
            "--outdir", self.outdir,
            "-work-dir", self.workdir,
            "--variant_name", "rs12913832",
            "--variant_chrom", "15",
            "--variant_pos", "28365618",
            "--variant_ref", "A",
            "--variant_alt", "G",
            "--variant_min_af", "0.05",
            "--variant_max_af", "1.0",
            "-ansi-log", "false",
        ]
        proc = subprocess.run(cmd, cwd=PIPELINE_DIR,
                              capture_output=True, text=True)
        self.assertEqual(
            proc.returncode, 0,
            f"nextflow failed\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}",
        )

        # QC + three downstream reports + carriers.tsv should be published.
        for name in ("qc_report.md", "metadata_report.md", "variant_report.md",
                     "carriers_report.md", "carriers.tsv"):
            self.assertTrue(
                os.path.isfile(os.path.join(self.outdir, name)),
                f"missing output: {name}",
            )

        # MultiQC report is best-effort: only assert it lands when
        # multiqc is importable on this host. Some CI runners (e.g. the
        # minimal stdlib-only one the existing pipeline assumes) won't
        # have it — skip the assertion rather than fail.
        try:
            import multiqc  # noqa: F401
        except ImportError:
            pass
        else:
            html = os.path.join(self.outdir, "multiqc_report.html")
            self.assertTrue(os.path.isfile(html),
                            f"missing MultiQC report: {html}")
            data_dir = os.path.join(self.outdir, "multiqc_data")
            self.assertTrue(os.path.isdir(data_dir),
                            f"missing MultiQC data dir: {data_dir}")
            # Custom-content sidecar should land in the parsed data — both
            # samples should appear under the vcf_qc table.
            with open(os.path.join(data_dir, "multiqc_data.json")) as fh:
                mqc_data = json.load(fh)
            saved = mqc_data.get("report_saved_raw_data", {})
            vcf_qc = saved.get("multiqc_vcf_qc_table", {})
            self.assertEqual(len(vcf_qc), 2,
                             f"expected 2 samples in vcf_qc table, got {vcf_qc}")
            # Bcftools-stats native module should also have ingested both.
            bcft = saved.get("multiqc_bcftools_stats", {})
            self.assertEqual(len(bcft), 2,
                             f"expected 2 samples in bcftools_stats, got {bcft}")

        with open(os.path.join(self.outdir, "qc_report.md")) as fh:
            qc_md = fh.read()
        # Both synthetic files should pass QC cleanly.
        self.assertIn("**Files scanned:** 2", qc_md)
        self.assertIn("**Passed:** 2", qc_md)
        self.assertIn("**Failed:** 0", qc_md)

        with open(os.path.join(self.outdir, "variant_report.md")) as fh:
            variant_md = fh.read()
        # chr15 file should be in-range, chr22 file should be not_applicable.
        self.assertIn("present_in_range", variant_md)
        self.assertIn("not_applicable", variant_md)
        self.assertIn("rs12913832", variant_md)

        with open(os.path.join(self.outdir, "carriers_report.md")) as fh:
            carrier_md = fh.read()
        self.assertIn("Allele-count integrity check", carrier_md)
        self.assertIn("Heterozygotes", carrier_md)

        with open(os.path.join(self.outdir, "carriers.tsv")) as fh:
            carriers_tsv = fh.read().splitlines()
        header = carriers_tsv[0].split("\t")
        self.assertEqual(header[:3], ["file", "sample", "variant_id"])
        data_rows = [l for l in carriers_tsv[1:] if l.strip()]
        # All 4 EUR samples should be homozygous alt → at minimum 4 carriers.
        self.assertGreaterEqual(len(data_rows), 4)
        samples_in_tsv = {l.split("\t")[1] for l in data_rows}
        self.assertTrue(
            {"HG00096", "HG00097", "HG00099", "HG00100"}.issubset(samples_in_tsv),
        )

        with open(os.path.join(self.outdir, "metadata_report.md")) as fh:
            meta_md = fh.read()
        self.assertIn("**Files scanned:** 2", meta_md)
        self.assertIn("phase3_shapeit2_mvncall_integrated_v5b", meta_md)


class PipelineStrictQcAbortsTest(unittest.TestCase):
    """Strict QC mode (default) aborts the workflow on bad input."""

    @classmethod
    def setUpClass(cls):
        require_tools("nextflow", "bcftools", "tabix", "bgzip")
        cls.tmpdir = tempfile.mkdtemp(prefix="pipeline_qc_abort_")
        cls.input_dir = os.path.join(cls.tmpdir, "vcfs")
        os.makedirs(cls.input_dir)

        vcf = write_vcf(
            standard_filename(cls.input_dir, "15"),
            "15", in_range_variants(), default_cohort(),
        )
        # Sabotage the index so QC hard-fails.
        os.remove(vcf + ".tbi")
        # Re-create a zero-byte .tbi so the workflow's own pre-check in main.nf
        # (`if (!tbi.exists()) error ...`) doesn't short-circuit ahead of us —
        # we want QC_VALIDATE to be the thing that fails the pipeline.
        open(vcf + ".tbi", "w").close()

        cls.outdir = os.path.join(cls.tmpdir, "results")
        cls.workdir = os.path.join(cls.tmpdir, "work")

    @classmethod
    def tearDownClass(cls):
        if not os.environ.get("NF_KEEP"):
            shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_strict_mode_aborts_on_hard_qc_failure(self):
        glob = os.path.join(self.input_dir, "*.vcf.gz")
        cmd = [
            "nextflow", "run", "main.nf",
            "--input", glob,
            "--outdir", self.outdir,
            "-work-dir", self.workdir,
            "--variant_name", "rs12913832",
            "--variant_chrom", "15", "--variant_pos", "28365618",
            "--variant_ref", "A", "--variant_alt", "G",
            "-ansi-log", "false",
        ]
        proc = subprocess.run(cmd, cwd=PIPELINE_DIR,
                              capture_output=True, text=True)
        self.assertNotEqual(
            proc.returncode, 0,
            "expected pipeline to abort on strict QC failure, but it "
            f"succeeded\nstdout:\n{proc.stdout}",
        )
        combined = proc.stdout + proc.stderr
        self.assertIn("QC_VALIDATE", combined)


if __name__ == "__main__":
    unittest.main()
