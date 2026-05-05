"""Test the BCFTOOLS_STATS module's symbolic-ALT / missing-ALT fallback.

`bcftools stats -s -` aborts with "Requested allele outside valid range"
on records where GT references an allele index that doesn't exist in
ALT — synthetic_people occasionally emits such records (GT=1|0 against
ALT="." after the rsID overlay), and structural-variant records
(`<DEL>` / `<DUP>` / `<INV>`) trip the same code path on some bcftools
versions.

The BCFTOOLS_STATS module's script block tries the full file first;
on failure it retries against ``bcftools view -e 'ALT="." || ALT~"<"'``.
This test verifies the two halves of that contract:

1. The malformed-record class really does break ``bcftools stats``
   (otherwise the fallback would never trigger and the test couldn't
   guard against future regressions).
2. The same filter that the module applies recovers a non-empty
   per-sample stats block for downstream MULTIQC consumption.
"""

import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fixtures import require_tools


# Minimal VCF carrying one good biallelic SNV plus one record where GT
# references allele index 1 but ALT is "." — exactly the synthetic_people
# rsID-overlay output that triggered the original report.
_MALFORMED_VCF = """\
##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=22,length=50818468,assembly=GRCh38>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
22\t1000\trs1\tA\tG\t100\tPASS\tAC=1;AN=4;AF=0.25\tGT:DP\t0|1:30\t0|0:30
22\t2000\trs2\tG\t.\t100\tPASS\tAC=1;AN=2;AF=0.5\tGT:DP\t1|0:25\t0|0:30
"""


class BcftoolsStatsFallbackTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        require_tools("bcftools", "tabix", "bgzip")
        cls.tmpdir = tempfile.mkdtemp(prefix="bcft_stats_fb_")
        plain = os.path.join(cls.tmpdir, "malformed.vcf")
        with open(plain, "w") as fh:
            fh.write(_MALFORMED_VCF)
        cls.vcf = plain + ".gz"
        # bgzip + tabix the fixture so bcftools can index-query it.
        subprocess.run(["bgzip", plain], check=True)
        subprocess.run(["tabix", "-p", "vcf", cls.vcf], check=True)

    def test_full_file_stats_fails(self):
        # Guard against future bcftools versions silently fixing this.
        # If bcftools ever stops complaining about the malformed record
        # the fallback becomes unnecessary — at which point we'd want to
        # know.
        proc = subprocess.run(
            ["bcftools", "stats", "-s", "-", self.vcf],
            capture_output=True, text=True,
        )
        self.assertNotEqual(
            proc.returncode, 0,
            msg="bcftools stats unexpectedly succeeded on a record with "
                "ALT=. and GT=1|0 — the BCFTOOLS_STATS fallback in "
                "modules/bcftools_stats.nf is no longer necessary and "
                "the test should be adapted accordingly.",
        )
        self.assertIn(
            "outside valid range", proc.stderr,
            msg=f"unexpected bcftools error: {proc.stderr[:300]}",
        )

    def test_filtered_pipeline_succeeds(self):
        # Same filter used in modules/bcftools_stats.nf.
        view = subprocess.Popen(
            ["bcftools", "view", "-e", 'ALT="." || ALT~"<"', self.vcf],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        stats = subprocess.run(
            ["bcftools", "stats", "-s", "-"],
            stdin=view.stdout, capture_output=True, text=True,
        )
        view.stdout.close()
        view.wait()
        self.assertEqual(stats.returncode, 0, msg=stats.stderr)
        # SN section should still report the one good biallelic SNV that
        # survived the filter.
        sn_lines = [l for l in stats.stdout.splitlines() if l.startswith("SN")]
        self.assertTrue(any("number of records:\t1" in l for l in sn_lines),
                        msg=f"unexpected SN section: {sn_lines}")
        self.assertTrue(any("number of SNPs:\t1" in l for l in sn_lines),
                        msg=f"unexpected SN section: {sn_lines}")


if __name__ == "__main__":
    unittest.main()
