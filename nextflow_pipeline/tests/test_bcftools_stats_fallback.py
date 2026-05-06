"""Test the BCFTOOLS_STATS module's missing-ALT fallback.

`bcftools stats -s -` aborts with "Requested allele outside valid range"
on records where GT references an allele index that doesn't exist in
ALT — synthetic_people occasionally emits such records (GT=1|0 against
ALT="." after the rsID overlay rewrites a cohort site to a dbSNP
reference-only entry while leaving the cohort GT block intact). Such
records are malformed by VCF spec and contribute zero information to
any stats; dropping them is loss-free.

The BCFTOOLS_STATS module's script block tries the full file first;
on failure it retries against ``bcftools view -e 'ALT="."'``. This
test verifies three pieces of the contract:

1. The malformed-record class really does break ``bcftools stats``
   (otherwise the fallback is dead code).
2. The same filter the module applies recovers a non-empty per-sample
   stats block for downstream MULTIQC consumption.
3. Symbolic-ALT structural variants (`<DEL>` / `<DUP>` / `<INV>`) are
   *not* dropped by the filter — they're loss-free to bcftools stats
   on their own, so keeping them preserves the SV count in the
   bcftools-stats MultiQC panel.
"""

import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fixtures import require_tools


# VCF carrying:
#   - one good biallelic SNV
#   - one structural-variant record (<DEL>) — should survive the filter,
#     contributes to the stats "others" count
#   - one malformed record: GT=1|0 against ALT="." — exactly the
#     synthetic_people rsID-overlay output that triggered the original
#     report; this is what the filter must drop
_MALFORMED_VCF = """\
##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=22,length=50818468,assembly=GRCh38>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=END,Number=1,Type=Integer,Description="SV end position">
##ALT=<ID=DEL,Description="Deletion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
22\t1000\trs1\tA\tG\t100\tPASS\tAC=1;AN=4;AF=0.25\tGT:DP\t0|1:30\t0|0:30
22\t2000\tsv1\tA\t<DEL>\t100\tPASS\tAC=1;AN=4;AF=0.25;SVTYPE=DEL;SVLEN=-500;END=2500\tGT:DP\t0|1:30\t0|0:30
22\t3000\trs2\tG\t.\t100\tPASS\tAC=1;AN=2;AF=0.5\tGT:DP\t1|0:25\t0|0:30
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

    def test_filtered_pipeline_succeeds_and_keeps_svs(self):
        # Same filter used in modules/bcftools_stats.nf — drops *only*
        # missing-ALT records, keeps symbolic ALTs in the stats input
        # so SVs are reflected in the "others" count.
        view = subprocess.Popen(
            ["bcftools", "view", "-e", 'ALT="."', self.vcf],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        stats = subprocess.run(
            ["bcftools", "stats", "-s", "-"],
            stdin=view.stdout, capture_output=True, text=True,
        )
        view.stdout.close()
        view.wait()
        self.assertEqual(stats.returncode, 0, msg=stats.stderr)
        sn_lines = [l for l in stats.stdout.splitlines() if l.startswith("SN")]
        # Filter drops the one malformed record (3 → 2). The biallelic
        # SNV stays as 1 SNP; the <DEL> stays in the "others" count.
        self.assertTrue(any("number of records:\t2" in l for l in sn_lines),
                        msg=f"unexpected SN section: {sn_lines}")
        self.assertTrue(any("number of SNPs:\t1" in l for l in sn_lines),
                        msg=f"unexpected SN section: {sn_lines}")
        self.assertTrue(any("number of others:\t1" in l for l in sn_lines),
                        msg=f"SVs were dropped by the filter — bcftools "
                            f"stats panel will underreport them. SN: {sn_lines}")


if __name__ == "__main__":
    unittest.main()
