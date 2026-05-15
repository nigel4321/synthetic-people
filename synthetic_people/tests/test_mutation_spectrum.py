"""Tests for the 96-channel mutation spectrum diagnostic (Tier 2 #5).

Covers:

- Channel indexing — every (sub, left, right) tuple maps to exactly
  one index in [0, 96), and the round-trip via ``channel_label`` is
  the canonical ``"<left>[<ref>><alt>]<right>"`` form.
- Pyrimidine normalisation — purine-context calls (G>A, etc.) fold
  to their pyrimidine equivalent via reverse-complement, including
  reversal of the flanking bases.
- N / IUPAC handling — out-of-alphabet flanks land in ``n_excluded``
  rather than poisoning a real channel.
- End-to-end on a hand-built FASTA + BCF: the spectrum's channel
  counts match what hand-counting predicts.
"""

from __future__ import annotations

import importlib.util
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from syntheticgen.mutation_spectrum import (  # noqa: E402
    MutationSpectrum,
    N_CHANNELS,
    all_channel_labels,
    channel_index,
    channel_label,
    compute_spectrum,
    pyrimidine_normalize,
    to_jsonable,
)

_HAVE_PYSAM = importlib.util.find_spec("pysam") is not None
_HAVE_BCFTOOLS = shutil.which("bcftools") is not None
# ``_write_vcf`` shells out to bgzip + tabix in addition to bcftools
# downstream of compute_spectrum. Gate the end-to-end tests on all
# four so a host with bcftools-only (or pysam-only) skips cleanly
# instead of hitting CalledProcessError mid-fixture. Same gating
# pattern as ``tests/test_cli_modes.py``.
_HAVE_BGZIP = shutil.which("bgzip") is not None
_HAVE_TABIX = shutil.which("tabix") is not None


# ---------------------------------------------------------------------------
# Pure-Python helpers — no pysam/bcftools needed
# ---------------------------------------------------------------------------


class ChannelIndexTest(unittest.TestCase):
    """Channel indexing is the foundation of the spectrum's
    interpretability. Lock the index ↔ label mapping so downstream
    signature-comparison tooling can rely on the order being stable.
    """

    def test_first_channel_is_a_c_to_a_a(self):
        # Index 0 = C>A substitution, left=A, right=A → "A[C>A]A".
        self.assertEqual(channel_index("C", "A", "A", "A"), 0)
        self.assertEqual(channel_label(0), "A[C>A]A")

    def test_last_channel_is_t_t_to_g_t(self):
        # Index 95 = T>G substitution, left=T, right=T → "T[T>G]T".
        self.assertEqual(channel_index("T", "G", "T", "T"), 95)
        self.assertEqual(channel_label(95), "T[T>G]T")

    def test_substitution_block_size(self):
        # Each substitution block is 16 channels (4 left × 4 right).
        # C>A spans [0, 16); C>G spans [16, 32); etc.
        self.assertEqual(channel_index("C", "G", "A", "A"), 16)
        self.assertEqual(channel_index("C", "T", "A", "A"), 32)
        self.assertEqual(channel_index("T", "A", "A", "A"), 48)
        self.assertEqual(channel_index("T", "C", "A", "A"), 64)
        self.assertEqual(channel_index("T", "G", "A", "A"), 80)

    def test_every_index_round_trips_via_label(self):
        # all_channel_labels() must produce 96 unique labels in
        # canonical order. Lock the count + uniqueness as a regression
        # guard against an accidental dedup or off-by-one in the
        # labelling helpers.
        labels = all_channel_labels()
        self.assertEqual(len(labels), N_CHANNELS)
        self.assertEqual(len(set(labels)), N_CHANNELS)

    def test_invalid_input_raises(self):
        # Purine REF must not reach channel_index — pyrimidine_normalize
        # is the gate. A bug that bypasses normalisation should raise,
        # not silently pick the wrong channel.
        with self.assertRaises(ValueError):
            channel_index("G", "A", "A", "A")
        # Out-of-alphabet flank:
        with self.assertRaises(ValueError):
            channel_index("C", "A", "N", "A")

    def test_label_out_of_range_raises(self):
        with self.assertRaises(ValueError):
            channel_label(-1)
        with self.assertRaises(ValueError):
            channel_label(N_CHANNELS)


class PyrimidineNormalizeTest(unittest.TestCase):
    """Folding purine-context calls to pyrimidine-context is what makes
    the spectrum a 96-channel space (not a 192-channel one). The
    reverse-complement must also reverse the flank order — a G>A in
    context A[G>A]C is the same mutation as G[C>T]T on the opposite
    strand, NOT T[C>T]G (which would be the result of complementing
    each base in place without swapping flank sides).
    """

    def test_pyrimidine_context_passthrough(self):
        # C-context: no flip.
        self.assertEqual(
            pyrimidine_normalize("C", "T", "A", "G"),
            ("C", "T", "A", "G"),
        )
        # T-context: no flip.
        self.assertEqual(
            pyrimidine_normalize("T", "C", "G", "A"),
            ("T", "C", "G", "A"),
        )

    def test_purine_context_reverse_complement_with_flank_reversal(self):
        # A[G>A]C on the + strand = G[C>T]T on the – strand.
        # - ref G → C, alt A → T
        # - old left A → complement T, BUT it becomes the new RIGHT.
        # - old right C → complement G, BUT it becomes the new LEFT.
        self.assertEqual(
            pyrimidine_normalize("G", "A", "A", "C"),
            ("C", "T", "G", "T"),
        )

    def test_a_to_t_purine_flips_to_t_to_a(self):
        # A[A>T]G on + = C[T>A]T on –.
        self.assertEqual(
            pyrimidine_normalize("A", "T", "A", "G"),
            ("T", "A", "C", "T"),
        )

    def test_n_flank_returns_none(self):
        self.assertIsNone(pyrimidine_normalize("C", "T", "N", "G"))
        self.assertIsNone(pyrimidine_normalize("C", "T", "A", "N"))

    def test_iupac_alt_returns_none(self):
        # An ALT outside ACGT should land in the exclusion bucket
        # rather than poison a channel.
        self.assertIsNone(pyrimidine_normalize("C", "Y", "A", "G"))


class MutationSpectrumDataclassTest(unittest.TestCase):
    """The dataclass is light, but the ``add`` accumulator and the
    ``fractions`` zero-vector handling are both used by the cohort
    aggregation path — lock both."""

    def test_add_accumulates_counts_and_totals(self):
        a = MutationSpectrum(
            counts=[0] * N_CHANNELS, n_total=10, n_excluded=2,
        )
        a.counts[0] = 5
        a.counts[10] = 3
        b = MutationSpectrum(
            counts=[0] * N_CHANNELS, n_total=4, n_excluded=1,
        )
        b.counts[0] = 1
        b.counts[20] = 2
        a.add(b)
        self.assertEqual(a.counts[0], 6)
        self.assertEqual(a.counts[10], 3)
        self.assertEqual(a.counts[20], 2)
        self.assertEqual(a.n_total, 14)
        self.assertEqual(a.n_excluded, 3)

    def test_fractions_zero_vector_when_no_binned_snvs(self):
        # An all-zero counts vector must not divide by zero.
        s = MutationSpectrum()
        self.assertEqual(s.fractions(), [0.0] * N_CHANNELS)

    def test_fractions_sum_to_one_when_any_binned(self):
        s = MutationSpectrum()
        s.counts[0] = 1
        s.counts[1] = 3
        fr = s.fractions()
        self.assertAlmostEqual(sum(fr), 1.0)
        self.assertAlmostEqual(fr[0], 0.25)
        self.assertAlmostEqual(fr[1], 0.75)


class ToJsonableTest(unittest.TestCase):
    def test_schema_shape(self):
        s = MutationSpectrum()
        s.n_total = 10
        s.n_excluded = 2
        s.counts[0] = 5
        s.counts[1] = 3
        out = to_jsonable(s)
        self.assertEqual(out["n_total"], 10)
        self.assertEqual(out["n_excluded"], 2)
        self.assertEqual(out["n_binned"], 8)
        self.assertEqual(len(out["channels"]), N_CHANNELS)
        first = out["channels"][0]
        self.assertEqual(first["label"], "A[C>A]A")
        self.assertEqual(first["count"], 5)
        self.assertAlmostEqual(first["fraction"], 5 / 8)


# ---------------------------------------------------------------------------
# End-to-end: compute_spectrum reads a real BCF + FASTA
# ---------------------------------------------------------------------------


def _write_fasta(dir_: Path, records: dict) -> Path:
    import pysam
    dir_.mkdir(parents=True, exist_ok=True)
    path = dir_ / "tiny.fa"
    path.write_text(
        "\n".join(f">{name}\n{seq}" for name, seq in records.items())
        + "\n"
    )
    pysam.faidx(str(path))
    return path


def _write_vcf(dir_: Path, chrom: str, records: list[tuple]) -> Path:
    """Write a tiny VCF with the given (pos, ref, alt) records, then
    bgzip + tabix-index it. Records must be in genomic order."""
    vcf_path = dir_ / "tiny.vcf"
    header = "\n".join([
        "##fileformat=VCFv4.2",
        f"##contig=<ID={chrom}>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ])
    body = "\n".join(
        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t."
        for pos, ref, alt in records
    )
    vcf_path.write_text(header + "\n" + body + "\n")
    subprocess.run(["bgzip", "-f", str(vcf_path)], check=True)
    subprocess.run(
        ["tabix", "-p", "vcf", str(vcf_path) + ".gz"], check=True,
    )
    return Path(str(vcf_path) + ".gz")


@unittest.skipUnless(
    _HAVE_PYSAM and _HAVE_BCFTOOLS and _HAVE_BGZIP and _HAVE_TABIX,
    "pysam + bcftools + bgzip + tabix needed for end-to-end test",
)
class ComputeSpectrumEndToEndTest(unittest.TestCase):
    """End-to-end: build a tiny FASTA + VCF, hand-compute the expected
    channel counts, run ``compute_spectrum``, and assert."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        # FASTA: 1-indexed positions 1..16 are A C G T A C G T C C G C T A T G.
        self.fasta_seq = "ACGTACGTCCGCTATG"
        self.fasta = _write_fasta(self.dir, {"22": self.fasta_seq})

    def tearDown(self):
        self.tmp.cleanup()

    def test_simple_c_to_t_at_acg_context(self):
        # POS=2, REF=C, ALT=T, context (1,2,3) = ACG → channel A[C>T]G.
        from syntheticgen.reference import load_fasta
        vcf = _write_vcf(self.dir, "22", [(2, "C", "T")])
        fa = load_fasta(self.fasta)
        try:
            spec = compute_spectrum(vcf, fa)
        finally:
            fa.close()
        self.assertEqual(spec.n_total, 1)
        self.assertEqual(spec.n_excluded, 0)
        idx = channel_index("C", "T", "A", "G")
        self.assertEqual(spec.counts[idx], 1)
        # Every other channel is zero.
        self.assertEqual(sum(spec.counts), 1)

    def test_purine_ref_folds_to_pyrimidine(self):
        # POS=3, REF=G, ALT=A. Context (2,3,4) = C G T.
        # G>A folds to C>T via RC; old left C → complement G (new
        # right), old right T → complement A (new left). So channel
        # is A[C>T]G.
        from syntheticgen.reference import load_fasta
        vcf = _write_vcf(self.dir, "22", [(3, "G", "A")])
        fa = load_fasta(self.fasta)
        try:
            spec = compute_spectrum(vcf, fa)
        finally:
            fa.close()
        idx = channel_index("C", "T", "A", "G")
        self.assertEqual(spec.counts[idx], 1)

    def test_pos_1_left_flank_off_chrom_end_excludes(self):
        # POS=1 has no left flank (pos 0 is off-end). fetch_ref_base
        # returns "N" there → pyrimidine_normalize returns None →
        # record lands in n_excluded.
        from syntheticgen.reference import load_fasta
        vcf = _write_vcf(self.dir, "22", [(1, "A", "C")])
        fa = load_fasta(self.fasta)
        try:
            spec = compute_spectrum(vcf, fa)
        finally:
            fa.close()
        self.assertEqual(spec.n_total, 1)
        self.assertEqual(spec.n_excluded, 1)
        self.assertEqual(sum(spec.counts), 0)

    def test_indels_excluded_at_iter_layer(self):
        # bcftools query's TYPE="snp" && N_ALT=1 filter drops indels
        # and multi-allelics at the source — they never reach the
        # channel-binning logic, so n_total counts only SNVs.
        from syntheticgen.reference import load_fasta
        vcf = _write_vcf(self.dir, "22", [
            (2, "C", "T"),     # SNV
            (5, "A", "ACG"),   # insertion
            (10, "CG", "C"),   # deletion
        ])
        fa = load_fasta(self.fasta)
        try:
            spec = compute_spectrum(vcf, fa)
        finally:
            fa.close()
        self.assertEqual(spec.n_total, 1)
        self.assertEqual(spec.n_excluded, 0)

    def test_bcftools_failure_raises_with_stderr_tail(self):
        # PR #92 review (Copilot): a corrupt / missing input must
        # surface as a RuntimeError with the bcftools stderr tail,
        # not silently produce an empty spectrum. Point at a path
        # that doesn't exist — bcftools query exits non-zero with
        # a clear error message.
        from syntheticgen.reference import load_fasta
        # Build a valid FASTA handle so the failure can only come
        # from the bcftools subprocess, not from fasta loading.
        fa = load_fasta(self.fasta)
        try:
            missing = self.dir / "does_not_exist.vcf.gz"
            with self.assertRaises(RuntimeError) as ctx:
                compute_spectrum(missing, fa)
            msg = str(ctx.exception)
            self.assertIn("bcftools query exited", msg)
            self.assertIn(str(missing), msg)
        finally:
            fa.close()


if __name__ == "__main__":
    unittest.main()
