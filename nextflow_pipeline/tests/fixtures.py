"""Shared test helpers: paths to bin/ scripts and standard synthetic cohorts."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import unittest

from synthetic_vcf import (
    SyntheticCohort,
    Variant,
    standard_filename,
    write_vcf,
)

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
PIPELINE_DIR = os.path.dirname(TESTS_DIR)
BIN_DIR = os.path.join(PIPELINE_DIR, "bin")


def bin_script(name: str) -> str:
    return os.path.join(BIN_DIR, name)


def require_tools(*tools: str) -> None:
    """Skip the current test module if a prerequisite binary is missing."""
    missing = [t for t in tools if not shutil.which(t)]
    if missing:
        raise unittest.SkipTest(f"missing required tools: {', '.join(missing)}")


def run_script(script: str, *args: str, cwd: str | None = None,
               check: bool = True) -> subprocess.CompletedProcess:
    """Invoke a pipeline bin/ script directly (bypasses Nextflow)."""
    cmd = [sys.executable, script, *args]
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if check and proc.returncode != 0:
        raise AssertionError(
            f"{script} failed (exit {proc.returncode})\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return proc


# -----------------------------------------------------------------------------
# Standard cohort + variants used across tests
# -----------------------------------------------------------------------------

# Target variant = rs12913832 (HERC2 eye-colour) at chr15:28365618 A>G.
# High EUR AF is biologically realistic and gives us a value firmly above
# variant_min_af=0.05 for the in-range test, plus a low AFR AF we can use
# to exercise the "below threshold" path by raising min_af.
TARGET_VARIANT = {
    "chrom": "15",
    "pos": 28365618,
    "ref": "A",
    "alt": "G",
    "id": "rs12913832",
}


def default_cohort() -> SyntheticCohort:
    return SyntheticCohort()


def bcftools_stats_padding_variants(
    chrom: str,
    n: int = 60,
    start: int = 30_000_000,
    step: int = 5_000,
    seed: int = 0,
) -> list[Variant]:
    """Synthesise ``n`` background variants with a realistic
    transition/transversion mix (~2:1) and varying allele frequencies.

    bcftools stats computes its per-sample Ts/Tv ratio from the PSC
    section; if the test cohort only has a handful of variants the
    per-sample counts collapse to all-zeros, which makes MultiQC's
    bcftools module raise ``ValueError: No datasets to plot`` on
    its Ts/Tv bargraph (multiqc 1.34 doesn't guard against the
    all-zero edge case before calling ``bargraph.plot``). Padding
    the input VCFs with this realistic background gives bcftools
    stats enough non-trivial PSC data for the plot to render.

    None of the pipeline-test assertions depend on these padding
    variants — they only assert on the named test variants
    (rs12913832 etc.) and chromosome-level metadata. Padding is
    purely there to exercise the production MULTIQC stage end-to-
    end with realistic data shape, rather than masking edge cases
    behind ``errorStrategy 'ignore'``.

    Deterministic in ``(chrom, n, start, step, seed)``.
    """
    import random
    rng = random.Random(f"bcftools-padding|{chrom}|{seed}")
    transitions = [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")]
    transversions = [
        ("A", "C"), ("C", "A"), ("A", "T"), ("T", "A"),
        ("G", "C"), ("C", "G"), ("G", "T"), ("T", "G"),
    ]
    out: list[Variant] = []
    for i in range(n):
        # ~2:1 Ti/Tv to match the human germline empirical ratio.
        if rng.random() < 2 / 3:
            ref, alt = rng.choice(transitions)
        else:
            ref, alt = rng.choice(transversions)
        # Spread AFs across the full common-variant range so different
        # samples carry different sites — needed for non-trivial PSC.
        af = rng.uniform(0.05, 0.50)
        out.append(Variant(
            pos=start + i * step,
            ref=ref, alt=alt,
            variant_id=".",
            af_by_pop={"ALL": af},
        ))
    return out


def in_range_variants() -> list[Variant]:
    """chr15 with rs12913832 present at a cohort-level AF around 0.25."""
    return [
        Variant(
            pos=28000000, ref="C", alt="T", variant_id=".",
            af_by_pop={"ALL": 0.05},
        ),
        Variant(
            pos=TARGET_VARIANT["pos"],
            ref=TARGET_VARIANT["ref"],
            alt=TARGET_VARIANT["alt"],
            variant_id=TARGET_VARIANT["id"],
            # 20 samples; EUR all-alt + modest elsewhere → cohort AF ≈ 0.25.
            af_by_pop={
                "EUR": 1.0, "AFR": 0.0, "EAS": 0.125,
                "SAS": 0.25, "AMR": 0.0,
            },
        ),
        Variant(
            pos=28500000, ref="G", alt="A", variant_id=".",
            af_by_pop={"ALL": 0.5},
        ),
    ]


def absent_variants() -> list[Variant]:
    """chr15 where the target position has a different ALT allele."""
    return [
        Variant(
            pos=TARGET_VARIANT["pos"], ref="A", alt="C",
            af_by_pop={"ALL": 0.1},
        ),
    ]


def position_empty_variants() -> list[Variant]:
    """chr15 with no record at the target position."""
    return [
        Variant(pos=27000000, ref="C", alt="T", af_by_pop={"ALL": 0.1}),
        Variant(pos=29000000, ref="G", alt="A", af_by_pop={"ALL": 0.1}),
    ]


def above_threshold_variants() -> list[Variant]:
    """All samples alt → AF=1.0, which will exceed any max_af < 1.0."""
    return [
        Variant(
            pos=TARGET_VARIANT["pos"],
            ref=TARGET_VARIANT["ref"],
            alt=TARGET_VARIANT["alt"],
            variant_id=TARGET_VARIANT["id"],
            genotypes=["1|1"] * len(SyntheticCohort().samples),
        ),
    ]


def below_threshold_variants() -> list[Variant]:
    """Exactly one alt allele in the cohort → AF = 1/40 = 0.025."""
    gts = ["0|0"] * 20
    gts[0] = "0|1"
    return [
        Variant(
            pos=TARGET_VARIANT["pos"],
            ref=TARGET_VARIANT["ref"],
            alt=TARGET_VARIANT["alt"],
            variant_id=TARGET_VARIANT["id"],
            genotypes=gts,
        ),
    ]


def af_unknown_variants() -> list[Variant]:
    """Variant record exists but all genotypes missing + AF stripped from INFO.

    fill-tags will compute AF=0 from the called genotypes (of which there are
    none), so the scan script emits the record with AF=0 — but if we also
    strip AF and make all genotypes missing, fill-tags produces no value and
    the pipeline should land on `present_af_unknown`.
    """
    return [
        Variant(
            pos=TARGET_VARIANT["pos"],
            ref=TARGET_VARIANT["ref"],
            alt=TARGET_VARIANT["alt"],
            variant_id=TARGET_VARIANT["id"],
            all_missing=True,
            strip_af_info=True,
        ),
    ]


def build_standard_cohort_vcf(path: str, chrom: str = "15") -> str:
    """Most-commonly-used fixture: chr15 VCF with rs12913832 in-range."""
    return write_vcf(path, chrom, in_range_variants(), default_cohort())
