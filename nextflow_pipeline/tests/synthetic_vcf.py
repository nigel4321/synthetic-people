"""Generate synthetic VCF files that mimic 1000 Genomes Phase 3 data.

Output is bgzipped + tabix-indexed and matches the conventions the pipeline
expects (VCFv4.1 header, ##reference, ##source=1000GenomesPhase3Pipeline,
per-super-pop AF INFO fields, phased diploid genotypes, filename stamped with
`phase3_shapeit2_mvncall_integrated_v5b.<date>`).

Stdlib only — no pysam/cyvcf2 dependency.
"""

from __future__ import annotations

import os
import random
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Iterable


# Real 1000G sample IDs grouped by super-population. Using real IDs keeps the
# synthetic data recognisable to anyone who works with this dataset.
SUPERPOPULATIONS = {
    "EUR": ["HG00096", "HG00097", "HG00099", "HG00100"],
    "AFR": ["NA18486", "NA18488", "NA18489", "NA18498"],
    "EAS": ["HG00403", "HG00404", "HG00406", "HG00407"],
    "SAS": ["HG03006", "HG03007", "HG03008", "HG03009"],
    "AMR": ["HG01112", "HG01113", "HG01119", "HG01121"],
}

DEFAULT_SAMPLES: list[str] = [s for pop in SUPERPOPULATIONS.values() for s in pop]

GRCH37_CONTIG_LENGTHS = {
    "1":  249250621, "2":  243199373, "3":  198022430, "4":  191154276,
    "5":  180915260, "6":  171115067, "7":  159138663, "8":  146364022,
    "9":  141213431, "10": 135534747, "11": 135006516, "12": 133851895,
    "13": 115169878, "14": 107349540, "15": 102531392, "16":  90354753,
    "17":  81195210, "18":  78077248, "19":  59128983, "20":  63025520,
    "21":  48129895, "22":  51304566, "X": 155270560, "Y":  59373566,
}


@dataclass
class Variant:
    """A single site to emit in the synthetic VCF.

    Either provide `genotypes` directly (one string per sample, e.g. "0|1")
    or provide `af_by_pop` and let the generator draw phased genotypes.
    """

    pos: int
    ref: str
    alt: str
    variant_id: str = "."
    # Direct-specification path: one genotype string per sample.
    genotypes: list[str] | None = None
    # Sampling path: per-super-pop alt-allele frequency. Unseeded callers
    # should construct the generator with a fixed seed for reproducibility.
    af_by_pop: dict[str, float] | None = None
    # If True, omit the AF INFO field entirely (to exercise recomputation /
    # `present_af_unknown` paths in the pipeline). AC/AN are still emitted.
    strip_af_info: bool = False
    # Force all genotypes to "./." — used to exercise AF-unknown edge cases.
    all_missing: bool = False
    extra_info: dict[str, str] = field(default_factory=dict)


@dataclass
class SyntheticCohort:
    samples: list[str] = field(default_factory=lambda: list(DEFAULT_SAMPLES))
    # Which super-pop each sample belongs to. Samples not listed here default
    # to "EUR" for computation convenience.
    pop_by_sample: dict[str, str] = field(default_factory=dict)

    def __post_init__(self):
        if not self.pop_by_sample:
            self.pop_by_sample = {
                s: pop
                for pop, members in SUPERPOPULATIONS.items()
                for s in members
            }

    def samples_in_pop(self, pop: str) -> list[str]:
        return [s for s in self.samples if self.pop_by_sample.get(s) == pop]


def _draw_phased_gt(rng: random.Random, af: float) -> str:
    """Draw a phased diploid genotype given an alt-allele frequency."""
    a = "1" if rng.random() < af else "0"
    b = "1" if rng.random() < af else "0"
    return f"{a}|{b}"


def _realize_genotypes(variant: Variant, cohort: SyntheticCohort,
                       rng: random.Random) -> list[str]:
    if variant.all_missing:
        return ["./." for _ in cohort.samples]
    if variant.genotypes is not None:
        if len(variant.genotypes) != len(cohort.samples):
            raise ValueError(
                f"variant at pos {variant.pos}: got {len(variant.genotypes)} "
                f"genotypes for {len(cohort.samples)} samples"
            )
        return list(variant.genotypes)
    af_by_pop = variant.af_by_pop or {}
    default_af = af_by_pop.get("ALL", 0.0)
    out: list[str] = []
    for s in cohort.samples:
        pop = cohort.pop_by_sample.get(s, "EUR")
        af = af_by_pop.get(pop, default_af)
        out.append(_draw_phased_gt(rng, af))
    return out


def _compute_pop_stats(gts: list[str], samples: list[str],
                       cohort: SyntheticCohort) -> tuple[int, int, dict[str, float]]:
    """Return (AC, AN, per-pop AF)."""
    ac = an = 0
    per_pop_ac: dict[str, int] = {}
    per_pop_an: dict[str, int] = {}
    for s, gt in zip(samples, gts):
        pop = cohort.pop_by_sample.get(s, "EUR")
        for allele in gt.replace("|", "/").split("/"):
            if allele in ("0", "1"):
                an += 1
                per_pop_an[pop] = per_pop_an.get(pop, 0) + 1
                if allele == "1":
                    ac += 1
                    per_pop_ac[pop] = per_pop_ac.get(pop, 0) + 1
    pop_af = {}
    for pop in ("EAS", "AMR", "AFR", "EUR", "SAS"):
        pac = per_pop_ac.get(pop, 0)
        pan = per_pop_an.get(pop, 0)
        pop_af[pop] = (pac / pan) if pan else 0.0
    return ac, an, pop_af


def _format_info(ac: int, an: int, pop_af: dict[str, float],
                 variant: Variant) -> str:
    af_overall = (ac / an) if an else 0.0
    pieces: list[str] = []
    if not variant.strip_af_info:
        pieces.append(f"AC={ac}")
        pieces.append(f"AF={af_overall:.6g}")
    pieces.append(f"AN={an}")
    pieces.append(f"NS={an // 2 if an else 0}")
    if not variant.strip_af_info:
        # Order deliberately matches real 1000G INFO line ordering.
        for pop in ("EAS", "AMR", "AFR", "EUR", "SAS"):
            pieces.append(f"{pop}_AF={pop_af.get(pop, 0.0):.6g}")
    pieces.append("VT=SNP")
    for k, v in variant.extra_info.items():
        pieces.append(f"{k}={v}")
    return ";".join(pieces)


DEFAULT_REFERENCE = (
    "ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/"
    "reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
)


def _format_per_call(gt_str: str, rng: random.Random) -> str:
    """Build the ``GT:DP:GQ:AD`` value for one sample's call.

    The pipeline's test VCFs were originally GT-only, which left
    bcftools stats reporting zero average depth and zero singletons
    per sample. MultiQC 1.28+ now raises ``ValueError("No datasets
    to plot")`` whenever a per-sample plot dataset is all zeros
    (e.g. the Sequencing-depth bargraph), so the e2e test would
    crash inside MultiQC even though the rest of the pipeline ran
    cleanly. Emitting plausible quality metrics here keeps the test
    fixtures within the data shape those downstream tools expect,
    without touching production code.

    Distributions are intentionally light: mean DP ~30 (Gaussian,
    clamped to [10, 60]); GQ 99 for confident hom calls, ~85 for
    heterozygous; AD allocated so the two counts sum to DP and
    proportions reflect the GT (50/50 for het with small noise,
    DP/0 or 0/DP for hom). Drawn from the same seeded RNG that
    produces the genotype so the entire VCF stays deterministic at
    a given seed.

    Missing calls (``./.`` or ``.|.``) emit ``./.:.:.:.,.`` — the
    VCF convention for "no data for this sample".
    """
    if gt_str in ("./.", ".|."):
        return f"{gt_str}:.:.:.,."

    sep = "|" if "|" in gt_str else "/"
    try:
        a1, a2 = gt_str.split(sep)
        a1_int = int(a1)
        a2_int = int(a2)
    except ValueError:
        return f"{gt_str}:.:.:.,."

    is_hom_ref = a1_int == 0 and a2_int == 0
    is_hom_alt = a1_int == a2_int and a1_int != 0
    is_het = a1_int != a2_int

    dp = max(10, min(60, int(rng.gauss(30, 5))))
    if is_het:
        gq = max(30, min(99, int(rng.gauss(85, 8))))
    else:
        gq = max(60, min(99, int(rng.gauss(99, 5))))

    if is_hom_ref:
        ref_d, alt_d = dp, 0
    elif is_hom_alt:
        ref_d, alt_d = 0, dp
    else:
        # Het: ~50/50 with small jitter, then re-clamp to sum to DP
        ref_d = max(0, min(dp, int(rng.gauss(dp / 2, dp * 0.1))))
        alt_d = dp - ref_d

    return f"{gt_str}:{dp}:{gq}:{ref_d},{alt_d}"


def _build_header(
    chrom: str,
    cohort: SyntheticCohort,
    file_date: str,
    reference: str = DEFAULT_REFERENCE,
    contigs_override: dict[str, int] | None = None,
    info_declarations: list[str] | None = None,
    declare_format_gt: bool = True,
) -> list[str]:
    lines = [
        "##fileformat=VCFv4.1",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        f"##fileDate={file_date}",
        f"##reference={reference}",
        "##source=1000GenomesPhase3Pipeline",
    ]
    # Full GRCh37 contig set keeps the header representative and means the
    # pipeline's chr-resolution logic (bare '15' vs 'chr15') is exercised.
    contig_map = contigs_override if contigs_override is not None else GRCH37_CONTIG_LENGTHS
    for c, length in contig_map.items():
        assembly = "b37" if contigs_override is None else "synthetic"
        lines.append(f"##contig=<ID={c},assembly={assembly},length={length}>")
    if info_declarations is None:
        info_declarations = [
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number '
            'of alternate alleles in called genotypes">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated '
            'allele frequency in the range (0,1)">',
            '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number '
            'of alleles in called genotypes">',
            '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of '
            'samples with data">',
            '##INFO=<ID=EAS_AF,Number=A,Type=Float,Description="Allele '
            'frequency in the EAS populations calculated from AC and AN">',
            '##INFO=<ID=EUR_AF,Number=A,Type=Float,Description="Allele '
            'frequency in the EUR populations calculated from AC and AN">',
            '##INFO=<ID=AFR_AF,Number=A,Type=Float,Description="Allele '
            'frequency in the AFR populations calculated from AC and AN">',
            '##INFO=<ID=AMR_AF,Number=A,Type=Float,Description="Allele '
            'frequency in the AMR populations calculated from AC and AN">',
            '##INFO=<ID=SAS_AF,Number=A,Type=Float,Description="Allele '
            'frequency in the SAS populations calculated from AC and AN">',
            '##INFO=<ID=VT,Number=.,Type=String,Description="indicates what '
            'type of variant the line represents">',
        ]
    lines.extend(info_declarations)
    if declare_format_gt:
        lines.extend([
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            # DP/GQ/AD declarations support the GT:DP:GQ:AD per-sample
            # block emitted by write_vcf. Real Phase 3 VCFs use just GT,
            # but downstream tools (bcftools stats Sequencing-depth
            # bargraph; MultiQC's bcftools module from 1.28+) require
            # non-trivial per-call quality data to render. Synthesising
            # plausible DP / GQ / AD per call keeps the test fixtures
            # within the data shape those tools expect.
            '##FORMAT=<ID=DP,Number=1,Type=Integer,'
            'Description="Approximate read depth">',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,'
            'Description="Genotype quality">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,'
            'Description="Allelic depths for the ref and alt alleles">',
        ])
    lines.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(cohort.samples)
    )
    return lines


def write_vcf(
    output_path: str,
    chrom: str,
    variants: Iterable[Variant],
    cohort: SyntheticCohort | None = None,
    file_date: str = "20130502",
    seed: int = 1000,
    reference: str = DEFAULT_REFERENCE,
    contigs_override: dict[str, int] | None = None,
    info_declarations: list[str] | None = None,
    declare_format_gt: bool = True,
) -> str:
    """Write a bgzipped, tabix-indexed VCF. Returns the .vcf.gz path.

    `output_path` should end in `.vcf.gz`. A matching `.tbi` is generated
    alongside it.

    Knobs relevant to QC-validation tests:
      reference          — string written in `##reference=...`. Set to
                           something non-human (e.g. "mm10") to trigger the
                           unknown-build warning.
      contigs_override   — replace the GRCh37 contig set with a custom dict
                           (e.g. {"chr1_alt": 1000} to trigger the
                           non-human-contig warning).
      info_declarations  — replace the default INFO lines (drop AF/AC/AN to
                           trigger the missing-allele-frequency warning).
      declare_format_gt  — set False to omit the GT FORMAT declaration (to
                           trigger the missing-GT warning).
    """
    if not output_path.endswith(".vcf.gz"):
        raise ValueError("output_path must end in .vcf.gz")
    if not shutil.which("bgzip") or not shutil.which("tabix"):
        raise RuntimeError("bgzip and tabix are required to build test VCFs")

    cohort = cohort or SyntheticCohort()
    rng = random.Random(seed)

    # Variants must be emitted in sorted order for tabix to index correctly.
    vs = sorted(variants, key=lambda v: v.pos)

    plain_path = output_path[:-len(".gz")]  # .vcf
    header_lines = _build_header(
        chrom, cohort, file_date,
        reference=reference,
        contigs_override=contigs_override,
        info_declarations=info_declarations,
        declare_format_gt=declare_format_gt,
    )
    with open(plain_path, "w") as fh:
        for line in header_lines:
            fh.write(line + "\n")
        for v in vs:
            gts = _realize_genotypes(v, cohort, rng)
            ac, an, pop_af = _compute_pop_stats(gts, cohort.samples, cohort)
            info = _format_info(ac, an, pop_af, v)
            per_call = [_format_per_call(g, rng) for g in gts]
            fields = [
                chrom, str(v.pos), v.variant_id, v.ref, v.alt,
                "100", "PASS", info, "GT:DP:GQ:AD", *per_call,
            ]
            fh.write("\t".join(fields) + "\n")

    # bgzip overwrites .vcf.gz if it exists (-f), leaves no .vcf behind.
    subprocess.run(["bgzip", "-f", plain_path], check=True,
                   capture_output=True)
    subprocess.run(["tabix", "-p", "vcf", "-f", output_path], check=True,
                   capture_output=True)
    return output_path


def standard_filename(dir_path: str, chrom: str,
                      date_stamp: str = "20130502") -> str:
    """Return a pipeline-representative filename inside `dir_path`."""
    base = (
        f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b."
        f"{date_stamp}.genotypes.vcf.gz"
    )
    return os.path.join(dir_path, base)


# -----------------------------------------------------------------------------
# CLI — run `python3 synthetic_vcf.py --help` for usage.
# -----------------------------------------------------------------------------

_BASES = ("A", "C", "G", "T")


def _random_variants(n: int, chrom: str, rng: random.Random,
                     pos_start: int | None = None,
                     pos_stop: int | None = None) -> list[Variant]:
    chrom_len = GRCH37_CONTIG_LENGTHS.get(chrom)
    lo = pos_start if pos_start is not None else 1_000_000
    hi = pos_stop if pos_stop is not None else (chrom_len or 50_000_000)
    if hi <= lo:
        raise ValueError(f"pos-stop ({hi}) must be > pos-start ({lo})")
    # Evenly spaced positions avoid accidental duplicates and sort order issues.
    step = max(1, (hi - lo) // (n + 1))
    variants = []
    for i in range(n):
        pos = lo + step * (i + 1)
        ref = rng.choice(_BASES)
        alt = rng.choice([b for b in _BASES if b != ref])
        af = rng.betavariate(0.5, 5.0)  # skewed low, mimicking real MAF spectrum
        variants.append(Variant(pos=pos, ref=ref, alt=alt,
                                af_by_pop={"ALL": af}))
    return variants


def _cli() -> int:
    import argparse
    import sys

    p = argparse.ArgumentParser(
        description=(
            "Generate a synthetic 1000G Phase 3-style VCF (bgzipped + "
            "tabix-indexed)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--chrom", required=True,
                   help="Chromosome (e.g. 15, 22, X). No 'chr' prefix.")
    p.add_argument("--output", "-o",
                   help="Output .vcf.gz path. Default: standard 1000G "
                        "filename in --outdir.")
    p.add_argument("--outdir", default=".",
                   help="Directory for the default filename (ignored if "
                        "--output is set).")
    p.add_argument("--n-variants", type=int, default=100,
                   help="Number of random SNPs to emit.")
    p.add_argument("--pos-start", type=int, default=None,
                   help="Lower bound for random variant positions.")
    p.add_argument("--pos-stop", type=int, default=None,
                   help="Upper bound for random variant positions.")
    p.add_argument("--seed", type=int, default=1000,
                   help="RNG seed — genotypes AND random positions are "
                        "derived from this.")
    p.add_argument("--date-stamp", default="20130502",
                   help="Date stamp used in the default filename.")
    # Optional: splice in a specific known variant at a fixed position/allele.
    p.add_argument("--target-name", default=None,
                   help="Include a named variant at --target-pos with the "
                        "given AFs. Pairs with --target-pos/--target-ref/"
                        "--target-alt/--target-af-*.")
    p.add_argument("--target-pos", type=int, default=None)
    p.add_argument("--target-ref", default=None)
    p.add_argument("--target-alt", default=None)
    p.add_argument("--target-af-eur", type=float, default=0.0)
    p.add_argument("--target-af-afr", type=float, default=0.0)
    p.add_argument("--target-af-eas", type=float, default=0.0)
    p.add_argument("--target-af-sas", type=float, default=0.0)
    p.add_argument("--target-af-amr", type=float, default=0.0)
    args = p.parse_args()

    out_path = args.output or standard_filename(
        args.outdir, args.chrom, args.date_stamp
    )
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    rng = random.Random(args.seed)
    variants = _random_variants(args.n_variants, args.chrom, rng,
                                args.pos_start, args.pos_stop)

    if args.target_name is not None:
        missing = [k for k, v in {
            "--target-pos": args.target_pos,
            "--target-ref": args.target_ref,
            "--target-alt": args.target_alt,
        }.items() if v in (None, "")]
        if missing:
            p.error("--target-name requires: " + ", ".join(missing))
        variants.append(Variant(
            pos=args.target_pos,
            ref=args.target_ref,
            alt=args.target_alt,
            variant_id=args.target_name,
            af_by_pop={
                "EUR": args.target_af_eur,
                "AFR": args.target_af_afr,
                "EAS": args.target_af_eas,
                "SAS": args.target_af_sas,
                "AMR": args.target_af_amr,
            },
        ))

    path = write_vcf(
        out_path, args.chrom, variants,
        cohort=SyntheticCohort(), file_date=args.date_stamp, seed=args.seed,
    )
    n = len(variants)
    print(f"wrote {path} ({n} variant{'s' if n != 1 else ''}, "
          f"{len(DEFAULT_SAMPLES)} samples)", file=sys.stderr)
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
