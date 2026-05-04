"""CLI entry point — wires the package together."""

from __future__ import annotations

import argparse
import json
import random
import shutil
import sys
from pathlib import Path

from .admixture import (
    DEFAULT_AFR_FRAC,
    DEFAULT_EUR_FRAC,
    DEFAULT_SAS_FRAC,
    ancestry_fractions,
    simulate_cohort as simulate_admixed_cohort,
    write_ancestry_bed,
)
from .background import load_background_pool, random_sample_id
from .builds import BUILDS
from .clinvar import (
    DEFAULT_CLINVAR_INJECT_DENSITY,
    DEFAULT_SIG_FILTER,
    annotate_clinvar,
    fetch_clinvar,
    inject_clinvar,
    load_clinvar_index,
    load_highlighted_candidates,
)
from .cosmic import (
    DEFAULT_COSMIC_INJECT_DENSITY,
    inject_cosmic,
    load_cosmic_records,
)
from .dbsnp import (
    DEFAULT_RSID_DENSITY,
    inject_rsids,
    load_rsid_pool,
)
from .errors import (
    DEFAULT_DROPOUT_RATE,
    DEFAULT_GT_ERROR_RATE,
    merge_stats,
    new_error_stats,
)
from .coalescent import (
    DEFAULT_CHR_LENGTH_MB,
    DEFAULT_DEMO_MODEL,
    DEFAULT_MU,
    DEFAULT_POPULATION,
    DEFAULT_REC_RATE,
    simulate_cohort,
)
from .cohort import draw_cohort_background, person_records_from_cohort
from .sfs import (
    DEFAULT_SFS_ALPHA,
    sfs_histogram,
    singleton_fraction,
    write_sfs_tsv,
)
from .sv import (
    DEFAULT_SV_LENGTH_MAX_BP,
    DEFAULT_SV_LENGTH_MIN_BP,
    DEFAULT_SVS_PER_PERSON,
    generate_person_svs,
)
from .truth import TruthBedWriter
from .writer import write_person_vcf


def parse_chromosomes(spec: str, build: str) -> list[str]:
    """Parse a --chromosomes spec into an ordered, deduplicated list.

    Accepts comma-separated tokens. Each token is either a single contig
    (e.g. ``22``, ``X``) or a numeric range like ``1-10`` (inclusive).
    Ranges only span autosomes; mix singletons with ranges freely:
    ``1-3,5,19-22,X``.

    Every resulting contig must exist in the chosen build's contig table.
    """
    contigs = BUILDS[build]["contigs"]
    seen: set = set()
    out: list = []
    for raw in spec.split(","):
        token = raw.strip()
        if not token:
            continue
        if "-" in token and not token.startswith("-"):
            lo_str, hi_str = token.split("-", 1)
            lo_str, hi_str = lo_str.strip(), hi_str.strip()
            if not (lo_str.isdigit() and hi_str.isdigit()):
                raise ValueError(
                    f"chromosome range {token!r} must be numeric "
                    "(e.g. '1-22')"
                )
            lo, hi = int(lo_str), int(hi_str)
            if lo > hi:
                raise ValueError(
                    f"chromosome range {token!r} is empty: "
                    f"start {lo} > end {hi}"
                )
            members = [str(i) for i in range(lo, hi + 1)]
        else:
            members = [token]
        for c in members:
            if c not in contigs:
                raise ValueError(
                    f"unknown chromosome {c!r} for build {build}; "
                    f"valid: {sorted(contigs)}"
                )
            if c not in seen:
                seen.add(c)
                out.append(c)
    if not out:
        raise ValueError("--chromosomes resolved to an empty list")
    return out


OPTIONAL_PYTHON_DEPS = [
    # These aren't used in M1 but will be required from M2+. `--check-deps`
    # reports their absence so a user can install ahead of time.
    ("numpy", "M2+ (quality metrics, sampling)"),
    ("msprime", "M5+ (coalescent simulation)"),
    ("tskit", "M5+ (tree-sequence handling)"),
    ("stdpopsim", "M6+ (human demographic models)"),
    ("matplotlib", "M10 (validation plots)"),
    ("allel", "M10 (scikit-allel PCA / LD decay)"),
]


def _check_deps(verbose: bool = True) -> int:
    """Check required binaries and (optionally) Python deps. Returns exit code."""
    missing_bins: list[str] = []
    for tool in ("bcftools", "tabix", "bgzip"):
        if not shutil.which(tool):
            missing_bins.append(tool)

    missing_py: list[tuple[str, str]] = []
    for mod, reason in OPTIONAL_PYTHON_DEPS:
        try:
            __import__(mod)
        except ImportError:
            missing_py.append((mod, reason))

    if verbose:
        if not missing_bins:
            print("htslib binaries: OK (bcftools, tabix, bgzip)",
                  file=sys.stderr)
        else:
            print("htslib binaries MISSING: " + ", ".join(missing_bins),
                  file=sys.stderr)
        if not missing_py:
            print("optional Python deps: all present", file=sys.stderr)
        else:
            print("optional Python deps (not required for M1):",
                  file=sys.stderr)
            for mod, reason in missing_py:
                print(f"  - {mod:<12} needed for {reason}", file=sys.stderr)
    # Only hard-fail on missing binaries. Python deps are optional in M1.
    return 1 if missing_bins else 0


def _parser(script_dir: Path) -> argparse.ArgumentParser:
    default_bg = str(script_dir.parent / "ALL.chr*.phase3_*.genotypes.vcf.gz")

    p = argparse.ArgumentParser(
        description=(
            "Generate a cohort of synthetic person VCFs. The default path "
            "simulates an LD-correct coalescent with stdpopsim demography; "
            "pass --legacy-background for the M4 1000G-pool + power-law "
            "SFS sampler. Each person receives a ClinVar-highlighted "
            "variant on top of the shared cohort background."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--n", type=int, default=10,
                   help="Cohort size: number of person VCFs to generate")
    p.add_argument("--output-dir", type=Path,
                   default=script_dir / "out",
                   help="Where to write person_<N>.vcf.gz")
    p.add_argument("--cache-dir", type=Path,
                   default=script_dir / "cache",
                   help="Where ClinVar is downloaded and cached")
    p.add_argument("--build", choices=list(BUILDS), default="GRCh38",
                   help="Reference build; must match background VCFs")
    p.add_argument("--seed", type=int, default=None,
                   help="RNG seed for deterministic output. Omit for "
                        "different people each run.")
    p.add_argument("--background-glob", action="append", default=None,
                   help="[legacy] Glob(s) for common-variant source VCFs. "
                        "Pass multiple times to combine sources. "
                        f"Default: {default_bg}")
    p.add_argument("--n-background", type=int, default=500,
                   help="[legacy] Number of shared background sites the "
                        "cohort is drawn over (hom-ref calls dropped "
                        "per person)")
    p.add_argument("--af-min", type=float, default=0.05,
                   help="[legacy] Minimum AF filter when loading the "
                        "coordinate pool from 1000G sources")
    p.add_argument("--sfs-alpha", type=float, default=DEFAULT_SFS_ALPHA,
                   help="[legacy path] Power-law exponent for the cohort "
                        "SFS: P(k) ∝ 1/k^α. α=1.0 is Watterson-neutral; "
                        "α=2.0 (default) biases toward singletons, "
                        "matching gnomAD-like spectra.")
    p.add_argument("--legacy-background", action="store_true",
                   help="Use the M4 1000G-pool + power-law SFS sampler "
                        "instead of the coalescent. No LD, no realistic "
                        "demography — retained for comparison and "
                        "offline-only use.")
    p.add_argument("--chromosomes", default="22",
                   help="[coalescent] Chromosomes to simulate. Accepts a "
                        "comma-separated list, numeric ranges, or both: "
                        "'22', '19,20,21,22', '1-22', '1-3,5,19-22,X'. "
                        "Shorter list = faster run.")
    p.add_argument("--chr-length-mb", type=float,
                   default=DEFAULT_CHR_LENGTH_MB,
                   help="[coalescent] Simulated prefix length per "
                        "chromosome in Mb. 0 = full length (slow on big "
                        "chromosomes, fine for chr22).")
    p.add_argument("--demo-model", default=DEFAULT_DEMO_MODEL,
                   help="[coalescent] stdpopsim demographic model id. "
                        "Pass 'none' for a constant-size Ne=10k "
                        "single-pop msprime draw (no real demography).")
    p.add_argument("--population", default=DEFAULT_POPULATION,
                   help="[coalescent] Sampling population within the demo "
                        "model (e.g. CEU / YRI / CHB for OutOfAfrica_3G09).")
    p.add_argument("--rec-rate", type=float, default=DEFAULT_REC_RATE,
                   help="[coalescent, --demo-model=none only] Uniform "
                        "recombination rate per bp per generation.")
    p.add_argument("--mu", type=float, default=DEFAULT_MU,
                   help="[coalescent, --demo-model=none only] Mutation "
                        "rate per bp per generation.")
    p.add_argument("--admixture", action="store_true",
                   help="UK-cohort admixture mode: simulate EUR + SAS + "
                        "AFR sources mixed via a single pulse and write "
                        "per-person local-ancestry BED truth tracks. "
                        "Overrides --demo-model / --population.")
    p.add_argument("--eur-frac", type=float, default=DEFAULT_EUR_FRAC,
                   help="[admixture] EUR ancestry proportion in the "
                        "admixture pulse.")
    p.add_argument("--sas-frac", type=float, default=DEFAULT_SAS_FRAC,
                   help="[admixture] SAS ancestry proportion in the "
                        "admixture pulse.")
    p.add_argument("--afr-frac", type=float, default=DEFAULT_AFR_FRAC,
                   help="[admixture] AFR ancestry proportion in the "
                        "admixture pulse. EUR+SAS+AFR must sum to 1.0.")
    p.add_argument("--clinvar-sig",
                   default=",".join(sorted(DEFAULT_SIG_FILTER)),
                   help="Comma-separated CLNSIG values to include when "
                        "drawing highlighted variants")
    p.add_argument("--clinvar-inject-density", type=float,
                   default=DEFAULT_CLINVAR_INJECT_DENSITY,
                   help="[M7] Fraction of cohort sites to overwrite "
                        "with random ClinVar pathogenic records, so "
                        "CLNSIG/CLNDN appear at realistic chromosome "
                        "coordinates. Set to 0 to skip injection (the "
                        "highlighted per-person variant still lands).")
    p.add_argument("--rsid-density", type=float,
                   default=DEFAULT_RSID_DENSITY,
                   help="[M7] Fraction of cohort sites to overwrite "
                        "with a dbSNP-known variant (real coordinates "
                        "and rsID). Set to 0 to skip rsID injection — "
                        "all background record IDs will then be '.'.")
    p.add_argument("--dbsnp-vcf", type=Path, default=None,
                   help="[M7] Path to a dbSNP-style VCF (rsIDs in the "
                        "ID column) for rsID injection. If omitted, "
                        "the cached ClinVar VCF is used (its INFO/RS "
                        "field carries dbSNP rs numbers).")
    p.add_argument("--somatic", action="store_true",
                   help="[M7] Enable COSMIC overlay/injection. Requires "
                        "--cosmic-vcf because COSMIC is registration-"
                        "gated and we never auto-fetch.")
    p.add_argument("--cosmic-vcf", type=Path, default=None,
                   help="[M7] Path to a COSMIC-format VCF. Required "
                        "when --somatic is set.")
    p.add_argument("--cosmic-inject-density", type=float,
                   default=DEFAULT_COSMIC_INJECT_DENSITY,
                   help="[M7] Fraction of cohort sites to overwrite "
                        "with COSMIC records when --somatic is set.")
    p.add_argument("--svs-per-person", type=int,
                   default=DEFAULT_SVS_PER_PERSON,
                   help="[M8] Number of structural variants (DEL / "
                        "DUP / INV) to emit per person, drawn at "
                        "random positions inside the simulated "
                        "region. Set to 0 to skip SV emission.")
    p.add_argument("--sv-length-min", type=int,
                   default=DEFAULT_SV_LENGTH_MIN_BP,
                   help="[M8] Minimum SV length in bp (log-uniform).")
    p.add_argument("--sv-length-max", type=int,
                   default=DEFAULT_SV_LENGTH_MAX_BP,
                   help="[M8] Maximum SV length in bp (log-uniform).")
    p.add_argument("--error-rate", type=float,
                   default=DEFAULT_GT_ERROR_RATE,
                   help="[M9] Per-call probability of a genotype flip "
                        "(false positive / negative). Applied after AD "
                        "is drawn from the truth, so flipped calls land "
                        "low-GQ. 0 disables.")
    p.add_argument("--dropout-rate", type=float,
                   default=DEFAULT_DROPOUT_RATE,
                   help="[M9] Per-call probability of a coverage "
                        "dropout (DP=0, GT=./., GQ=0). 0 disables.")
    p.add_argument("--art", action="store_true",
                   help="[M9, heavy] Use ART read simulation + "
                        "bcftools call instead of the lightweight "
                        "noise model. Requires a reference FASTA on "
                        "disk and the `art_illumina` binary on PATH; "
                        "wired in M11 alongside the GRCh38 reference.")
    p.add_argument("--check-deps", action="store_true",
                   help="Check htslib binaries and optional Python deps, "
                        "then exit")
    p.set_defaults(_default_bg=default_bg)
    return p


def main(argv: list[str] | None = None) -> int:
    script_dir = Path(__file__).resolve().parent.parent
    args = _parser(script_dir).parse_args(argv)

    if args.check_deps:
        return _check_deps()

    # Hard-fail on missing htslib binaries even without --check-deps.
    for tool in ("bcftools", "tabix", "bgzip"):
        if not shutil.which(tool):
            sys.exit(f"required tool not on PATH: {tool}")

    if args.art:
        sys.exit(
            "--art (ART read simulation + bcftools call) requires the "
            "GRCh38 reference FASTA, which is wired in M11. Use the "
            "default lightweight noise model (--error-rate / "
            "--dropout-rate) for now."
        )

    try:
        chromosomes = parse_chromosomes(args.chromosomes, args.build)
    except ValueError as exc:
        sys.exit(f"--chromosomes: {exc}")

    rng = random.Random(args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reference build: {args.build}", file=sys.stderr)
    print("Fetching ClinVar (cached across runs)...", file=sys.stderr)
    clinvar_vcf = fetch_clinvar(args.cache_dir, args.build)

    sig_filter = {s.strip() for s in args.clinvar_sig.split(",") if s.strip()}
    print(
        f"Loading highlighted candidates (CLNSIG in {sorted(sig_filter)})...",
        file=sys.stderr,
    )
    candidates = load_highlighted_candidates(clinvar_vcf, sig_filter)
    print(f"  {len(candidates)} candidates matched", file=sys.stderr)
    if not candidates:
        sys.exit("no ClinVar variants matched the CLNSIG filter — widen "
                 "--clinvar-sig")

    person_ancestry: list = []  # admixture path fills per-person segments
    if args.legacy_background:
        bg_globs = args.background_glob or [args._default_bg]
        print(f"[legacy] Sampling background coordinates from: {bg_globs}",
              file=sys.stderr)
        background_pool = load_background_pool(
            bg_globs, args.af_min, per_source_limit=5000, rng=rng,
        )
        print(f"  coordinate pool: {len(background_pool)} variants",
              file=sys.stderr)
        if not background_pool:
            print("  [warn] no background sources matched — output VCFs "
                  "will contain only the highlighted variant",
                  file=sys.stderr)
        print(
            f"[legacy] Drawing cohort background: {args.n_background} "
            f"shared sites across {args.n} people (α={args.sfs_alpha})",
            file=sys.stderr,
        )
        cohort_sites = draw_cohort_background(
            background_pool, args.n, args.n_background, args.sfs_alpha,
            rng,
        )
    elif args.admixture:
        proportions = (args.eur_frac, args.sas_frac, args.afr_frac)
        print(
            f"Simulating UK admixed cohort: {args.n} people across "
            f"chroms {chromosomes} "
            f"(EUR={args.eur_frac:.2f}, SAS={args.sas_frac:.2f}, "
            f"AFR={args.afr_frac:.2f}, length={args.chr_length_mb} Mb)",
            file=sys.stderr,
        )
        cohort_sites, person_ancestry = simulate_admixed_cohort(
            chromosomes=chromosomes, build=args.build,
            n_people=args.n, length_mb=args.chr_length_mb,
            proportions=proportions,
            rec_rate=args.rec_rate, mu=args.mu,
            rng=rng, verbose=True,
        )
    else:
        demo_model = None if args.demo_model.lower() == "none" \
            else args.demo_model
        print(
            f"Simulating coalescent cohort: {args.n} people across "
            f"chroms {chromosomes} "
            f"(model={demo_model or 'uniform-constant-Ne'}, "
            f"pop={args.population}, length={args.chr_length_mb} Mb)",
            file=sys.stderr,
        )
        cohort_sites = simulate_cohort(
            chromosomes=chromosomes, build=args.build,
            n_people=args.n, length_mb=args.chr_length_mb,
            demo_model=demo_model, population=args.population,
            rec_rate=args.rec_rate, mu=args.mu,
            rng=rng, verbose=True,
        )

    overlay_stats = {
        "clinvar_annotated": 0,
        "clinvar_injected": 0,
        "rsid_injected": 0,
        "cosmic_injected": 0,
    }
    if not args.legacy_background:
        # ClinVar annotation: catch any natural collision before injection
        # rewrites coordinates. With msprime simulating positions
        # 1..sim_length and ClinVar at real chromosome coordinates the
        # natural collision rate is essentially zero, but we still run
        # annotate_clinvar so a future demography that lands on real
        # coordinates picks up CLNSIG for free.
        print("Loading ClinVar overlay index for "
              f"{chromosomes}...", file=sys.stderr)
        clinvar_index = load_clinvar_index(
            clinvar_vcf, chromosomes, sig_filter=sig_filter,
            max_per_chrom=20_000,
        )
        print(f"  {len(clinvar_index)} ClinVar pathogenic records "
              "available for overlay", file=sys.stderr)
        overlay_stats["clinvar_annotated"] = annotate_clinvar(
            cohort_sites, clinvar_index)
        if overlay_stats["clinvar_annotated"]:
            print(f"  annotated {overlay_stats['clinvar_annotated']} "
                  "natural ClinVar collisions", file=sys.stderr)

        if args.clinvar_inject_density > 0 and clinvar_index:
            overlay_stats["clinvar_injected"] = inject_clinvar(
                cohort_sites, clinvar_index,
                args.clinvar_inject_density, rng,
            )
            print(f"  injected {overlay_stats['clinvar_injected']} "
                  "ClinVar pathogenic records "
                  f"(density={args.clinvar_inject_density:.3f})",
                  file=sys.stderr)

        # Reserve ClinVar-injected rows so subsequent overlays don't
        # overwrite their (pos, ref, alt) and break the CLNSIG↔variant
        # correspondence. inject_clinvar re-sorts sites on exit, so we
        # rebuild the index set against the post-sort list.
        clinvar_reserved = {i for i, s in enumerate(cohort_sites)
                            if s.get("clnsig")}

        # rsID injection: prefer a user-supplied dbSNP VCF, fall back to
        # ClinVar's own INFO/RS field (no extra download required).
        if args.rsid_density > 0:
            rsid_source = args.dbsnp_vcf or clinvar_vcf
            print(f"Loading rsID pool from {rsid_source}...",
                  file=sys.stderr)
            rsid_pool = load_rsid_pool(rsid_source, chromosomes,
                                       max_per_chrom=20_000)
            print(f"  {len(rsid_pool)} rsID-bearing records available",
                  file=sys.stderr)
            overlay_stats["rsid_injected"] = inject_rsids(
                cohort_sites, rsid_pool, args.rsid_density, rng,
                reserve_indices=clinvar_reserved,
            )
            print(f"  injected {overlay_stats['rsid_injected']} "
                  f"rsIDs (density={args.rsid_density:.3f})",
                  file=sys.stderr)

        # Re-derive reservations to include rsID-injected rows for the
        # COSMIC pass so all three overlays land on disjoint sites.
        all_overlay_reserved = {i for i, s in enumerate(cohort_sites)
                                if s.get("clnsig") or
                                s["id"].startswith("rs")}

        if args.somatic:
            if args.cosmic_vcf is None:
                sys.exit("--somatic requires --cosmic-vcf (COSMIC is "
                         "registration-gated; supply a local VCF path)")
            if not args.cosmic_vcf.exists():
                sys.exit(f"--cosmic-vcf not found: {args.cosmic_vcf}")
            print(f"Loading COSMIC pool from {args.cosmic_vcf}...",
                  file=sys.stderr)
            cosmic_pool = load_cosmic_records(
                args.cosmic_vcf, chromosomes, max_per_chrom=20_000)
            print(f"  {len(cosmic_pool)} COSMIC records available",
                  file=sys.stderr)
            overlay_stats["cosmic_injected"] = inject_cosmic(
                cohort_sites, cosmic_pool,
                args.cosmic_inject_density, rng,
                reserve_indices=all_overlay_reserved,
            )
            print(f"  injected {overlay_stats['cosmic_injected']} "
                  "COSMIC records "
                  f"(density={args.cosmic_inject_density:.3f})",
                  file=sys.stderr)

    hist = sfs_histogram(cohort_sites)
    total_alts = sum(hist.values())
    singletons = hist.get(1, 0)
    sfrac = singleton_fraction(hist)
    print(
        f"  cohort sites: {len(cohort_sites)}; alt observations: "
        f"{total_alts}; singletons: {singletons} ({sfrac:.1%})",
        file=sys.stderr,
    )

    summary_dir = args.output_dir / "summary"
    sfs_path = summary_dir / "sfs.tsv"
    write_sfs_tsv(sfs_path, hist)
    print(f"  SFS histogram written to {sfs_path}", file=sys.stderr)

    sample_ids = [random_sample_id(rng) for _ in range(args.n)]

    # SV bounds: each SV occupies [pos, pos + svlen]; we draw POS up to
    # `chrom_length_bp - sv_length_max` to keep the END inside the
    # simulated region. For the legacy path we just use the full
    # contig length (no sim window).
    if args.legacy_background:
        sv_chrom_span = max(
            (BUILDS[args.build]["contigs"][c] for c in chromosomes),
            default=0,
        )
    else:
        sv_chrom_span = int(args.chr_length_mb * 1_000_000) \
            if args.chr_length_mb > 0 else max(
                BUILDS[args.build]["contigs"].values())
    sv_chromosomes = chromosomes
    sv_total = 0

    print(f"Writing {args.n} person VCFs into {args.output_dir}",
          file=sys.stderr)
    if args.error_rate > 0 or args.dropout_rate > 0:
        print(
            f"Sequencing-error model: error_rate={args.error_rate:.4f}, "
            f"dropout_rate={args.dropout_rate:.4f}",
            file=sys.stderr,
        )
    error_stats_total = new_error_stats()
    truth_dir = args.output_dir / "truth"
    contig_order = {c: i for i, c
                    in enumerate(BUILDS[args.build]["contigs"])}
    manifest_people: list = []
    for i, sid in enumerate(sample_ids):
        hi = dict(rng.choice(candidates))
        hi["gt"] = rng.choices(("0|1", "1|1"), weights=(0.7, 0.3))[0]
        background = person_records_from_cohort(cohort_sites, i)
        person_svs: list = []
        if args.svs_per_person > 0 and sv_chrom_span > args.sv_length_max:
            person_svs = generate_person_svs(
                rng, sv_chromosomes, sv_chrom_span,
                n_svs=args.svs_per_person,
                length_min_bp=args.sv_length_min,
                length_max_bp=args.sv_length_max,
            )
            background.extend(person_svs)
            sv_total += len(person_svs)
        person = {
            "sample_id": sid,
            "highlighted": hi,
            "background": background,
        }
        out = args.output_dir / f"person_{i+1:04d}.vcf.gz"
        person_stats = new_error_stats()
        golden_path = truth_dir / f"person_{i+1:04d}.golden.bed"
        noise_path = truth_dir / f"person_{i+1:04d}.noise.bed"
        tw = TruthBedWriter(golden_path, noise_path,
                            contig_order=contig_order)
        write_person_vcf(
            out, person, args.build, rng,
            error_rate=args.error_rate,
            dropout_rate=args.dropout_rate,
            stats=person_stats,
            truth_writer=tw,
        )
        tw.close()
        merge_stats(error_stats_total, person_stats)

        person_entry = {
            "index": i + 1,
            "sample_id": sid,
            "vcf": out.name,
            "highlighted": {
                "id": hi["id"],
                "chrom": hi["chrom"],
                "pos": hi["pos"],
                "ref": hi["ref"],
                "alt": ",".join(hi["alts"]),
                "gt": hi["gt"],
            },
            "n_background_records": len(background),
            "n_svs": len(person_svs),
            "errors": dict(person_stats),
            "golden_bed": f"truth/{golden_path.name}",
            "noise_bed": f"truth/{noise_path.name}",
            "n_golden": tw.golden_count,
            "n_noise": tw.noise_count,
        }

        if person_ancestry:
            bed_path = args.output_dir / "ancestry" / \
                f"person_{i+1:04d}.bed"
            write_ancestry_bed(bed_path, person_ancestry[i])
            fracs = ancestry_fractions(person_ancestry[i])
            person_entry["ancestry_bed"] = \
                f"ancestry/{bed_path.name}"
            person_entry["ancestry_fractions"] = {
                p: round(v, 4) for p, v in fracs.items()
            }
            print(
                f"  [{i+1:>4}/{args.n}] {out.name} — {sid} — "
                f"highlighted {hi['id']} at {hi['chrom']}:{hi['pos']} "
                f"{hi['ref']}>{','.join(hi['alts'])} ({hi['gt']}), "
                f"{len(background)} background records, "
                f"ancestry "
                + ",".join(f"{p}={v:.2f}"
                           for p, v in person_entry["ancestry_fractions"]
                           .items()),
                file=sys.stderr,
            )
        else:
            print(
                f"  [{i+1:>4}/{args.n}] {out.name} — {sid} — "
                f"highlighted {hi['id']} at {hi['chrom']}:{hi['pos']} "
                f"{hi['ref']}>{','.join(hi['alts'])} ({hi['gt']}), "
                f"{len(background)} background records",
                file=sys.stderr,
            )

        manifest_people.append(person_entry)

    mode = ("legacy-background" if args.legacy_background
            else "admixture-uk" if args.admixture
            else "coalescent")
    manifest = {
        "build": args.build,
        "n_people": args.n,
        "mode": mode,
        "chromosomes": chromosomes,
        "seed": args.seed,
        "people": manifest_people,
    }
    if args.admixture:
        manifest["ancestry_proportions"] = {
            "EUR": args.eur_frac,
            "SAS": args.sas_frac,
            "AFR": args.afr_frac,
        }
    manifest["svs"] = {
        "per_person": args.svs_per_person,
        "length_min_bp": args.sv_length_min,
        "length_max_bp": args.sv_length_max,
        "total": sv_total,
    }
    realised_fdr = (
        (error_stats_total["flipped"] + error_stats_total["dropped"]) /
        error_stats_total["total_calls"]
    ) if error_stats_total["total_calls"] else 0.0
    manifest["errors"] = {
        "mode": "art" if args.art else "lightweight",
        "error_rate": args.error_rate,
        "dropout_rate": args.dropout_rate,
        "stats": dict(error_stats_total),
        "realised_fdr": round(realised_fdr, 6),
    }
    if args.error_rate > 0 or args.dropout_rate > 0:
        print(
            f"Sequencing-error stats: "
            f"flipped={error_stats_total['flipped']}, "
            f"dropped={error_stats_total['dropped']}, "
            f"total_calls={error_stats_total['total_calls']}, "
            f"realised_fdr={realised_fdr:.4%}",
            file=sys.stderr,
        )
    if not args.legacy_background:
        manifest["overlays"] = {
            "clinvar_inject_density": args.clinvar_inject_density,
            "rsid_density": args.rsid_density,
            "dbsnp_source": (str(args.dbsnp_vcf) if args.dbsnp_vcf
                             else "clinvar:INFO/RS"),
            "somatic": bool(args.somatic),
            "cosmic_inject_density": (args.cosmic_inject_density
                                      if args.somatic else 0.0),
            "stats": overlay_stats,
        }
    manifest_path = args.output_dir / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"Manifest written to {manifest_path}", file=sys.stderr)

    print("Done.", file=sys.stderr)
    return 0
