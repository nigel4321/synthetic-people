"""CLI entry point — wires the package together."""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import random
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
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
from .bcf_writer import CohortBcfWriter
from .truth import TruthBedWriter
from .writer import write_person_vcf


# Shared state that the per-person worker reads. Populated in main()
# before the worker pool is created so the children inherit it via the
# fork start method's copy-on-write semantics. Keep it module-level
# (not closure-captured) so workers under fork pick it up directly,
# without pickling the cohort_sites payload across the task boundary.
_PERSON_WORKER_STATE: dict = {}


def _person_worker(i: int, sid: str, seed: int) -> tuple:
    """Run the per-person work for cohort index ``i``.

    Returns a tuple of ``(person_entry, person_stats, n_svs)``. The
    function reads shared inputs from :data:`_PERSON_WORKER_STATE` so
    they are inherited via fork rather than pickled per task.
    """
    state = _PERSON_WORKER_STATE
    rng = random.Random(seed)

    candidates = state["candidates"]
    cohort_sites = state["cohort_sites"]
    build = state["build"]
    output_dir = state["output_dir"]
    truth_dir = state["truth_dir"]
    contig_order = state["contig_order"]
    svs_per_person = state["svs_per_person"]
    sv_length_max = state["sv_length_max"]
    sv_length_min = state["sv_length_min"]
    sv_chrom_span = state["sv_chrom_span"]
    sv_chromosomes = state["sv_chromosomes"]
    error_rate = state["error_rate"]
    dropout_rate = state["dropout_rate"]
    person_ancestry = state["person_ancestry"]

    hi = dict(rng.choice(candidates))
    hi["gt"] = rng.choices(("0|1", "1|1"), weights=(0.7, 0.3))[0]
    background = person_records_from_cohort(cohort_sites, i)
    person_svs: list = []
    if svs_per_person > 0 and sv_chrom_span > sv_length_max:
        person_svs = generate_person_svs(
            rng, sv_chromosomes, sv_chrom_span,
            n_svs=svs_per_person,
            length_min_bp=sv_length_min,
            length_max_bp=sv_length_max,
        )
        background.extend(person_svs)
    person = {
        "sample_id": sid,
        "highlighted": hi,
        "background": background,
    }
    out = output_dir / f"person_{i+1:04d}.vcf.gz"
    person_stats = new_error_stats()
    golden_path = truth_dir / f"person_{i+1:04d}.golden.bed"
    noise_path = truth_dir / f"person_{i+1:04d}.noise.bed"
    tw = TruthBedWriter(golden_path, noise_path,
                        contig_order=contig_order)
    write_person_vcf(
        out, person, build, rng,
        error_rate=error_rate,
        dropout_rate=dropout_rate,
        stats=person_stats,
        truth_writer=tw,
    )
    tw.close()

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
        bed_path = output_dir / "ancestry" / f"person_{i+1:04d}.bed"
        write_ancestry_bed(bed_path, person_ancestry[i])
        fracs = ancestry_fractions(person_ancestry[i])
        person_entry["ancestry_bed"] = f"ancestry/{bed_path.name}"
        person_entry["ancestry_fractions"] = {
            p: round(v, 4) for p, v in fracs.items()
        }

    return person_entry, person_stats, len(person_svs)


def _format_person_log(entry: dict, n_total: int) -> str:
    """Format the one-line progress log emitted per person."""
    hi = entry["highlighted"]
    base = (
        f"  [{entry['index']:>4}/{n_total}] {entry['vcf']} — "
        f"{entry['sample_id']} — "
        f"highlighted {hi['id']} at {hi['chrom']}:{hi['pos']} "
        f"{hi['ref']}>{hi['alt']} ({hi['gt']}), "
        f"{entry['n_background_records']} background records"
    )
    if "ancestry_fractions" in entry:
        base += ", ancestry " + ",".join(
            f"{p}={v:.2f}" for p, v in entry["ancestry_fractions"].items()
        )
    return base


def submit_overlays(args, chromosomes: list, clinvar_vcf: Path,
                    sig_filter: set,
                    executor: ThreadPoolExecutor) -> dict:
    """Submit the ClinVar / rsID / COSMIC loader futures.

    All three loaders are bcftools subprocess + I/O: they release the
    GIL, so they can run on a thread pool concurrently with the
    compute-heavy coalescent simulation.

    Returns a dict with keys ``clinvar_index``, ``rsid_pool``,
    ``cosmic_pool``. Each value is either a ``Future`` that resolves
    to the loaded structure, or ``None`` if that overlay is skipped
    (rsID with ``--rsid-density 0``; COSMIC unless ``--somatic``).

    The caller must validate ``--somatic`` / ``--cosmic-vcf`` before
    invoking this so a bad path doesn't get silently scheduled. The
    helper trusts its inputs.
    """
    futures: dict = {
        "clinvar_index": None,
        "rsid_pool": None,
        "cosmic_pool": None,
    }
    futures["clinvar_index"] = executor.submit(
        load_clinvar_index, clinvar_vcf, chromosomes,
        sig_filter=sig_filter, max_per_chrom=20_000,
    )
    if args.rsid_density > 0:
        rsid_source = args.dbsnp_vcf or clinvar_vcf
        futures["rsid_pool"] = executor.submit(
            load_rsid_pool, rsid_source, chromosomes,
            max_per_chrom=20_000,
        )
    if args.somatic:
        futures["cosmic_pool"] = executor.submit(
            load_cosmic_records, args.cosmic_vcf, chromosomes,
            max_per_chrom=20_000,
        )
    return futures


def resolve_workers(requested: int) -> int:
    """Return the effective worker count.

    `requested == 0` means auto: use ``os.cpu_count()``. `requested == 1`
    runs serially. On non-Linux hosts we always fall back to 1, because
    Phase-1 parallelism uses ``mp.get_context("fork")`` which is unsafe
    or unsupported elsewhere.
    """
    if requested < 0:
        raise ValueError("--workers must be >= 0")
    import os as _os
    import sys as _sys
    if _sys.platform != "linux":
        return 1
    if requested == 0:
        return max(1, _os.cpu_count() or 1)
    return requested


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
    p.add_argument("--workers", type=int, default=0,
                   help="[perf] Worker processes for per-chromosome "
                        "msprime simulations and per-person VCF writes. "
                        "0 = auto (os.cpu_count()), 1 = serial. Linux "
                        "only — non-Linux hosts fall back to serial. "
                        "Output is deterministic for a given --seed "
                        "regardless of --workers, but differs from a "
                        "pre-Phase-1 run at the same seed because the "
                        "rng consumption pattern changed.")
    p.add_argument("--mode", choices=("per-person", "cohort", "both"),
                   default="per-person",
                   help="Output shape. per-person (default): emit one "
                        "VCF per person, identical to today's behaviour. "
                        "cohort: emit a single multi-sample cohort BCF "
                        "and skip per-person fan-out — derive per-person "
                        "later via `bcftools view -s SAMPLE`. both: emit "
                        "both deliverables. The cohort BCF is the "
                        "scaling-friendly format for large --n; "
                        "per-person stays the default so existing users "
                        "see no behaviour change unless they opt in.")
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

    try:
        workers = resolve_workers(args.workers)
    except ValueError as exc:
        sys.exit(f"--workers: {exc}")
    if args.workers != 0 and args.workers != workers:
        # Requested non-zero workers but resolve_workers downgraded
        # (e.g. non-Linux host). Tell the user.
        print(
            f"  [warn] --workers={args.workers} downgraded to "
            f"{workers} (parallelism is Linux-only in Phase 1)",
            file=sys.stderr,
        )

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

    # Validate --somatic / --cosmic-vcf up front so the overlay
    # prefetch below can trust its inputs and so failures land before
    # the (long) simulation runs.
    if args.somatic and not args.legacy_background:
        if args.cosmic_vcf is None:
            sys.exit("--somatic requires --cosmic-vcf (COSMIC is "
                     "registration-gated; supply a local VCF path)")
        if not args.cosmic_vcf.exists():
            sys.exit(f"--cosmic-vcf not found: {args.cosmic_vcf}")

    # Phase 2: prefetch the overlay loaders concurrently with the
    # simulation. ClinVar / rsID / COSMIC each kick off a bcftools
    # subprocess that releases the GIL, so a thread pool overlaps
    # neatly with the msprime / process-pool work below. Skipped on
    # the legacy path because that path doesn't run overlays.
    overlay_executor: ThreadPoolExecutor | None = None
    overlay_futures: dict = {
        "clinvar_index": None,
        "rsid_pool": None,
        "cosmic_pool": None,
    }
    if not args.legacy_background:
        overlay_executor = ThreadPoolExecutor(
            max_workers=3, thread_name_prefix="overlay")
        overlay_futures = submit_overlays(
            args, chromosomes, clinvar_vcf, sig_filter,
            overlay_executor,
        )
        scheduled = [k for k, v in overlay_futures.items() if v is not None]
        print(
            f"  prefetching overlay loaders in background: "
            f"{', '.join(scheduled)}",
            file=sys.stderr,
        )

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
            rng=rng, verbose=True, workers=workers,
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
            rng=rng, verbose=True, workers=workers,
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
        print(f"Awaiting ClinVar overlay index for {chromosomes}...",
              file=sys.stderr)
        clinvar_index = overlay_futures["clinvar_index"].result()
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
        if overlay_futures["rsid_pool"] is not None:
            rsid_source = args.dbsnp_vcf or clinvar_vcf
            print(f"Awaiting rsID pool from {rsid_source}...",
                  file=sys.stderr)
            rsid_pool = overlay_futures["rsid_pool"].result()
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

        if overlay_futures["cosmic_pool"] is not None:
            print(f"Awaiting COSMIC pool from {args.cosmic_vcf}...",
                  file=sys.stderr)
            cosmic_pool = overlay_futures["cosmic_pool"].result()
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

    if overlay_executor is not None:
        # All futures have been .result()'d above; shutting down here
        # frees the threads while the cohort is still in scope.
        overlay_executor.shutdown(wait=True)

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
    # Pre-derive per-person seeds from the master rng before any
    # per-person work runs. This is what makes the output deterministic
    # for a given --seed regardless of --workers: the per-person rng
    # depends only on its seed, and the seeds are sampled from the
    # master rng in a fixed order.
    person_seeds = [rng.randint(1, 2**31 - 1) for _ in range(args.n)]

    # Phase 5a: write the cohort BCF when --mode is `cohort` or `both`.
    # The BCF carries the truth-state cohort GTs (no per-call DP/GQ/AD
    # noise — those are layered in at per-person derivation time, the
    # same way today's writer.py does). At large --n this is the
    # scaling-friendly deliverable; per-person VCFs can be derived from
    # it later via `bcftools view -s SAMPLE`.
    cohort_bcf_path: Path | None = None
    if args.mode in ("cohort", "both"):
        cohort_dir = args.output_dir / "cohort"
        cohort_bcf_path = cohort_dir / "cohort.bcf"
        contig_index = {c: i for i, c
                        in enumerate(BUILDS[args.build]["contigs"])}
        bcf_t0 = time.monotonic()
        print(
            f"Writing cohort BCF: {args.n} samples × "
            f"{len(cohort_sites)} sites → {cohort_bcf_path}",
            file=sys.stderr,
        )
        # Sort by (contig_order, pos) so the BCF is in genome order
        # regardless of how chromosomes were emitted by the simulator.
        cohort_sites_sorted = sorted(
            cohort_sites,
            key=lambda s: (contig_index.get(s["chrom"],
                                            len(contig_index)), s["pos"]),
        )
        last_log = bcf_t0
        per_chrom_count: dict = {}
        with CohortBcfWriter(cohort_bcf_path, args.build,
                             sample_ids) as bcfw:
            for s in cohort_sites_sorted:
                bcfw.write_site(s)
                per_chrom_count[s["chrom"]] = \
                    per_chrom_count.get(s["chrom"], 0) + 1
                now = time.monotonic()
                # Throttled progress: every ~20 s on long runs, plus a
                # final summary at end. Tells the user "still alive,
                # making progress" without spamming the log on small
                # cohorts that finish in milliseconds.
                if now - last_log > 20.0:
                    written = sum(per_chrom_count.values())
                    rate = written / (now - bcf_t0) if now > bcf_t0 else 0
                    print(
                        f"  cohort BCF: {written:,}/"
                        f"{len(cohort_sites_sorted):,} sites "
                        f"({rate:,.0f}/s)",
                        file=sys.stderr,
                    )
                    last_log = now
        bcf_secs = time.monotonic() - bcf_t0
        for chrom in chromosomes:
            n = per_chrom_count.get(chrom, 0)
            if n:
                print(f"  cohort BCF chrom {chrom}: {n} sites",
                      file=sys.stderr)
        print(f"  cohort BCF complete in {bcf_secs:.1f}s",
              file=sys.stderr)

    if args.mode == "cohort":
        # No per-person fan-out — cohort BCF is the only deliverable.
        # Drop the in-memory cohort_sites payload now and write a slim
        # manifest pointing at the BCF.
        del cohort_sites
        manifest_mode = (
            "legacy-background" if args.legacy_background
            else "admixture-uk" if args.admixture
            else "coalescent"
        )
        manifest = {
            "build": args.build,
            "n_people": args.n,
            "mode": manifest_mode,
            "shape": "cohort",
            "chromosomes": chromosomes,
            "seed": args.seed,
            "samples": sample_ids,
            "cohort_bcfs": [
                str(cohort_bcf_path.relative_to(args.output_dir))
            ],
        }
        if args.admixture:
            manifest["ancestry_proportions"] = {
                "EUR": args.eur_frac,
                "SAS": args.sas_frac,
                "AFR": args.afr_frac,
            }
        if not args.legacy_background:
            manifest["overlays"] = {
                "clinvar_inject_density": args.clinvar_inject_density,
                "rsid_density": args.rsid_density,
                "dbsnp_source": (
                    str(args.dbsnp_vcf) if args.dbsnp_vcf
                    else "clinvar:INFO/RS"),
                "somatic": bool(args.somatic),
                "cosmic_inject_density": (
                    args.cosmic_inject_density if args.somatic else 0.0),
                "stats": overlay_stats,
            }
        manifest_path = args.output_dir / "manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2)
        print(
            f"--mode cohort: skipped per-person VCF writes; cohort BCF "
            f"is the deliverable. Derive per-person VCFs via "
            f"`bcftools view -s SAMPLE_ID {cohort_bcf_path.name}`.",
            file=sys.stderr,
        )
        print(f"Manifest written to {manifest_path}", file=sys.stderr)
        print("Done.", file=sys.stderr)
        return 0

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

    _PERSON_WORKER_STATE.update({
        "candidates": candidates,
        "cohort_sites": cohort_sites,
        "build": args.build,
        "output_dir": args.output_dir,
        "truth_dir": truth_dir,
        "contig_order": contig_order,
        "svs_per_person": args.svs_per_person,
        "sv_length_max": args.sv_length_max,
        "sv_length_min": args.sv_length_min,
        "sv_chrom_span": sv_chrom_span,
        "sv_chromosomes": sv_chromosomes,
        "error_rate": args.error_rate,
        "dropout_rate": args.dropout_rate,
        "person_ancestry": person_ancestry,
    })

    use_pool = workers > 1 and args.n > 1
    fanout_t0 = time.monotonic()
    last_progress_log = fanout_t0
    progress_log_interval = 20.0   # seconds — same cadence as cohort BCF

    def _maybe_log_progress(done: int) -> None:
        nonlocal last_progress_log
        now = time.monotonic()
        # First / last person always log; intermediate progress is
        # throttled so the output stays readable on small cohorts.
        if (now - last_progress_log) < progress_log_interval \
                and 0 < done < args.n:
            return
        elapsed = now - fanout_t0
        rate = done / elapsed if elapsed > 0 else 0.0
        remaining = args.n - done
        eta = remaining / rate if rate > 0 else float("inf")
        eta_str = f"{eta:.0f}s" if eta < float("inf") else "?"
        print(
            f"  person VCFs: {done:,}/{args.n:,} written "
            f"({rate:.1f}/s, elapsed {elapsed:.0f}s, eta {eta_str})",
            file=sys.stderr,
        )
        last_progress_log = now

    if use_pool:
        print(f"  fan-out: {workers} worker processes "
              f"(fork start method, --n={args.n})",
              file=sys.stderr)
        ctx = mp.get_context("fork")
        results: list = [None] * args.n
        with ProcessPoolExecutor(max_workers=workers,
                                 mp_context=ctx) as ex:
            futures = [
                ex.submit(_person_worker, i, sample_ids[i],
                          person_seeds[i])
                for i in range(args.n)
            ]
            # Iterate in submission order so progress logs come out in
            # person order even if workers complete out of order.
            for i, fut in enumerate(futures):
                results[i] = fut.result()
                _maybe_log_progress(i + 1)
    else:
        results = []
        for i in range(args.n):
            results.append(_person_worker(i, sample_ids[i], person_seeds[i]))
            _maybe_log_progress(i + 1)

    # Drop the shared payload now that workers have finished, so the
    # cohort_sites reference can be GC'd before manifest writing.
    _PERSON_WORKER_STATE.clear()

    manifest_people: list = []
    for entry, person_stats, n_svs in results:
        merge_stats(error_stats_total, person_stats)
        sv_total += n_svs
        manifest_people.append(entry)
        print(_format_person_log(entry, args.n), file=sys.stderr)

    mode = ("legacy-background" if args.legacy_background
            else "admixture-uk" if args.admixture
            else "coalescent")
    manifest = {
        "build": args.build,
        "n_people": args.n,
        "mode": mode,
        "shape": args.mode,
        "chromosomes": chromosomes,
        "seed": args.seed,
        # Always emit the flat list of sample IDs at the top level so a
        # downstream tool wanting "give me every sample" doesn't need a
        # different code path per --mode value (cohort returns
        # `samples`; per-person/both used to require iterating
        # `people[*].sample_id`).
        "samples": sample_ids,
        "people": manifest_people,
    }
    if cohort_bcf_path is not None:
        manifest["cohort_bcfs"] = [
            str(cohort_bcf_path.relative_to(args.output_dir))
        ]
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
