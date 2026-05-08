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
from .background import draw_sample_ids, load_background_pool
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
    auto_derate_workers,
    auto_pick_chunk_size_mb,
    simulate_cohort,
    simulate_cohort_iter,
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
from .bcf_writer import CohortBcfWriter, write_cohort_bcf_parallel
from .cohort_derivation import derive_person_records, derive_persons_batch
from .memprofile import mark as memprofile_mark
from .resume import ResumeMismatch, load_or_create_meta
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
    # Phase 5b2: workers source their per-person background either from
    # an in-memory cohort_sites list (legacy / admixture paths) or by
    # querying the streamed per-chrom cohort BCFs from disk
    # (non-legacy non-admixture coalescent path). The flow is otherwise
    # identical — same per-person dict shape, same downstream writer
    # call, same truth-BED writer.
    cohort_sites = state.get("cohort_sites")
    cohort_bcfs = state.get("cohort_bcfs")
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
    # Phase 5g batched-extraction path: when the parent has pre-staged
    # this batch's per-person record dicts in ``batch_backgrounds``,
    # workers pick up their share via fork-inherited state and skip
    # the per-person bcftools fan-out entirely. Falls back to the
    # per-person ``derive_person_records`` path for callers that
    # haven't been migrated (legacy admixture etc.).
    batch_backgrounds = state.get("batch_backgrounds")
    if batch_backgrounds is not None and sid in batch_backgrounds:
        background = batch_backgrounds[sid]
    elif cohort_bcfs is not None:
        background = derive_person_records(cohort_bcfs, sid)
    else:
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


def _format_duration(seconds: float) -> str:
    """Render a duration in ``Hh Mm Ss`` / ``Mm Ss`` / ``Ss`` form.

    Used by progress logs so multi-hour ETAs read at a glance — a
    raw ``eta 27970s`` is harder to skim than ``eta 7h 46m 10s``.
    Returns ``"?"`` for ``inf`` (the rate-is-zero case) so callers
    can drop their own infinity-guard branches.
    """
    if seconds == float("inf") or seconds != seconds:  # inf or NaN
        return "?"
    total = int(round(seconds))
    if total < 0:
        total = 0
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours}h {minutes}m {secs}s"
    if minutes:
        return f"{minutes}m {secs}s"
    return f"{secs}s"


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
    p.add_argument("--chr-chunk-mb", type=float, default=0.0,
                   help="[coalescent, perf] Split each chromosome's "
                        "msprime simulation into independent sub-chunks "
                        "of this size (Mb). Bounds per-chunk peak RAM "
                        "for large cohorts. 0 (default) = auto-pick "
                        "from `psutil.virtual_memory().available` and "
                        "`--workers`; positive value overrides. Larger "
                        "chunks preserve cross-chunk LD better but use "
                        "more peak RAM during simulation. See "
                        "PERFORMANCE_PLAN §\"Phase 5f\" for the "
                        "boundary-smoothing semantics.")
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
                   help="[perf] Worker processes for parallel cohort "
                        "BCF writes (sample-slice, post-Phase-5e) and "
                        "per-person VCF writes. 0 = auto "
                        "(os.cpu_count()), 1 = serial. Linux only — "
                        "non-Linux hosts fall back to serial. "
                        "Note: msprime simulation itself is single-"
                        "threaded and runs serially across "
                        "chromosomes regardless of --workers — "
                        "increasing --workers parallelises only the "
                        "cohort BCF write and the per-person fan-out, "
                        "not the simulation. Output is deterministic "
                        "for a given --seed across any --workers "
                        "value.")
    p.add_argument("--fanout-batch-size", type=int, default=4,
                   help="[perf] On the streamed coalescent path, group "
                        "this many sample IDs into one bcftools query "
                        "per cohort BCF during the per-person fan-out. "
                        "The legacy per-person path issued one bcftools "
                        "subprocess per (person, chrom) pair (~66k for "
                        "n=3000 × 22 chroms); batching cuts that by ~B× "
                        "while bounding the parent process's batch RSS "
                        "at B × per-person record list (empirically "
                        "~280 MB/person at n=3000 × 22 chroms). The "
                        "binding memory ceiling is "
                        "(parent_baseline + B × per_person) × workers, "
                        "because each fork-spawned worker's apparent "
                        "RSS includes the parent's batch via COW and "
                        "the kernel's OOM scoring is per-process. "
                        "Default 4 keeps the worst case under ~14 GB "
                        "with --workers 8 on a 32 GB host. Raise it "
                        "for hosts with more RAM or for runs with "
                        "fewer workers; lower it (or drop --workers) "
                        "if the fan-out OOMs.")
    p.add_argument("--profile-memory", type=Path, default=None,
                   metavar="TSV_PATH",
                   help="[diagnostic] Spawn a background thread that "
                        "samples this process's RSS once per second to "
                        "TSV_PATH (one row per sample, plus labelled "
                        "checkpoint rows at key code transitions). "
                        "Flushed + fsynced after every write so an OOM "
                        "kill preserves data up to the kernel reap. "
                        "Needs `pip install psutil`. Useful for "
                        "diagnosing where peak RAM lands during a "
                        "failing run.")
    p.add_argument("--no-resume", action="store_true",
                   help="On the streamed coalescent path: ignore any "
                        "existing out/cohort/cohort.meta.json + cohort "
                        "BCFs and start a fresh simulation. Default "
                        "behaviour is to resume a prior run when its "
                        "params match — useful for multi-hour cohort "
                        "runs that get interrupted by OOM, SIGINT, or "
                        "node failure. Param mismatches surface a "
                        "clear error rather than silently re-using "
                        "incompatible state.")
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


def _run_cohort_streamed(args, chromosomes: list, rng: random.Random,
                         overlay_executor, overlay_futures: dict,
                         candidates: list) -> int:
    """Streaming cohort flow for ``--mode cohort`` (Phase 5b1).

    Replaces the in-memory cohort_sites accumulator with a chromosome-
    by-chromosome pipeline: simulate one chromosome → apply overlays
    in-place on that chunk → write to ``cohort/cohort.chr<N>.bcf`` →
    free the chunk → repeat. Peak RAM is bounded by one chromosome's
    working set rather than the whole cohort, so cohort sizes that
    OOM the in-memory path (anywhere past `n ≈ 1 000-5 000` on a
    typical workstation, depending on `--chr-length-mb`) finish
    cleanly here.

    Determinism: same `--seed` produces the same per-chromosome BCFs
    in the streamed path; rng consumption order is fixed (sample IDs
    drawn first, then overlays applied per chunk in chromosome order).
    Output at a given seed differs from the 5a in-memory cohort path
    because the overlay rng is consumed per-chunk, not globally —
    same caveat as Phase 1 for parallel chromosome simulation.

    Overlay-density semantics: ``--rsid-density 0.20`` and
    ``--clinvar-inject-density 0.01`` apply *per chromosome* in the
    streamed path rather than over the whole cohort site list. Net
    counts are roughly equal; per-chrom is arguably more correct
    biologically since chromosomes have different lengths and ClinVar
    densities.
    """
    # Resolve overlay futures up front so a malformed source fails
    # fast (before the first chromosome is simulated).
    clinvar_index: list = []
    rsid_pool: list = []
    cosmic_pool: list = []
    if overlay_futures.get("clinvar_index") is not None:
        print(f"Awaiting ClinVar overlay index for {chromosomes}...",
              file=sys.stderr)
        clinvar_index = overlay_futures["clinvar_index"].result()
        memprofile_mark(
            f"clinvar overlay index resolved ({len(clinvar_index)})")
        print(f"  {len(clinvar_index)} ClinVar pathogenic records "
              "available for overlay", file=sys.stderr)
    if overlay_futures.get("rsid_pool") is not None:
        print(f"Awaiting rsID pool...", file=sys.stderr)
        rsid_pool = overlay_futures["rsid_pool"].result()
        memprofile_mark(
            f"rsid pool resolved ({len(rsid_pool)})")
        print(f"  {len(rsid_pool)} rsID-bearing records available",
              file=sys.stderr)
    if overlay_futures.get("cosmic_pool") is not None:
        print(f"Awaiting COSMIC pool from {args.cosmic_vcf}...",
              file=sys.stderr)
        cosmic_pool = overlay_futures["cosmic_pool"].result()
        memprofile_mark(
            f"cosmic pool resolved ({len(cosmic_pool)})")
        print(f"  {len(cosmic_pool)} COSMIC records available",
              file=sys.stderr)
    if overlay_executor is not None:
        overlay_executor.shutdown(wait=True)

    workers = resolve_workers(args.workers)

    cohort_dir = args.output_dir / "cohort"

    # Resume contract (Phase 5b2): load or freshly derive the cohort-
    # identity state (sample_ids, person_seeds, per-chrom overlay
    # seeds). When an existing cohort.meta.json matches the current
    # run's params, we re-use its samples + seeds + completed-chrom
    # list rather than re-drawing from rng — so a resumed run
    # produces the same final output as the original would have on a
    # non-interrupted machine.
    try:
        resume = load_or_create_meta(
            args, chromosomes, cohort_dir, rng,
            force_fresh=args.no_resume,
        )
    except ResumeMismatch as exc:
        sys.exit(f"--resume mismatch: {exc}")
    sample_ids = resume.samples
    if resume.completed_chromosomes:
        print(
            f"  resuming run: skipping already-complete "
            f"chromosomes {resume.completed_chromosomes} "
            f"(use --no-resume to start fresh)",
            file=sys.stderr,
        )

    # demo_model gets normalised once up front because both the
    # chunk-size auto-pick (below) and the streaming loop (further
    # below) reference it.
    demo_model = (None if args.demo_model.lower() == "none"
                  else args.demo_model)

    # Phase 5f auto-pick (with the post-mortem refinements from the
    # 16 GB-host failure trace): query free RAM at run start, then
    # — when `--workers` is on auto — possibly derate the worker
    # count so the per-worker chunk size stays at or above the
    # configured floor. Without the derate, a high cpu_count host
    # forces the auto-pick to shrink chunks below 1 Mb just to fit
    # all workers in memory simultaneously, which trades useful
    # parallelism for chunk-startup overhead.
    chunk_size_mb = args.chr_chunk_mb
    if chunk_size_mb <= 0 and not args.legacy_background:
        try:
            import psutil
            available = psutil.virtual_memory().available
            # Auto-derate workers for this run if the user accepted
            # the default cpu_count and the chunk that fits at that
            # parallelism is sub-floor. An explicit --workers N
            # bypasses the derate so the user's choice is honoured.
            if args.workers == 0:
                derated = auto_derate_workers(
                    n_people=args.n, length_mb=args.chr_length_mb,
                    demo_model=demo_model, available_bytes=available,
                    requested_workers=workers,
                )
                if derated < workers:
                    print(
                        f"  --workers auto-derated from {workers} to "
                        f"{derated} to keep per-chunk RAM under the "
                        f"safety target on this host (available RAM "
                        f"{available / (1024**3):.1f} GB, n={args.n})",
                        file=sys.stderr,
                    )
                    workers = derated
            chunk_size_mb = auto_pick_chunk_size_mb(
                n_people=args.n, length_mb=args.chr_length_mb,
                demo_model=demo_model, available_bytes=available,
                workers=workers,
            )
            picked_msg = (
                f"  --chr-chunk-mb auto-picked {chunk_size_mb:.2f} Mb "
                f"(available RAM {available / (1024**3):.1f} GB, "
                f"--workers {workers}, n={args.n})"
            )
            if chunk_size_mb >= args.chr_length_mb:
                picked_msg += " — full-chromosome sim fits, no chunking"
            print(picked_msg, file=sys.stderr)
        except ImportError:
            print(
                "  --chr-chunk-mb auto-pick needs `psutil` "
                "(`pip install psutil`); falling back to "
                "full-chromosome simulation. Pass `--chr-chunk-mb N` "
                "to set a chunk size explicitly.",
                file=sys.stderr,
            )
            chunk_size_mb = 0.0
    elif chunk_size_mb > 0:
        print(
            f"  --chr-chunk-mb override: {chunk_size_mb:.2f} Mb",
            file=sys.stderr,
        )

    print(
        f"Simulating coalescent cohort (streamed): {args.n} people "
        f"across chroms {chromosomes} "
        f"(model={args.demo_model or 'uniform-constant-Ne'}, "
        f"pop={args.population}, length={args.chr_length_mb} Mb)",
        file=sys.stderr,
    )

    cohort_bcf_paths: list[Path] = []

    # Aggregate stats across chunks. SFS histograms add as Counters;
    # overlay counters add as ints.
    from collections import Counter
    sfs_total: Counter = Counter()
    overlay_stats = {
        "clinvar_annotated": 0,
        "clinvar_injected": 0,
        "rsid_injected": 0,
        "cosmic_injected": 0,
    }

    bcf_t0 = time.monotonic()
    last_progress_log = bcf_t0

    chromosomes_to_simulate = [
        c for c in chromosomes
        if not resume.is_chromosome_done(c)
    ]
    # Pre-record paths for already-complete chromosomes so the
    # manifest carries the full list whether or not we re-simulated
    # them on this invocation.
    for chrom in chromosomes:
        if resume.is_chromosome_done(chrom):
            cohort_bcf_paths.append(
                cohort_dir / f"cohort.chr{chrom}.bcf")

    memprofile_mark(
        f"streaming sim start ({len(chromosomes_to_simulate)} chroms)")
    for chrom, sites in simulate_cohort_iter(
        chromosomes=chromosomes_to_simulate, build=args.build,
        n_people=args.n, length_mb=args.chr_length_mb,
        demo_model=demo_model, population=args.population,
        rec_rate=args.rec_rate, mu=args.mu,
        rng=rng, verbose=True, workers=workers,
        chunk_size_mb=chunk_size_mb,
    ):
        memprofile_mark(f"chrom {chrom} sites yielded ({len(sites)})")
        # Each chromosome's overlays use a per-chrom rng seeded from
        # the resume record so they're independent of streaming order
        # and reproducible across resumes.
        chrom_overlay_rng = random.Random(resume.overlay_seeds[chrom])

        # Apply overlays to this chunk only. Density semantics are
        # per-chrom; reservation tracking is local to the chunk.
        if clinvar_index:
            overlay_stats["clinvar_annotated"] += annotate_clinvar(
                sites, clinvar_index)
            if args.clinvar_inject_density > 0:
                overlay_stats["clinvar_injected"] += inject_clinvar(
                    sites, clinvar_index,
                    args.clinvar_inject_density, chrom_overlay_rng,
                )
        clinvar_reserved = {i for i, s in enumerate(sites)
                            if s.get("clnsig")}
        if rsid_pool and args.rsid_density > 0:
            overlay_stats["rsid_injected"] += inject_rsids(
                sites, rsid_pool, args.rsid_density,
                chrom_overlay_rng,
                reserve_indices=clinvar_reserved,
            )
        all_overlay_reserved = {i for i, s in enumerate(sites)
                                if s.get("clnsig") or
                                s["id"].startswith("rs")}
        if cosmic_pool:
            overlay_stats["cosmic_injected"] += inject_cosmic(
                sites, cosmic_pool, args.cosmic_inject_density,
                chrom_overlay_rng,
                reserve_indices=all_overlay_reserved,
            )

        # SFS contribution from this chunk.
        sfs_total.update(sfs_histogram(sites))

        # Per-chrom BCF — `cohort.chr<N>.bcf`. The relative path goes
        # straight into the manifest's cohort_bcfs list.
        chrom_bcf = cohort_dir / f"cohort.chr{chrom}.bcf"
        # Sort within-chrom by pos for genome-order output (overlays
        # may have re-sorted but injection helpers do that too — be
        # defensive).
        sites.sort(key=lambda s: s["pos"])
        # Phase 5e Phase A: parallelise the cohort BCF write across
        # ``workers`` sample-slice writers. At workers <= 1 this
        # collapses to the original serial CohortBcfWriter path.
        # Cohort simulation itself is now serial across chromosomes
        # (the parallel-chromosome ProcessPoolExecutor path was
        # removed in this phase), so ``workers`` here controls
        # within-chromosome parallelism only.
        write_cohort_bcf_parallel(
            chrom_bcf, args.build, sample_ids, sites,
            workers=workers,
        )
        # Record the BCF in the path list and persist completion to
        # the resume record. From this point a re-run with the same
        # params will skip this chromosome.
        if chrom_bcf not in cohort_bcf_paths:
            cohort_bcf_paths.append(chrom_bcf)
        resume.mark_chromosome_done(chrom)
        memprofile_mark(f"chrom {chrom} BCF written")

        # Drop the chunk's reference so the next chromosome's
        # simulation has the previous chunk's RAM available.
        del sites

        now = time.monotonic()
        if now - last_progress_log > 20.0:
            elapsed = now - bcf_t0
            done = len(cohort_bcf_paths)
            remaining = len(chromosomes) - done
            rate_chr = done / elapsed if elapsed > 0 else 0.0
            eta = (remaining / rate_chr) if rate_chr > 0 else float("inf")
            print(
                f"  cohort BCFs: {done}/{len(chromosomes)} chromosomes "
                f"written (elapsed {_format_duration(elapsed)}, "
                f"eta {_format_duration(eta)})",
                file=sys.stderr,
            )
            last_progress_log = now

    print(
        f"  cohort sites: {sum(sfs_total.values())} alt observations "
        f"across {len(chromosomes)} chromosomes; singletons: "
        f"{sfs_total.get(1, 0)} ({singleton_fraction(sfs_total):.1%})",
        file=sys.stderr,
    )
    summary_dir = args.output_dir / "summary"
    sfs_path = summary_dir / "sfs.tsv"
    write_sfs_tsv(sfs_path, dict(sfs_total))
    print(f"  SFS histogram written to {sfs_path}", file=sys.stderr)

    if overlay_stats["clinvar_annotated"]:
        print(f"  annotated {overlay_stats['clinvar_annotated']} "
              "natural ClinVar collisions across cohort",
              file=sys.stderr)
    if overlay_stats["clinvar_injected"]:
        print(f"  injected {overlay_stats['clinvar_injected']} "
              "ClinVar pathogenic records (per-chrom density="
              f"{args.clinvar_inject_density:.3f})",
              file=sys.stderr)
    if overlay_stats["rsid_injected"]:
        print(f"  injected {overlay_stats['rsid_injected']} rsIDs "
              f"(per-chrom density={args.rsid_density:.3f})",
              file=sys.stderr)
    if overlay_stats["cosmic_injected"]:
        print(f"  injected {overlay_stats['cosmic_injected']} COSMIC "
              "records (per-chrom density="
              f"{args.cosmic_inject_density:.3f})",
              file=sys.stderr)

    manifest_mode = (
        "legacy-background" if args.legacy_background
        else "admixture-uk" if args.admixture
        else "coalescent"
    )
    cohort_bcf_rels = [
        str(p.relative_to(args.output_dir)) for p in cohort_bcf_paths
    ]
    overlay_block = None
    if not args.legacy_background:
        overlay_block = {
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

    if args.mode == "cohort":
        # Cohort-only deliverable: write the slim manifest pointing at
        # the per-chrom BCFs and stop. Per-person VCFs can be derived
        # later via `bcftools view -s SAMPLE_ID <cohort.chrN.bcf>`.
        manifest = {
            "build": args.build,
            "n_people": args.n,
            "mode": manifest_mode,
            "shape": "cohort",
            "chromosomes": chromosomes,
            "seed": args.seed,
            "samples": sample_ids,
            "cohort_bcfs": cohort_bcf_rels,
        }
        if overlay_block is not None:
            manifest["overlays"] = overlay_block
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2)
        print(
            f"--mode cohort (streamed): {len(cohort_bcf_paths)} "
            f"per-chrom BCFs written under "
            f"{cohort_dir.relative_to(args.output_dir)}/. Derive "
            f"per-person VCFs via `bcftools view -s SAMPLE_ID "
            f"<cohort.chrN.bcf>`.",
            file=sys.stderr,
        )
        print(f"Manifest written to {manifest_path}", file=sys.stderr)
        print("Done.", file=sys.stderr)
        return 0

    # --mode per-person or both: derive per-person VCFs from the
    # streamed cohort BCFs (Phase 5b2). Each worker spawns its own
    # `bcftools view -s SAMPLE | bcftools view -e 'GT="ref"'` pipeline
    # to read its assigned sample's records back, then writes a
    # per-person VCF the same way the in-memory path does (highlighted
    # variant + SVs + DP/GQ/AD + sequencing noise + truth BEDs all
    # layered identically).
    #
    # Person seeds come from the resume record, not freshly drawn
    # from rng — the master rng has been consumed by the streaming
    # loop's pre-derived per-chrom overlay seeds, so its state at
    # this point depends on the chromosome count and isn't a stable
    # starting point. The resume record's person_seeds were locked
    # in at the very start of the run.
    person_seeds = resume.person_seeds

    print(f"Writing {args.n} person VCFs into {args.output_dir} "
          f"(streamed: deriving from cohort BCFs)",
          file=sys.stderr)
    if args.error_rate > 0 or args.dropout_rate > 0:
        print(
            f"Sequencing-error model: error_rate={args.error_rate:.4f}, "
            f"dropout_rate={args.dropout_rate:.4f}",
            file=sys.stderr,
        )

    truth_dir = args.output_dir / "truth"
    contig_order = {c: i for i, c
                    in enumerate(BUILDS[args.build]["contigs"])}
    sv_chrom_span = (int(args.chr_length_mb * 1_000_000)
                     if args.chr_length_mb > 0
                     else max(BUILDS[args.build]["contigs"].values()))
    sv_chromosomes = chromosomes
    error_stats_total = new_error_stats()
    sv_total = 0

    # Workers source their per-person background by querying the
    # cohort BCFs on disk (cohort_bcfs key). The legacy / admixture
    # in-memory path uses cohort_sites instead — _person_worker checks
    # for whichever key is present. Fork-shared state inheritance
    # avoids pickling the (small) BCF path list per task.
    _PERSON_WORKER_STATE.update({
        "candidates": candidates,
        "cohort_bcfs": cohort_bcf_paths,
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
        "person_ancestry": [],   # admixture not on the streamed path
    })

    fanout_workers = resolve_workers(args.workers)
    use_pool = fanout_workers > 1 and args.n > 1
    fanout_t0 = time.monotonic()
    last_progress_log = fanout_t0
    progress_log_interval = 20.0
    memprofile_mark(
        f"per-person fanout start ({args.n} people, "
        f"{fanout_workers} workers)")

    def _maybe_log_progress(done: int) -> None:
        nonlocal last_progress_log
        now = time.monotonic()
        if (now - last_progress_log) < progress_log_interval \
                and 0 < done < args.n:
            return
        elapsed = now - fanout_t0
        rate = done / elapsed if elapsed > 0 else 0.0
        remaining = args.n - done
        eta = remaining / rate if rate > 0 else float("inf")
        print(
            f"  person VCFs: {done:,}/{args.n:,} written "
            f"({rate:.1f}/s, elapsed {_format_duration(elapsed)}, "
            f"eta {_format_duration(eta)})",
            file=sys.stderr,
        )
        last_progress_log = now

    # Phase 5g batched extraction. The legacy per-person path spawned
    # ``bcftools view -s SID`` once per (person, chrom) — at
    # ``n=3000 × 22 chroms`` that was 66,000 subprocesses each scanning
    # a few hundred MB of cohort BCF, projecting to ~38 hours. Here we
    # group sample IDs into batches and run one ``bcftools query`` per
    # chrom that emits all batch members' GT columns in a single
    # decode pass; the parent dispatches them into per-person record
    # lists, then forks a worker pool for the per-person VCF writes.
    # Batch size trades extraction-pass count against parent peak
    # memory (each batch holds B × per-person records simultaneously
    # before workers consume them).
    batch_size = max(1, args.fanout_batch_size)
    if batch_size < args.n:
        n_batches = (args.n + batch_size - 1) // batch_size
        print(
            f"  fan-out: batched extraction (batch_size={batch_size}, "
            f"{n_batches} batches, {fanout_workers} worker processes)",
            file=sys.stderr,
        )
    elif use_pool:
        print(f"  fan-out: {fanout_workers} worker processes "
              f"(fork start method, --n={args.n})",
              file=sys.stderr)

    results: list = [None] * args.n
    ctx = mp.get_context("fork") if use_pool else None
    for batch_start in range(0, args.n, batch_size):
        batch_end = min(batch_start + batch_size, args.n)
        batch_indices = range(batch_start, batch_end)
        batch_sids = sample_ids[batch_start:batch_end]
        batch_num = batch_start // batch_size + 1
        # Stage A: one bcftools subprocess per chrom decodes the
        # whole batch's columns at once. Records live on the parent
        # heap and are inherited via copy-on-write fork into each
        # worker that needs them.
        memprofile_mark(
            f"batch {batch_num} stage A start "
            f"(B={len(batch_sids)})")
        batch_backgrounds = derive_persons_batch(
            cohort_bcf_paths, batch_sids)
        _PERSON_WORKER_STATE["batch_backgrounds"] = batch_backgrounds
        memprofile_mark(
            f"batch {batch_num} stage A extracted")
        # Stage B: workers consume per-person records and write per-
        # person VCFs. A fresh pool per batch is needed so each
        # batch's children fork-inherit the right
        # ``batch_backgrounds`` snapshot — fork is millisecond-cheap,
        # and the alternative (one long-lived pool + per-task pickle
        # of the records) would dominate runtime at these list sizes.
        # The trace gets a checkpoint right after the pool spawns
        # (where COW-inherited fork attribution shows up as
        # children_rss_mb) so peak memory at the worker fork-out
        # is the easiest spot to read off the next TSV.
        if use_pool:
            with ProcessPoolExecutor(max_workers=fanout_workers,
                                     mp_context=ctx) as ex:
                futures = [
                    (i, ex.submit(_person_worker, i, sample_ids[i],
                                  person_seeds[i]))
                    for i in batch_indices
                ]
                memprofile_mark(
                    f"batch {batch_num} stage B pool spawned")
                for done_count, (i, fut) in enumerate(
                        futures, start=batch_start + 1):
                    results[i] = fut.result()
                    _maybe_log_progress(done_count)
        else:
            for i in batch_indices:
                results[i] = _person_worker(
                    i, sample_ids[i], person_seeds[i])
                _maybe_log_progress(i + 1)
        memprofile_mark(f"batch {batch_num} done")
        # Drop the parent's reference so the next batch's allocations
        # don't pile on top of the previous one.
        _PERSON_WORKER_STATE["batch_backgrounds"] = None
        del batch_backgrounds

    _PERSON_WORKER_STATE.clear()
    memprofile_mark("per-person fanout done")

    manifest_people: list = []
    for entry, person_stats, n_svs in results:
        merge_stats(error_stats_total, person_stats)
        sv_total += n_svs
        manifest_people.append(entry)
        print(_format_person_log(entry, args.n), file=sys.stderr)

    realised_fdr = (
        (error_stats_total["flipped"] + error_stats_total["dropped"]) /
        error_stats_total["total_calls"]
    ) if error_stats_total["total_calls"] else 0.0

    manifest = {
        "build": args.build,
        "n_people": args.n,
        "mode": manifest_mode,
        "shape": args.mode,
        "chromosomes": chromosomes,
        "seed": args.seed,
        "samples": sample_ids,
        "people": manifest_people,
    }
    # cohort BCFs always present in the streamed path (5b2 derives
    # per-person from them, so they exist on disk regardless of
    # whether --mode is per-person or both). The user can rm them if
    # they only want per-person; the manifest notes them so the
    # downstream pipeline can find them either way.
    manifest["cohort_bcfs"] = cohort_bcf_rels
    manifest["svs"] = {
        "per_person": args.svs_per_person,
        "length_min_bp": args.sv_length_min,
        "length_max_bp": args.sv_length_max,
        "total": sv_total,
    }
    manifest["errors"] = {
        "mode": "art" if args.art else "lightweight",
        "error_rate": args.error_rate,
        "dropout_rate": args.dropout_rate,
        "stats": dict(error_stats_total),
        "realised_fdr": round(realised_fdr, 6),
    }
    if overlay_block is not None:
        manifest["overlays"] = overlay_block

    if args.error_rate > 0 or args.dropout_rate > 0:
        print(
            f"Sequencing-error stats: "
            f"flipped={error_stats_total['flipped']}, "
            f"dropped={error_stats_total['dropped']}, "
            f"total_calls={error_stats_total['total_calls']}, "
            f"realised_fdr={realised_fdr:.4%}",
            file=sys.stderr,
        )

    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"Manifest written to {manifest_path}", file=sys.stderr)
    print("Done.", file=sys.stderr)
    return 0


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

    # Install the memory profiler before any heavy work so the TSV
    # captures the full run. The flag is opt-in; when unset, all the
    # `memprofile.mark(...)` calls scattered through this function
    # are zero-cost no-ops.
    profiler = None
    if args.profile_memory is not None:
        try:
            from .memprofile import MemoryProfiler, install
        except ImportError as exc:
            sys.exit(
                f"--profile-memory requires psutil but it failed to "
                f"import: {exc}. Install with `pip install psutil` "
                f"or drop the flag."
            )
        profiler = MemoryProfiler(args.profile_memory)
        profiler.start()
        install(profiler)
        print(
            f"  --profile-memory: sampling RSS to "
            f"{args.profile_memory} (1 s cadence + labelled marks)",
            file=sys.stderr,
        )

    print(f"Reference build: {args.build}", file=sys.stderr)
    print("Fetching ClinVar (cached across runs)...", file=sys.stderr)
    clinvar_vcf = fetch_clinvar(args.cache_dir, args.build)
    memprofile_mark("clinvar fetched")

    sig_filter = {s.strip() for s in args.clinvar_sig.split(",") if s.strip()}
    print(
        f"Loading highlighted candidates (CLNSIG in {sorted(sig_filter)})...",
        file=sys.stderr,
    )
    candidates = load_highlighted_candidates(clinvar_vcf, sig_filter)
    memprofile_mark(f"candidates loaded ({len(candidates)} records)")
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

    # Phase 5b1 — streaming cohort flow.
    # When --mode cohort is set on the standard coalescent path, divert
    # before the in-memory cohort_sites accumulator. Each chromosome's
    # sites are simulated, overlaid, written to its own per-chrom BCF,
    # and freed before the next chromosome is simulated, so peak RAM
    # stays bounded by one chromosome's working set rather than the
    # whole cohort. Legacy and admixture paths keep today's behaviour
    # because their data shapes don't fit the streaming model: the
    # legacy 1000G-pool path samples sites globally up front, and the
    # admixture path emits per-person ancestry segments alongside
    # cohort sites.
    if not args.legacy_background and not args.admixture:
        # All three --mode values now flow through the streamed
        # pipeline on the standard coalescent path. Phase 5b2 added
        # per-person derivation from the streamed cohort BCFs so
        # `--mode per-person` and `both` benefit from the streaming
        # RAM bound too — not just `--mode cohort` (which 5b1
        # delivered first). The legacy and admixture paths keep the
        # in-memory accumulator below because their data shapes
        # don't fit the streaming model (legacy samples coordinates
        # globally up front; admixture emits per-person ancestry
        # segments alongside cohort sites).
        return _run_cohort_streamed(
            args, chromosomes, rng,
            overlay_executor, overlay_futures, candidates,
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

    sample_ids = draw_sample_ids(args.n, rng)
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
        print(
            f"  person VCFs: {done:,}/{args.n:,} written "
            f"({rate:.1f}/s, elapsed {_format_duration(elapsed)}, "
            f"eta {_format_duration(eta)})",
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
    memprofile_mark("per-person fanout done")

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
