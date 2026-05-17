"""CLI entry point — wires the package together."""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import random
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
from typing import Any

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
from .bcf_writer import (
    CohortBcfWriter,
    write_cohort_bcf_parallel,
    write_cohort_bcf_parallel_from_arrow,
)
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
    # M13.3: per-person sex enables haploid emission for chrX non-PAR
    # in males, chrY non-PAR in males, chrY drops for females, and MT
    # collapses to single-allele for everyone.
    person_sexes = state.get("person_sexes")
    sex = person_sexes[i] if person_sexes else None

    # M13.3 review (Copilot): the highlighted candidate must be a
    # variant the person's per-record-ploidy will actually emit.
    # Without this filter a female could draw a chrY ClinVar
    # candidate that ``write_person_vcf`` then drops (ploidy=0) —
    # the manifest + golden BED would still list a HIGHLIGHTED row
    # for a variant that's absent from the VCF. Filter the pool to
    # ``ploidy_for(...) != 0`` for the person's sex.
    candidate_pool = candidates
    if sex is not None:
        from .builds import ploidy_for as _ploidy_for
        candidate_pool = [
            c for c in candidates
            if _ploidy_for(c["chrom"], sex, build, c["pos"]) != 0
        ]
        if not candidate_pool:
            # PR #108 review (Copilot): if every candidate is on a
            # chrom this person's ploidy will drop (e.g. every
            # candidate on chrY for a female), fail fast with a
            # clear error rather than reverting to the full list —
            # the silent fallback would reintroduce the exact bug
            # this filter was added to prevent (highlighted record
            # absent from VCF but still referenced in manifest +
            # golden BED). The user can widen --clinvar-sig or pass
            # a broader candidate source if they hit this.
            raise RuntimeError(
                f"highlighted-candidate pool empty for sex={sex!r} "
                f"on build={build!r}: every one of "
                f"{len(candidates)} ClinVar candidates lands on a "
                f"chromosome this person's ploidy drops (e.g. "
                f"every candidate on chrY for a female cohort "
                f"member). Widen --clinvar-sig or supply a broader "
                f"candidate source so at least one candidate is "
                f"emitted-by-this-sex."
            )
    hi = dict(rng.choice(candidate_pool))
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
        sex=sex,
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


# ---------------------------------------------------------------------------
# Phase 5d.1 — cohort-mode resolution + pre-flight disk check
# ---------------------------------------------------------------------------

# n threshold above which `--cohort-mode auto` picks the Arrow path. The
# Phase B sites-list path was measured-stable up to ~30k and was the
# ceiling at which CPython refcount-COW divergence started to bite; we
# pick the auto-flip well below the failure boundary.
_ARROW_AUTO_THRESHOLD = 100_000

# Rough variants-per-Mb yield for the coalescent default (mu, rec_rate,
# population). Used only for the pre-flight disk-space estimate, so a
# 2x error here just produces a slightly noisy warning, not a wrong
# decision.
_VARIANTS_PER_MB = 5_000


def _effective_chr_length_mb(
    chr_length_mb: float,
    chromosomes: list | None,
    build: str | None,
) -> float:
    """Per-chromosome length the auto-pick heuristics should use.

    ``--chr-length-mb`` semantics: ``> 0`` caps each chromosome's
    simulated prefix to that many Mb; ``0`` (default) means simulate
    the full contig length. For RAM / disk estimates the BINDING
    constraint is one chromosome's working set (we process them
    serially), so when ``chr_length_mb <= 0`` we need to look at the
    *longest* contig in the requested set — typically chr1 at ~249 Mb
    on GRCh38.

    Returns the resolved length in Mb. Falls through to ``0.0`` when
    neither an explicit length nor enough info to derive one is
    available (the caller's heuristics should then skip
    chr_length-driven checks rather than misfire on stale data).
    """
    if chr_length_mb > 0:
        return chr_length_mb
    if not chromosomes or not build:
        return 0.0
    from .builds import BUILDS
    contigs = BUILDS[build]["contigs"]
    max_bp = max(
        (contigs[c] for c in chromosomes if c in contigs),
        default=0,
    )
    return max_bp / 1_000_000 if max_bp > 0 else 0.0


def _chunking_would_split(chunk_size_mb: float,
                          effective_chr_length_mb: float) -> bool:
    """True iff ``chunk_size_mb`` would actually split the per-chrom
    simulation into multiple chunks. Mirrors the predicate inside
    :func:`coalescent.simulate_chromosome_ts` / ``simulate_chromosome``:
    a chunk size of 0 or one ≥ the chrom length is a no-op.
    """
    if chunk_size_mb <= 0:
        return False
    if effective_chr_length_mb <= 0:
        return False
    return chunk_size_mb < effective_chr_length_mb


def _estimate_materialised_parent_peak_bytes(
    n_samples: int, chr_length_mb: float
) -> int:
    """Estimate the parent-process peak RSS during the cohort-write
    phase of the *materialised* sites-list path (``sites_list`` or
    ``arrow`` mode). Used by ``--cohort-mode auto`` to decide when to
    fall through to the streaming path.

    Calibration anchor: memprof28 measured ~17 GB at n=3000 ×
    chr_length=70 Mb. That's ~17000 MB / (3000 × 70) ≈ 0.081 MB per
    (sample × Mb), i.e. ~80 KB per (sample × Mb). The dominant
    contribution is the per-site dict + sparse carriers list. Pad the
    coefficient slightly upward to be conservative — better to
    auto-pick streaming a bit early than OOM.
    """
    bytes_per_sample_per_mb = 90_000  # ~90 KB
    return int(n_samples * chr_length_mb * bytes_per_sample_per_mb)


def _check_cohort_mode_chunking_compat(
    *,
    cli_cohort_mode: str,
    resolved_cohort_mode: str,
    chunk_size_mb: float,
    chromosomes: list | None,
    build: str | None,
    chr_length_mb: float,
) -> tuple[str, str | None]:
    """Validate streaming-mode + chunking compatibility.

    Used at two cli call sites:

    - **pre-flight** (before chunk auto-pick): catches the explicit-
      explicit combination ``--cohort-mode arrow-streaming`` +
      ``--chr-chunk-mb`` that would split. Caller passes
      ``chunk_size_mb=args.chr_chunk_mb`` — i.e. the explicit value
      (0 if user defaulted, which makes ``_chunking_would_split``
      return False and this helper a no-op).
    - **second pass** (after chunk auto-pick lands): catches the
      auto-picked-chunk-would-split case. Caller passes the resolved
      ``chunk_size_mb``.

    Returns ``(final_mode, status_message)`` where ``status_message``
    is one of:

    - ``None`` — nothing to do; ``final_mode == resolved_cohort_mode``.
    - prefix ``"ERROR: "`` — fatal; caller should ``sys.exit`` the
      message body. Fires when user explicitly asked for
      ``--cohort-mode arrow-streaming`` AND chunking would split.
    - prefix ``"INFO: "`` — informational demote; caller should
      print to stderr. Fires when ``--cohort-mode auto`` picked
      streaming but chunking would split, and the helper falls
      back to ``arrow`` (which preserves mmap-share for workers).

    Pulled out of cli's main flow so the user-facing guard rails
    around streaming-+-chunking are unit-testable without needing
    the full ``cli.main`` integration test stack. PR #55 Copilot
    review caught the gap.
    """
    if resolved_cohort_mode != "arrow-streaming":
        return (resolved_cohort_mode, None)
    eff_len_mb = _effective_chr_length_mb(
        chr_length_mb, chromosomes, build,
    )
    if not _chunking_would_split(chunk_size_mb, eff_len_mb):
        return (resolved_cohort_mode, None)
    # Streaming + chunking-would-split is unsupported. Branch on
    # whether the user picked streaming explicitly or via auto.
    if cli_cohort_mode == "arrow-streaming":
        # NB: ``--chr-chunk-mb 0`` is the auto-pick sentinel, NOT a
        # no-chunking flag, so we can't recommend it as a remedy —
        # an explicit-streaming user who passes 0 would just trip
        # this same check on whatever the auto-picker lands at.
        # Recommend an explicit chunk size at or above the effective
        # chrom length (which guarantees no split), or switching
        # modes. ``ceil`` (not ``round`` or ``.0f``) because a
        # recommendation of ``round(70.1) = 70`` would still split
        # a 70.1 Mb chromosome, and ``round(0.4) = 0`` would
        # reintroduce the forbidden ``--chr-chunk-mb 0`` suggestion
        # for sub-Mb effective lengths. PR #60 review caught this.
        # We reach this branch only when
        # ``_chunking_would_split`` returned True, which requires
        # ``eff_len_mb > 0``, so ``ceil`` always yields >= 1.
        suggested_chunk = math.ceil(eff_len_mb)
        return (
            resolved_cohort_mode,
            f"ERROR: --cohort-mode arrow-streaming does not yet "
            f"support a chunked simulation that splits chromosomes "
            f"(chunk_size_mb={chunk_size_mb:.2f} Mb < effective "
            f"chrom length {eff_len_mb:.1f} Mb). Pass an explicit "
            f"--chr-chunk-mb {suggested_chunk} (or larger) to keep "
            f"the chromosome unsplit, or use --cohort-mode arrow / "
            f"sites_list.",
        )
    return (
        "arrow",
        f"INFO: --cohort-mode auto demoted arrow-streaming → arrow "
        f"because chunk_size_mb={chunk_size_mb:.2f} would split the "
        f"simulation (streaming does not yet support chunked sim).",
    )


def _resolve_cohort_mode(cli_mode: str, n: int,
                         chr_length_mb: float = 0.0,
                         chromosomes: list | None = None,
                         build: str | None = None,
                         chunk_size_mb: float = 0.0,
                         host_ram_bytes: int | None = None) -> str:
    """Map ``--cohort-mode {sites_list,arrow,arrow-streaming,auto}``
    + ``--n`` + chromosome size + chunking to the concrete mode used
    by the cohort-write loop.

    Returns one of ``"sites_list"``, ``"arrow"``, or
    ``"arrow-streaming"`` (never ``"auto"``).

    ``auto`` resolution:

    - First, check whether the predicted materialised parent peak
      exceeds ~50% of host RAM. If so, *and* chunking would not
      actually split the simulation (streaming doesn't yet support
      chunked sim), pick ``arrow-streaming``.
    - Otherwise, pick ``arrow`` for ``n >= 100000``.
    - Otherwise, pick ``sites_list`` (Phase B; bounded at n~30k by
      the refcount-COW divergence).

    For the predicted-peak check the effective per-chromosome length
    is needed. When ``chr_length_mb > 0`` use it directly; when it's
    0 (= full contig length), derive the longest contig from
    ``(chromosomes, build)`` since we process chromosomes serially
    and chr1 is typically the binding constraint.

    ``chunk_size_mb`` is consulted only for the streaming branch:
    streaming doesn't yet support a chunked simulation that would
    actually split a chromosome, so auto avoids streaming when
    chunking would split. The cli still has a separate explicit-
    mode check that hard-errors on
    ``--cohort-mode arrow-streaming --chr-chunk-mb > 0`` when the
    chunk would split.
    """
    if cli_mode == "sites_list":
        return "sites_list"
    if cli_mode == "arrow":
        return "arrow"
    if cli_mode == "arrow-streaming":
        return "arrow-streaming"
    # auto
    if host_ram_bytes is None:
        try:
            import psutil
            host_ram_bytes = psutil.virtual_memory().total
        except ImportError:
            host_ram_bytes = 32 * 1024**3  # 32 GB conservative fallback
    effective_length_mb = _effective_chr_length_mb(
        chr_length_mb, chromosomes, build,
    )
    if effective_length_mb > 0:
        predicted = _estimate_materialised_parent_peak_bytes(
            n, effective_length_mb,
        )
        if predicted > 0.5 * host_ram_bytes:
            # Streaming wins on memory — *unless* chunking would
            # split the simulation, in which case the chunked
            # materialised+arrow path is the right call (chunked
            # streaming isn't implemented yet). When chunking
            # blocks streaming we still prefer ``arrow`` over
            # ``sites_list`` regardless of the n>=100k threshold:
            # the workload that triggered the predicted-peak gate
            # is large enough that workers benefit from mmap-share
            # even if parent still materialises the sites list.
            if not _chunking_would_split(
                chunk_size_mb, effective_length_mb,
            ):
                return "arrow-streaming"
            return "arrow"
    return "arrow" if n >= _ARROW_AUTO_THRESHOLD else "sites_list"


def _estimate_arrow_chrom_scratch_bytes(
    n_samples: int, chr_length_mb: float,
    cohort_mode: str = "arrow",
) -> int:
    """Estimate per-chromosome scratch disk for the Arrow intermediate
    (and, for the ``arrow-streaming`` mode only, the Fix B.1 carriers
    sidecar).

    Used by :func:`_preflight_arrow_disk_check`. Two terms:

    1. **Arrow file** (all arrow modes): ``n_haplotypes ×
       variants_per_chrom × 1 byte`` (dense int8 genotypes; INFO +
       offsets add ~5 %). This is the on-disk size of the cohort
       intermediate after Arrow's IPC writer encodes it.
    2. **Carriers sidecar** (``arrow-streaming`` only): packed
       ``np.int32 × 2`` per non-zero haplotype. Each site averages
       ~θ × log(n) carriers under coalescent SFS — empirically sums
       to roughly the same byte total as the Arrow file at canonical
       AFs. Budget 1× the Arrow term as a conservative estimate.

    Only the ``arrow-streaming`` path creates a sidecar; the
    materialised ``arrow`` mode keeps carriers in the in-RAM sites
    list and never spills. Adding the sidecar term unconditionally
    would falsely fail tight-disk runs that use ``--cohort-mode
    arrow``. PR #77 review caught this.

    Rough accuracy is OK — the check produces a warning when free
    disk is tight, not a hard refusal.
    """
    variants = max(1, int(chr_length_mb * _VARIANTS_PER_MB))
    arrow_raw = 2 * n_samples * variants
    arrow_with_overhead = int(arrow_raw * 1.05)
    if cohort_mode == "arrow-streaming":
        # Sidecar ~same order of magnitude as the Arrow file at
        # canonical AFs; budget 1× as a conservative estimate.
        return arrow_with_overhead + arrow_with_overhead
    return arrow_with_overhead


def _preflight_arrow_disk_check(
    cohort_dir: Path, n_samples: int, chr_length_mb: float,
    cohort_mode: str = "arrow",
) -> None:
    """Warn (or fail) if the cohort directory's filesystem doesn't have
    enough free space to hold one chromosome's Arrow scratch + the
    final BCF concurrently.

    The Arrow file is created and deleted per chromosome so the steady-
    state requirement is one Arrow file plus the BCF being written. We
    require ``2x`` the per-chrom estimate free; below ``1x`` we fail
    with a clear message because the Arrow write itself wouldn't
    complete.

    ``cohort_mode`` is plumbed through to
    :func:`_estimate_arrow_chrom_scratch_bytes` so the sidecar term
    is only added for ``arrow-streaming`` runs (PR #77 review).
    """
    cohort_dir = Path(cohort_dir)
    cohort_dir.mkdir(parents=True, exist_ok=True)
    per_chrom = _estimate_arrow_chrom_scratch_bytes(
        n_samples, chr_length_mb, cohort_mode=cohort_mode,
    )
    try:
        free = shutil.disk_usage(cohort_dir).free
    except OSError:
        # Best-effort — if disk_usage fails we don't block the run.
        return

    # Sidecar-aware scratch note for arrow-streaming users — the
    # ~2× budget vs plain ``arrow`` is intentional, not a bug.
    if cohort_mode == "arrow-streaming":
        scratch_note = (
            " (Arrow IPC file + carriers sidecar)"
        )
    else:
        scratch_note = ""

    if free < per_chrom:
        raise SystemExit(
            f"--cohort-mode {cohort_mode} needs "
            f"~{per_chrom / 1e9:.1f} GB scratch{scratch_note} per "
            f"chromosome (n={n_samples}, length={chr_length_mb} Mb); "
            f"only {free / 1e9:.1f} GB free under {cohort_dir}. "
            f"Free up disk or pass --cohort-mode sites_list "
            f"(supported up to n~30000)."
        )
    if free < 2 * per_chrom:
        print(
            f"  WARNING: --cohort-mode {cohort_mode} estimates "
            f"~{per_chrom / 1e9:.1f} GB scratch{scratch_note} per "
            f"chromosome and only {free / 1e9:.1f} GB is free under "
            f"{cohort_dir}. The first chromosome will fit, but the "
            f"final BCF + Arrow file may not coexist. Free up disk "
            f"or set `--cohort-mode sites_list` to be safe.",
            file=sys.stderr,
        )


def _write_cohort_chrom_bcf(
    mode: str,
    chrom: str,
    chrom_bcf: Path,
    sites: list,
    sample_ids: list,
    build: str,
    workers: int,
    cohort_dir: Path,
    arrow_batch_size: int,
) -> None:
    """Dispatch one chromosome's cohort BCF write between the
    sites-list (Phase B) and Arrow (Phase 5d.1) paths.

    Arrow path: writes ``cohort_dir / .arrow / cohort.chr<N>.arrow``,
    runs the parallel BCF writers against it, deletes the Arrow file
    on success. On failure the Arrow file is left in place for
    postmortem (matches the ``.partials/`` cleanup behaviour).
    """
    if mode in ("arrow", "arrow-streaming"):
        arrow_dir = cohort_dir / ".arrow"
        arrow_dir.mkdir(parents=True, exist_ok=True)
        arrow_path = arrow_dir / f"cohort.chr{chrom}.arrow"
        from .cohort_arrow import write_arrow_file
        # ``sites`` here is either a materialised list (mode=="arrow")
        # or a streaming generator (mode=="arrow-streaming"). Both
        # iterate identically; ``iter()`` is a no-op on generators.
        write_arrow_file(
            arrow_path, chrom, len(sample_ids), iter(sites),
            batch_size=arrow_batch_size,
        )
        # Clean up the Arrow file only on success — on failure leave it
        # in place so a postmortem can mmap it and inspect what the
        # writer was about to consume.
        write_cohort_bcf_parallel_from_arrow(
            arrow_path, chrom_bcf, build, sample_ids, workers,
        )
        if arrow_path.exists():
            arrow_path.unlink()
        try:
            arrow_dir.rmdir()
        except OSError:
            pass
        return

    write_cohort_bcf_parallel(
        chrom_bcf, build, sample_ids, sites, workers=workers,
    )


# ---------------------------------------------------------------------------
# Per-chromosome iterators (Phase 5d.1 streaming refactor — PR 3)
# ---------------------------------------------------------------------------
#
# Two iterators with the same ``(chrom, sites)`` contract so the cli's
# main cohort-phase loop can stay simple. ``sites`` is either a fully-
# materialised list (legacy path) or a streaming generator that the
# Arrow writer drains directly (PR 3 path). Both versions handle:
#
#   - msprime simulation (per chrom)
#   - the four overlays (annotate_clinvar, inject_clinvar, inject_rsids,
#     inject_cosmic) with per-chrom rng derived from the resume record
#   - SFS histogram accumulation
#
# What differs:
#
#   - Materialised: builds the full sites list, applies overlays in
#     place, runs sfs_histogram on the list, sorts by pos.
#   - Streaming: builds only the metadata + overlay plans, returns a
#     generator that pulls full sites from the tree on demand. SFS
#     accumulation runs inline (per yielded site).


def _resolve_reference_fasta(
    args, chromosomes: list,
) -> tuple[Any, "Path | None"]:
    """Resolve ``--reference-fasta`` / ``--no-reference-fasta`` /
    auto-fetch into an ``(open_handle, path)`` pair.

    Three modes, in priority order:

    1. ``--no-reference-fasta`` flag set → returns ``(None, None)``;
       caller falls back to fabricated REF (legacy
       ``rng.choice("ACGT")`` path).
    2. ``--reference-fasta <path>`` set → opens that path.
    3. Otherwise (the default) → auto-discovers / auto-downloads
       to ``<cache_dir>/reference/<build>.fa`` via
       :func:`reference.fetch_reference_fasta`. Mirrors the
       ClinVar caching pattern: first run downloads + indexes,
       subsequent runs find the cache and return immediately.

    The two-tuple return shape exists because the streaming
    coalescent path keeps the open handle (in-process use, fastest)
    while the admixture path needs the path string (its
    ``ProcessPoolExecutor`` workers re-open from the path because
    ``FastaFile`` doesn't pickle).

    Pre-flight ``validate_fasta`` runs against the requested
    chromosomes + ``chr_length_mb`` so a missing chrom / too-short
    contig fails at startup, not after msprime has been running
    for hours.
    """
    # PR #86 review (Copilot): --no-reference-fasta + an explicit
    # --reference-fasta path is contradictory. Reject loudly rather
    # than silently letting the opt-out win.
    if (
        getattr(args, "no_reference_fasta", False)
        and args.reference_fasta is not None
    ):
        sys.exit(
            "--no-reference-fasta and --reference-fasta are mutually "
            "exclusive: --no-reference-fasta opts out of real REF "
            "entirely, --reference-fasta <path> opts in to a specific "
            "FASTA. Pick one.",
        )

    if getattr(args, "no_reference_fasta", False):
        print(
            "  --no-reference-fasta: skipping reference FASTA; "
            "emitted REF will be drawn from rng.choice('ACGT') "
            "(legacy fabricated-REF path — VCFs will not pass "
            "`bcftools norm --check-ref`)",
            file=sys.stderr,
        )
        return None, None

    # PR #86 + #87 review (Copilot): default-on means non-bioinformatician
    # users hit this path too. Every reference-resolution failure mode —
    # auto-fetch (``fetch_reference_fasta``) AND open (``load_fasta``) —
    # needs an actionable message, not a bare stack trace. Covers:
    #   - ImportError (no pysam; raised by both fetch and load)
    #   - FileNotFoundError (stale --reference-fasta path)
    #   - ValueError (unknown build, missing reference_fasta_url,
    #     unindexed/unreadable FASTA from load_fasta)
    try:
        if args.reference_fasta is not None:
            fasta_path = Path(args.reference_fasta)
        else:
            # Auto-discover / auto-fetch. ``fetch_reference_fasta`` is
            # the cached form: re-runs against the same ``cache_dir``
            # are a no-op once the FASTA is on disk.
            print(
                f"  --reference-fasta unset; resolving cached FASTA "
                f"at <cache-dir>/reference/{args.build}.fa "
                f"(pass --no-reference-fasta to skip)",
                file=sys.stderr,
            )
            from .reference import fetch_reference_fasta
            fasta_path = fetch_reference_fasta(args.cache_dir, args.build)
        from .reference import load_fasta, validate_fasta
        print(f"  loading reference FASTA: {fasta_path}", file=sys.stderr)
        fa = load_fasta(fasta_path)
    except (ImportError, FileNotFoundError, ValueError) as exc:
        sys.exit(
            f"--reference-fasta: {exc}\n"
            "Hint: install pysam (`pip install pysam`) or rerun with "
            "`--no-reference-fasta` to fall back to the legacy "
            "fabricated-REF path.",
        )
    # PR #86 review (claude review): validate_fasta raise must close
    # the just-opened handle before sys.exit so test harnesses that
    # invoke main() many times don't leak file descriptors.
    try:
        validate_fasta(fa, chromosomes, args.chr_length_mb, args.build)
    except ValueError as exc:
        fa.close()
        sys.exit(str(exc))
    print(
        f"  reference FASTA OK ({len(fa.references)} contigs)",
        file=sys.stderr,
    )
    return fa, fasta_path


def _iter_chrom_sites_materialised(
    *, chromosomes_to_simulate, args, demo_model, rng, workers,
    resume, clinvar_index, rsid_pool, cosmic_pool,
    overlay_stats, sfs_total, chunk_size_mb, fasta=None,
):
    """Legacy materialised per-chromosome iterator. Yields
    ``(chrom, sites_list)`` after applying overlays + sorting in place.

    ``fasta`` (an open ``pysam.FastaFile`` handle, or ``None``)
    forwards M12's reference-aware REF picking to the inner
    producers — when present, REF bases come from the FASTA;
    when ``None`` the legacy ``rng.choice('ACGT')`` path applies.
    """
    from .coalescent import simulate_cohort_iter
    for chrom, sites in simulate_cohort_iter(
        chromosomes=chromosomes_to_simulate, build=args.build,
        n_people=args.n, length_mb=args.chr_length_mb,
        demo_model=demo_model, population=args.population,
        rec_rate=args.rec_rate, mu=args.mu,
        rng=rng, verbose=True, workers=workers,
        chunk_size_mb=chunk_size_mb,
        fasta=fasta,
    ):
        memprofile_mark(f"chrom {chrom} sites yielded ({len(sites)})")
        chrom_overlay_rng = random.Random(resume.overlay_seeds[chrom])
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
        sfs_total.update(sfs_histogram(sites))
        # Sort within-chrom by pos (overlays re-sort, but be defensive).
        sites.sort(key=lambda s: s["pos"])
        yield chrom, sites


def _iter_chrom_sites_streaming(
    *, chromosomes_to_simulate, args, demo_model, rng,
    resume, clinvar_index, rsid_pool, cosmic_pool,
    overlay_stats, sfs_total, chunk_size_mb, fasta=None,
):
    """Streaming per-chromosome iterator (Phase 5d.1 PR 3). Yields
    ``(chrom, sites_generator)`` where ``sites_generator`` is the
    output of :func:`coalescent.stream_cohort_sites` — sites flow
    through the Arrow writer without ever being materialised in a
    Python list. Parent peak drops from ~17 GB to a few hundred MB
    at WGS scale.

    Byte-identical to the materialised iterator at every fixed seed
    combination (master rng + per-chrom overlay rng) — see the PR 2
    parity tests in ``test_streaming_cohort.py``.

    SFS histogram accumulation happens inline (a tee-style generator
    wraps the stream so each yielded site updates ``sfs_total`` once).
    Overlay counts are accumulated against ``overlay_stats`` once
    per chromosome, before the generator is consumed, by inspecting
    the plans before pass 2 runs.
    """
    from .coalescent import (
        simulate_cohort_ts_iter,
        _tree_sequence_to_sites_meta,
        _build_annotation_map,
        _plan_cohort_overlays,
        _stream_cohort_pass2,
    )
    from .carriers_sidecar import CarriersSidecar
    from .titv import DEFAULT_TARGET_TITV

    # The carriers spill sidecar lives next to the per-chrom Arrow
    # scratch file. Resolving the cohort dir from args mirrors the
    # cohort BCF path resolution elsewhere in the cli. Spill is
    # enabled unconditionally for the streaming path because the
    # disk-IO overhead at small n (~few ms/yield) is dwarfed by the
    # parent-RSS ceiling it removes at WGS scale — see Fix B.1 in
    # PERFORMANCE_BUDGETS.md § "Known scaling ceiling".
    spill_dir = Path(args.output_dir) / "cohort" / "carriers_scratch"

    for chrom, ts, chrom_rng in simulate_cohort_ts_iter(
        chromosomes=chromosomes_to_simulate, build=args.build,
        n_people=args.n, length_mb=args.chr_length_mb,
        demo_model=demo_model, population=args.population,
        rec_rate=args.rec_rate, mu=args.mu,
        rng=rng, verbose=True, chunk_size_mb=chunk_size_mb,
    ):
        memprofile_mark(f"chrom {chrom} ts ready ({ts.num_sites} sites)")
        chrom_overlay_rng = random.Random(resume.overlay_seeds[chrom])

        # Pre-compute the overlay plans once per chrom (so overlay_stats
        # can be incremented before pass 2 drains the generator) while
        # bracketing chrom_rng across the meta walk so pass 2's rng
        # reproduces the same ref/alt picks. The chained-overlay
        # planner mirrors the materialised path's sort-between-injects
        # semantics for byte-identical output — see
        # :func:`coalescent._plan_cohort_overlays`.
        state_before = chrom_rng.getstate()
        sites_meta = _tree_sequence_to_sites_meta(
            ts, chrom, args.n, chrom_rng, DEFAULT_TARGET_TITV,
            fasta=fasta,
        )
        annotation_map = _build_annotation_map(
            sites_meta, clinvar_index,
        )
        clinvar_plan, rsid_plan, cosmic_plan = _plan_cohort_overlays(
            sites_meta, annotation_map,
            clinvar_records=clinvar_index,
            clinvar_inject_density=args.clinvar_inject_density,
            rsid_pool=rsid_pool, rsid_density=args.rsid_density,
            cosmic_records=cosmic_pool,
            cosmic_inject_density=args.cosmic_inject_density,
            overlay_rng=chrom_overlay_rng,
        )
        overlay_stats["clinvar_annotated"] += len(annotation_map)
        overlay_stats["clinvar_injected"] += len(clinvar_plan)
        overlay_stats["rsid_injected"] += len(rsid_plan)
        overlay_stats["cosmic_injected"] += len(cosmic_plan)
        chrom_rng.setstate(state_before)

        inject_map: dict = {}
        inject_map.update(clinvar_plan)
        inject_map.update(rsid_plan)
        inject_map.update(cosmic_plan)

        # Per-chrom carriers sidecar — spilled to disk so the
        # safe-yield heap holds (offset, length) refs instead of
        # full per-site carriers arrays. Cleanup via try/finally so
        # the sidecar is unlinked even when the downstream Arrow
        # writer raises mid-chrom; the sidecar's own ``close()`` is
        # idempotent.
        sidecar = CarriersSidecar(
            spill_dir / f"carriers.chr{chrom}.spill",
        )
        try:
            stream = _stream_cohort_pass2(
                ts, chrom, args.n, chrom_rng, DEFAULT_TARGET_TITV,
                sites_meta, inject_map, annotation_map,
                carriers_sidecar=sidecar,
                fasta=fasta,
            )

            # Tee the stream so we update sfs_total per site as it
            # flows through to the Arrow writer.
            def _tee_sfs(site_iter, sfs_total_ref):
                for site in site_iter:
                    for ac in site.get("acs") or []:
                        sfs_total_ref[ac] = sfs_total_ref.get(ac, 0) + 1
                    yield site

            yield chrom, _tee_sfs(stream, sfs_total)
        finally:
            sidecar.close()


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


def _male_fraction_arg(value: str) -> float:
    """Argparse ``type=`` for ``--male-fraction``.

    Mirrors the Pydantic ``CohortConfig.male_fraction`` constraint
    (``ge=0.0, le=1.0``) so the CLI rejects out-of-range and non-
    finite inputs (``-1``, ``2``, ``nan``, ``inf``) instead of
    funnelling them into ``_draw_sexes`` and silently producing
    all-male / all-female / NaN-degenerate cohorts.
    """
    try:
        f = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(
            f"--male-fraction must be a float in [0.0, 1.0]; "
            f"got {value!r}"
        ) from exc
    if not math.isfinite(f) or not 0.0 <= f <= 1.0:
        raise argparse.ArgumentTypeError(
            f"--male-fraction must be a finite float in [0.0, 1.0]; "
            f"got {value!r}"
        )
    return f


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
    p.add_argument("--config", type=str, default=None,
                   help="Load values from this YAML config file. CLI "
                        "flags still override config values; config "
                        "values still override built-in defaults. If "
                        "omitted, the tool looks for "
                        "`generate_people_config.yaml` in the current "
                        "directory and uses it automatically when "
                        "present. Pass --no-config to disable that "
                        "auto-discovery.")
    p.add_argument("--no-config", action="store_true",
                   help="Skip the auto-discovery of "
                        "`generate_people_config.yaml` in the current "
                        "directory. Has no effect when --config is set.")
    p.add_argument("--print-config", action="store_true",
                   help="Print a starter `generate_people_config.yaml` "
                        "to stdout (every field at its built-in "
                        "default, with a leading comment per field) "
                        "and exit. Redirect to a file to bootstrap "
                        "the YAML config: "
                        "`python synthetic_people/generate_people.py "
                        "--print-config > "
                        "generate_people_config.yaml`. The emitted "
                        "config is a valid no-op as-is — running with "
                        "it changes nothing — so edit just the values "
                        "you want to change. See TUTORIAL.md §10.")
    p.add_argument("--n", type=int, default=10,
                   help="Cohort size: number of person VCFs to generate")
    p.add_argument("--output-dir", type=Path,
                   default=script_dir / "out",
                   help="Where to write person_<N>.vcf.gz")
    p.add_argument("--cache-dir", type=Path,
                   default=script_dir / "cache",
                   help="Cache directory: holds the ClinVar VCF "
                        "(~70 MB) plus, since M12 (2026-05-14), the "
                        "build's reference FASTA (~3 GB decompressed "
                        "under <cache-dir>/reference/<build>.fa). "
                        "Both are downloaded on first run and reused "
                        "thereafter. Pass --no-reference-fasta to "
                        "skip the FASTA fetch.")
    p.add_argument("--build", choices=list(BUILDS), default="GRCh38",
                   help="Reference build; must match background VCFs")
    p.add_argument("--reference-fasta", type=Path, default=None,
                   help="[M12] Path to the reference FASTA matching "
                        "--build (e.g. GRCh38 primary assembly). When "
                        "set, REF bases come from this FASTA instead "
                        "of `rng.choice('ACGT')`, so the emitted VCFs "
                        "validate cleanly with `bcftools norm "
                        "--check-ref`. The FASTA must have a `.fai` "
                        "index (run `samtools faidx <fasta>` once if "
                        "missing). When unset (the default), the cli "
                        "auto-discovers a cached FASTA at "
                        "`<cache-dir>/reference/<build>.fa`; if absent "
                        "it downloads the Ensembl primary assembly "
                        "(~900 MB compressed → ~3 GB decompressed) "
                        "into that path on first run. Pass "
                        "--no-reference-fasta to opt out of both the "
                        "auto-discover and the auto-download.")
    p.add_argument("--no-reference-fasta", action="store_true",
                   help="[M12] Skip the reference-FASTA "
                        "auto-discover / auto-download. Emitted REF "
                        "is then drawn from `rng.choice('ACGT')` "
                        "(legacy fabricated-REF path), so the output "
                        "won't pass `bcftools norm --check-ref`. "
                        "Useful for smoke tests / dev runs that "
                        "don't need real REF and don't want to pay "
                        "the cache cost.")
    p.add_argument("--seed", type=int, default=None,
                   help="RNG seed for deterministic output. Omit for "
                        "different people each run.")
    p.add_argument("--male-fraction", type=_male_fraction_arg,
                   default=0.5,
                   help="[M13] Probability each person is drawn as "
                        "male. 0.5 (default) = balanced cohort; 0.2 = "
                        "~20%% male, ~80%% female; 0.8 = ~80%% male. "
                        "Must be a finite float in [0.0, 1.0]. "
                        "Mirrors the YAML field `cohort.male_fraction`. "
                        "Per-person sex assignment is recorded at the "
                        "top level of `manifest.json` as `sex: [\"m\", "
                        "\"f\", ...]`, parallel-indexed to `samples`.")
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
                        "noise model. Currently rejected with a "
                        "clear message: M12 ships the reference "
                        "FASTA cache so the dep is in place, but "
                        "the ART pipeline itself (art_illumina + "
                        "bcftools call wiring) hasn't been hooked "
                        "up yet.")
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
    p.add_argument("--cohort-mode",
                   choices=("sites_list", "arrow", "arrow-streaming",
                            "auto"),
                   default="auto",
                   help="[perf, Phase 5d.1] Cohort intermediate "
                        "between simulation and the BCF write. "
                        "`sites_list` (Phase B) keeps the cohort in "
                        "RAM and fork-shares to workers; works up to "
                        "n~30k. `arrow` materialises the sites list "
                        "and streams it to an Arrow IPC scratch file "
                        "per chromosome; workers mmap-read. "
                        "`arrow-streaming` (option 3) skips the "
                        "materialised sites list entirely — parent "
                        "streams directly from the msprime tree "
                        "sequence to the Arrow file, dropping parent "
                        "peak from ~17 GB to a few hundred MB at WGS "
                        "scale. `auto` picks `arrow-streaming` when "
                        "predicted parent peak > 50%% host RAM, "
                        "`arrow` for --n>=100000 below that, "
                        "`sites_list` for smaller cohorts. The Arrow "
                        "paths need `pip install pyarrow`; without "
                        "it, `auto` falls back to `sites_list` with a "
                        "warning.")
    p.add_argument("--cohort-arrow-batch-size", type=int, default=256,
                   help="[perf, Phase 5d.1] Sites-per-batch when "
                        "streaming the cohort Arrow IPC file. The "
                        "default of 256 was identified by Spike 2b "
                        "as the knee of the parent-RSS / write-"
                        "throughput trade-off (predicted parent peak "
                        "~5 GB at n=1M; +7%% throughput cost vs the "
                        "fastest cell). Only consulted in "
                        "`--cohort-mode arrow`.")
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

    # Phase 5d.1 — resolve the cohort intermediate mode for this run
    # and (if it ends up Arrow) pre-flight the disk before any
    # simulation work. Doing this before resume / sample draw means
    # an insufficient-disk run aborts cheaply.
    #
    # Auto-pick considers full-contig length (when --chr-length-mb 0)
    # by deriving the longest contig in the requested chromosome set
    # — chr1 at ~249 Mb on GRCh38 is what drives the WGS-n=3000
    # streaming pick. Chunking is also taken into account because
    # streaming doesn't yet support a chunked simulation that would
    # split a chromosome; in that case auto stays on arrow.
    cohort_mode = _resolve_cohort_mode(
        args.cohort_mode, args.n,
        chr_length_mb=args.chr_length_mb,
        chromosomes=chromosomes, build=args.build,
        chunk_size_mb=args.chr_chunk_mb,
    )
    # Pre-flight: explicit ``--cohort-mode arrow-streaming`` combined
    # with an explicit ``--chr-chunk-mb`` that would split the
    # simulation isn't supported (simulate_chromosome_ts raises
    # NotImplementedError if it gets there). Catch it here so the
    # user sees a clear message immediately rather than after several
    # minutes of msprime work. Passing ``args.chr_chunk_mb`` here
    # (default 0) makes the helper a no-op for auto users — the
    # auto-picked-chunk case is caught by the second-pass check below.
    _, _msg = _check_cohort_mode_chunking_compat(
        cli_cohort_mode=args.cohort_mode,
        resolved_cohort_mode=cohort_mode,
        chunk_size_mb=args.chr_chunk_mb,
        chromosomes=chromosomes,
        build=args.build,
        chr_length_mb=args.chr_length_mb,
    )
    if _msg and _msg.startswith("ERROR: "):
        sys.exit(_msg[len("ERROR: "):])
    if cohort_mode in ("arrow", "arrow-streaming"):
        try:
            from . import cohort_arrow  # noqa: F401
        except ImportError as exc:
            if args.cohort_mode == "auto":
                print(
                    f"  WARNING: --cohort-mode auto chose "
                    f"`{cohort_mode}` for --n={args.n}, but pyarrow "
                    f"is not installed; falling back to `sites_list`. "
                    f"Install with `pip install pyarrow` to use the "
                    f"Arrow path. ({exc})",
                    file=sys.stderr,
                )
                cohort_mode = "sites_list"
            else:
                sys.exit(
                    f"--cohort-mode {cohort_mode} requires pyarrow; "
                    f"install with `pip install pyarrow` or pass "
                    f"--cohort-mode sites_list (supported up to "
                    f"n~30000)."
                )
    if cohort_mode in ("arrow", "arrow-streaming"):
        # Use the EFFECTIVE chrom length, not the raw flag value:
        # when --chr-length-mb is 0 (full contig), the raw value would
        # make the scratch-bytes estimate collapse to ~0, silently
        # turning the disk pre-flight into a no-op for exactly the
        # WGS-scale runs that need it most. PR #55 Copilot review
        # caught this; resolve to chr1-equivalent here too.
        # ``cohort_mode`` is passed so the sidecar term in the
        # scratch estimate is only added for ``arrow-streaming``
        # — materialised ``arrow`` doesn't spill (PR #77 review).
        _preflight_arrow_disk_check(
            cohort_dir, args.n,
            _effective_chr_length_mb(
                args.chr_length_mb, chromosomes, args.build,
            ),
            cohort_mode=cohort_mode,
        )
    print(
        f"  cohort intermediate mode: {cohort_mode} "
        f"(--cohort-mode={args.cohort_mode}, --n={args.n})",
        file=sys.stderr,
    )

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

    # Second-pass cohort-mode adjustment: if the resolved
    # ``chunk_size_mb`` (auto-picked or explicit) would actually
    # split the simulation, and we previously resolved cohort_mode to
    # arrow-streaming, we have to handle the unsupported combination.
    # The helper routes both paths: ERROR for explicit
    # ``--cohort-mode arrow-streaming``, INFO + demote → arrow for
    # ``--cohort-mode auto``.
    cohort_mode, _msg = _check_cohort_mode_chunking_compat(
        cli_cohort_mode=args.cohort_mode,
        resolved_cohort_mode=cohort_mode,
        chunk_size_mb=chunk_size_mb,
        chromosomes=chromosomes,
        build=args.build,
        chr_length_mb=args.chr_length_mb,
    )
    if _msg and _msg.startswith("ERROR: "):
        sys.exit(_msg[len("ERROR: "):])
    elif _msg and _msg.startswith("INFO: "):
        print(f"  {_msg[len('INFO: '):]}", file=sys.stderr)

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

    # M12: resolve the reference FASTA. The helper handles three
    # modes (--no-reference-fasta, --reference-fasta <path>, and
    # auto-fetch via cache_dir/reference/<build>.fa); see
    # ``_resolve_reference_fasta`` for the contract. The streaming
    # path uses the in-process handle directly; ``fasta_path`` is
    # threaded through to the materialised iterator's fasta arg
    # for parity but isn't used downstream of the open handle.
    fasta, _fasta_path = _resolve_reference_fasta(
        args, chromosomes_to_simulate,
    )

    memprofile_mark(
        f"streaming sim start ({len(chromosomes_to_simulate)} chroms)")
    # PR #86 review (Copilot + claude review): close the FASTA handle
    # on every exit path. The iterator helpers below pass ``fasta``
    # into generators that hold a reference until they're exhausted,
    # so it's only safe to close after the loop completes (or unwinds
    # via exception). Tests run ``main()`` repeatedly in-process —
    # without this finally an n-run test would accumulate one open fd
    # per run.
    try:
        if cohort_mode == "arrow-streaming":
            chrom_iter = _iter_chrom_sites_streaming(
                chromosomes_to_simulate=chromosomes_to_simulate,
                args=args, demo_model=demo_model, rng=rng,
                resume=resume,
                clinvar_index=clinvar_index, rsid_pool=rsid_pool,
                cosmic_pool=cosmic_pool,
                overlay_stats=overlay_stats, sfs_total=sfs_total,
                chunk_size_mb=chunk_size_mb,
                fasta=fasta,
            )
        else:
            chrom_iter = _iter_chrom_sites_materialised(
                chromosomes_to_simulate=chromosomes_to_simulate,
                args=args, demo_model=demo_model, rng=rng,
                workers=workers,
                resume=resume,
                clinvar_index=clinvar_index, rsid_pool=rsid_pool,
                cosmic_pool=cosmic_pool,
                overlay_stats=overlay_stats, sfs_total=sfs_total,
                chunk_size_mb=chunk_size_mb,
                fasta=fasta,
            )

        for chrom, sites in chrom_iter:
            # Per-chrom BCF — `cohort.chr<N>.bcf`. The relative path
            # goes straight into the manifest's cohort_bcfs list.
            chrom_bcf = cohort_dir / f"cohort.chr{chrom}.bcf"
            _write_cohort_chrom_bcf(
                mode=cohort_mode,
                chrom=chrom,
                chrom_bcf=chrom_bcf,
                sites=sites,
                sample_ids=sample_ids,
                build=args.build,
                workers=workers,
                cohort_dir=cohort_dir,
                arrow_batch_size=args.cohort_arrow_batch_size,
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
                eta = (
                    (remaining / rate_chr)
                    if rate_chr > 0 else float("inf")
                )
                print(
                    f"  cohort BCFs: {done}/{len(chromosomes)} "
                    f"chromosomes written (elapsed "
                    f"{_format_duration(elapsed)}, "
                    f"eta {_format_duration(eta)})",
                    file=sys.stderr,
                )
                last_progress_log = now
    finally:
        if fasta is not None:
            fasta.close()

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
            # M13.1: per-person sex, parallel-indexed to ``samples``.
            "sex": resume.sexes,
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
        # M13.3: per-person sex for haploid emission. ``resume.sexes``
        # is the canonical source (deterministic per-seed, persisted
        # across resumes) on the streamed coalescent path.
        "person_sexes": resume.sexes,
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
        # M13.1: per-person sex, parallel-indexed to ``samples``.
        "sex": resume.sexes,
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
    parser = _parser(script_dir)
    args = parser.parse_args(argv)

    # ``--print-config`` short-circuits before any config discovery
    # or simulation work: emit the starter YAML to stdout and exit
    # so the user can redirect it to a file as a starting point.
    if getattr(args, "print_config", False):
        from .config import render_default_config_yaml
        sys.stdout.write(render_default_config_yaml())
        return 0

    # ------------------------------------------------------------
    # Optional YAML config layer (Phase config-file). CLI flags
    # already win; config layer fills in anything not on the CLI;
    # defaults apply where neither did. Discovery is cwd-only.
    # ------------------------------------------------------------
    from .config import (
        DEFAULT_CONFIG_FILENAME,
        discover_config_file,
        load_and_validate_config,
        merge_config_into_args,
        parse_explicit_cli_args,
        format_effective_values,
    )

    config_path = None
    if getattr(args, "config", None):
        config_path = Path(args.config)
        if not config_path.is_file():
            sys.exit(f"--config: file not found: {config_path}")
    elif not getattr(args, "no_config", False):
        config_path = discover_config_file(Path.cwd())

    config_obj = None
    if config_path is not None:
        print(
            f"  Loading values from config file: {config_path}",
            file=sys.stderr,
        )
        config_obj = load_and_validate_config(config_path)
        explicit_cli = parse_explicit_cli_args(parser, argv)
        # Re-parse so the parser's defaults are restored if
        # parse_explicit_cli_args left anything inconsistent.
        args = merge_config_into_args(
            args, config_obj, explicit_cli, parser=parser,
        )
        parser_defaults = {
            a.dest: a.default for a in parser._actions
            if a.dest not in ("help", "config", "no_config")
        }
        effective = format_effective_values(
            args, parser_defaults, config_obj, explicit_cli,
        )
        if effective:
            print("  Effective non-default values:", file=sys.stderr)
            for line in effective:
                print(line, file=sys.stderr)

    if args.check_deps:
        return _check_deps()

    # Hard-fail on missing htslib binaries even without --check-deps.
    for tool in ("bcftools", "tabix", "bgzip"):
        if not shutil.which(tool):
            sys.exit(f"required tool not on PATH: {tool}")

    if args.art:
        sys.exit(
            "--art (ART read simulation + bcftools call) is gated "
            "with a clear rejection: M12 (2026-05-14) ships the "
            "reference FASTA cache so the dep is in place, but the "
            "ART pipeline itself (art_illumina + bcftools call "
            "wiring) hasn't been hooked up yet. Use the default "
            "lightweight noise model (--error-rate / --dropout-rate) "
            "for now."
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
        # M12: resolve the reference FASTA (handles
        # ``--no-reference-fasta``, explicit ``--reference-fasta``,
        # AND auto-fetch from ``cache_dir/reference/<build>.fa``).
        # The open handle is discarded here — the admixture path
        # uses a ``ProcessPoolExecutor`` so workers re-open from
        # the path inside ``simulate_chromosome``; the parent-side
        # ``load + validate`` still runs to fail fast at startup
        # on bad FASTA / missing chrom (PR #84 review).
        _fa_admix, _fa_path = _resolve_reference_fasta(args, chromosomes)
        try:
            cohort_sites, person_ancestry = simulate_admixed_cohort(
                chromosomes=chromosomes, build=args.build,
                n_people=args.n, length_mb=args.chr_length_mb,
                proportions=proportions,
                rec_rate=args.rec_rate, mu=args.mu,
                rng=rng, verbose=True, workers=workers,
                fasta_path=(str(_fa_path) if _fa_path else None),
            )
        finally:
            if _fa_admix is not None:
                _fa_admix.close()
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
    # M13.1: per-person sex assignment from a SEPARATE rng so the
    # master rng state is unchanged. Drawing sex from ``rng`` here
    # would shift every downstream consumer (overlay seeds, error
    # model, etc.) and a fixed-seed run would no longer reproduce
    # pre-M13.1 output. ``_draw_sexes`` is deterministic given
    # ``(args.seed, args.n, args.male_fraction)``.
    from .resume import _draw_sexes
    person_sexes = _draw_sexes(args.seed, args.n, args.male_fraction)

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
            # M13.1: per-person sex, parallel-indexed to ``samples``.
            "sex": person_sexes,
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
        # M13.3: per-person sex from ``person_sexes`` drawn earlier
        # in main() — same _draw_sexes helper that resume.sexes
        # uses on the streamed path, so the two code paths emit
        # consistent sex assignments at the same --seed.
        "person_sexes": person_sexes,
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
        # M13.1: per-person sex assignment, parallel-indexed to
        # ``samples`` (i.e. ``sex[i]`` is the sex of ``samples[i]``).
        # Top-level rather than per-person so the field is available
        # in cohort-only mode where ``people`` is empty. Each entry
        # is ``"m"`` or ``"f"``. The simulation itself doesn't yet
        # use this (M13.3+ wires ploidy / PAR / MT clonality through)
        # — but it's recorded now so consumers can already reason
        # about the cohort's sex composition.
        "sex": person_sexes,
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
