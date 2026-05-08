"""Coalescent backbone — msprime + stdpopsim driver for LD-correct cohorts.

Replaces the M4 1000G-pool + power-law SFS path with a proper coalescent
simulation. `stdpopsim` supplies a standard human demographic model
(default: `OutOfAfrica_3G09`, sampling from CEU) and chromosome lengths
from the GRCh38 assembly; `msprime` does the ancestry + mutation draw.
Mutations are binary on the tree, so we synthesise REF/ALT base pairs
from the M3+ Ti/Tv calibrator (`titv.choose_alt`) to hit the human ~2.1
ratio on the de-novo SNVs.

The output is a list of sites in exactly the same shape M4's cohort
sampler returns (`chrom`, `pos`, `id`, `ref`, `alts`, `afs`, `acs`,
`gts`), so the CLI's per-person writer loop is unchanged.
"""

from __future__ import annotations

import random
import sys

from .builds import BUILDS
from .titv import DEFAULT_TARGET_TITV, choose_alt


DEFAULT_DEMO_MODEL = "OutOfAfrica_3G09"
DEFAULT_POPULATION = "CEU"
DEFAULT_REC_RATE = 1e-8        # uniform recombination rate (per bp per gen)
DEFAULT_MU = 1.29e-8           # human average mutation rate
DEFAULT_CHR_LENGTH_MB = 5.0    # simulated span per chromosome (0 = full)


# ---------------------------------------------------------------------------
# Phase 5f — chunked simulation memory model
# ---------------------------------------------------------------------------
# Auto-pick a per-chunk simulation length such that msprime's working
# memory fits the host. Calibrated against the user-provided
# ``--profile-memory`` trace from a 16 GB host running the failing
# config (``--n 3000 --chromosomes 1-22 --chr-length-mb 70``,
# OutOfAfrica_3G09, 4 auto-workers). The trace showed children RSS
# climbing to ~16 GB total (4 workers × ~4 GB each) at the
# auto-picked chunk size of 8.7 Mb. That gives ~153 KiB of working
# memory per (sample × Mb) at OOA scale — twice my original
# 80 KiB calibration which was extrapolated from a single
# full-chromosome OOM observation rather than from a ratio'd
# measurement at a known chunk size. The previous coefficient
# under-picked, the host hit its RAM ceiling, and workers
# swap-thrashed for 46 minutes without progressing. Round 153 up to
# 160 KiB for a safety margin against demographic models heavier
# than OOA_3G09.
#
# Constant-Ne msprime (`--demo-model none`) is still roughly 5×
# cheaper because there are no bottleneck spikes in the active-
# lineage count; the constant scales proportionally.
#
# These constants are intentionally pessimistic — better to pick a
# slightly smaller chunk than the host strictly needs and waste a bit
# of throughput than to OOM mid-simulation. Users with surprising
# working memory can override via `--chr-chunk-mb N`.
CHUNK_RAM_BYTES_PER_SAMPLE_PER_MB_OOA = 160 * 1024
CHUNK_RAM_BYTES_PER_SAMPLE_PER_MB_CONSTANT_NE = 32 * 1024

# Default safety target for the auto-pick. The previous 50% was the
# total budget across all workers and the parent process; the trace
# showed that at 50% the host's RAM ceiling was being hit when model
# error stacked with parent overhead + ClinVar pool + bcftools
# subprocesses. 25% leaves headroom for those plus any residual
# model error.
DEFAULT_AUTO_PICK_TARGET_FRACTION = 0.25

# When a user requests ``--workers W`` (or it auto-resolves to
# cpu_count), the auto-pick may have to shrink chunks below a useful
# size to fit the budget. Below this floor, msprime's per-chunk
# startup cost starts to dominate the per-chunk simulation time, and
# the boundary-smoothing benefit erodes. If chunk size would drop
# below this threshold, we'd rather drop a worker than pick smaller
# chunks. The auto-derate path checks this.
CHUNK_AUTO_DERATE_FLOOR_MB = 2.0

# Each chunk simulates a small overlap margin past its declared end so
# the central region's coalescent context isn't truncated abruptly at
# the boundary. Variants in the trailing overlap are dropped at write
# time so the per-chrom BCF stays duplicate-free. Documented as
# *boundary smoothing*, not true cross-chunk LD recovery — chunks are
# still independent simulations.
CHUNK_OVERLAP_FRACTION = 0.10   # 10% of chunk size
CHUNK_OVERLAP_MIN_BP = 500_000  # 0.5 Mb floor — too small and the
                                # overlap doesn't help the boundary
CHUNK_OVERLAP_MAX_BP = 5_000_000  # 5 Mb ceiling — past that we'd be
                                  # paying for simulation we won't use


def estimate_chunk_ram_bytes(n_people: int, chunk_size_mb: float,
                             demo_model: str | None) -> int:
    """Estimate msprime's per-chunk peak working memory in bytes.

    Used by the auto-pick path to pick the largest chunk size whose
    estimated working set fits the host. Pessimistic by design — the
    cost of a slight under-pick is wasted throughput; the cost of an
    over-pick is OOM mid-run.
    """
    if demo_model is None or (
        isinstance(demo_model, str) and demo_model.lower() == "none"
    ):
        rate = CHUNK_RAM_BYTES_PER_SAMPLE_PER_MB_CONSTANT_NE
    else:
        rate = CHUNK_RAM_BYTES_PER_SAMPLE_PER_MB_OOA
    return int(rate * n_people * chunk_size_mb)


def auto_pick_chunk_size_mb(
    n_people: int, length_mb: float,
    demo_model: str | None,
    available_bytes: int,
    workers: int = 1,
    target_fraction: float = DEFAULT_AUTO_PICK_TARGET_FRACTION,
) -> float:
    """Pick the largest chunk size whose estimated working set fits
    in ``target_fraction × available_bytes / workers``.

    Workers are accounted for because each parallel worker holds one
    chunk's tree sequence simultaneously — peak system RAM is
    ``workers × per-chunk RAM``.

    Returns the chunk size in megabases. If a single full-chromosome
    simulation already fits the budget, returns ``length_mb`` (no
    chunking needed).
    """
    if length_mb <= 0:
        return 0.0
    per_worker_target = max(
        1, int(available_bytes * target_fraction // max(1, workers))
    )
    full_estimate = estimate_chunk_ram_bytes(n_people, length_mb, demo_model)
    if full_estimate <= per_worker_target:
        return length_mb
    # Scale the chunk size linearly with the budget overshoot. The
    # estimate is linear in chunk_size_mb (rate × n × mb), so the
    # right chunk size is just length_mb × (target / full_estimate).
    factor = per_worker_target / full_estimate
    chunk_mb = max(1.0, length_mb * factor)
    return chunk_mb


def auto_derate_workers(
    n_people: int, length_mb: float,
    demo_model: str | None,
    available_bytes: int,
    requested_workers: int,
    target_fraction: float = DEFAULT_AUTO_PICK_TARGET_FRACTION,
    floor_chunk_mb: float = CHUNK_AUTO_DERATE_FLOOR_MB,
) -> int:
    """Possibly downgrade ``requested_workers`` so the auto-picked
    chunk size stays at or above ``floor_chunk_mb``.

    The user-reported failure mode was 4 parallel workers each
    holding a ~4 GB tree sequence at the auto-picked chunk size,
    saturating the host's 16 GB RAM and stalling in swap. Symmetric
    fix: when a chunk size drops below the floor, reduce workers
    instead of accepting tiny chunks. Tiny chunks pay msprime's
    per-chunk startup cost too often and erode the 5f boundary-
    smoothing benefit; reducing workers preserves chunk size at the
    cost of less parallelism (which is the right trade-off when
    RAM is the bound, not CPU).

    Caps at the requested worker count and a floor of 1 — this
    function never *increases* parallelism. Callers that pass
    ``requested_workers=1`` get 1 back unconditionally.
    """
    if requested_workers <= 1:
        return 1
    # Walk down from the requested count; pick the largest worker
    # count where the per-worker chunk size auto-picks at or above
    # the floor. The walk is short (≤ cpu_count) so a linear scan
    # is fine.
    for w in range(requested_workers, 0, -1):
        chunk_mb = auto_pick_chunk_size_mb(
            n_people, length_mb, demo_model, available_bytes, w,
            target_fraction,
        )
        if chunk_mb >= floor_chunk_mb:
            return w
    # Even at workers=1 the chunk would be sub-floor. Honour the
    # request anyway — at workers=1 a smaller chunk is the only way
    # to make progress without a host RAM upgrade. The CLI prints a
    # warning at flag-resolution time so the user knows what hit
    # them.
    return 1


def _require_deps():
    """Import msprime/stdpopsim lazily so the package works without them."""
    try:
        import msprime  # noqa: F401
        import stdpopsim  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "Coalescent path requires msprime + stdpopsim. Install via "
            "`pip install -r synthetic_people/requirements.txt` or pass "
            "--legacy-background to use the M4 1000G-pool path."
        ) from exc


def simulate_chromosome(chrom: str, build: str, n_people: int,
                        length_mb: float, demo_model: str | None,
                        population: str, rec_rate: float, mu: float,
                        rng: random.Random,
                        titv_target: float = DEFAULT_TARGET_TITV,
                        chunk_size_mb: float = 0.0) -> list:
    """Simulate one chromosome and return a list of cohort site dicts.

    With `demo_model` set, we drive the simulation through stdpopsim so
    the model's own demographic history and (optionally) recombination
    map are applied. With `demo_model=None` we fall back to a constant-
    size Ne=10_000 single-population coalescent driven directly by
    `msprime.sim_ancestry` — faster and dependency-lighter, but no LD
    map and no realistic demography.

    Phase 5f: when ``chunk_size_mb > 0`` and smaller than the simulated
    chromosome length, the chromosome is split into independent
    sub-chunks. Each chunk is its own msprime simulation against a
    contig of length ``chunk_size_mb + overlap_margin``; per-chunk
    seeds derive deterministically from a chromosome seed plus the
    chunk index, so re-running the same configuration produces
    byte-identical chunk outputs. Variants are emitted with positions
    offset by ``chunk_index × chunk_size_bp`` in cohort coordinates.
    The trailing overlap region of each chunk (positions ≥
    ``chunk_size_bp``) is dropped at write time so the per-chrom site
    list stays duplicate-free.

    Tradeoff: cross-chunk LD is lost — chunks are independent
    simulations, so haplotypes are uncorrelated across chunk
    boundaries. Documented as *boundary smoothing* in
    PERFORMANCE_PLAN §"Phase 5f".
    """
    _require_deps()

    if chrom not in BUILDS[build]["contigs"]:
        raise ValueError(f"unknown chromosome {chrom!r} for build {build}")
    chrom_length = BUILDS[build]["contigs"][chrom]
    sim_length_bp = int(chrom_length if length_mb <= 0
                        else min(chrom_length, length_mb * 1_000_000))

    chunk_size_bp = int(chunk_size_mb * 1_000_000) if chunk_size_mb > 0 else 0
    if chunk_size_bp <= 0 or chunk_size_bp >= sim_length_bp:
        # Single-chunk path: today's behaviour. msprime simulates the
        # full requested length in one pass.
        seed = rng.randint(1, 2**31 - 1)
        ts = _simulate_one(chrom, build, n_people, sim_length_bp,
                           demo_model, population, rec_rate, mu, seed)
        return _tree_sequence_to_sites(
            ts, chrom, n_people, rng, titv_target, position_offset=0)

    return _simulate_chromosome_chunked(
        chrom, build, n_people, sim_length_bp, chunk_size_bp,
        demo_model, population, rec_rate, mu, rng, titv_target,
    )


def _simulate_one(chrom: str, build: str, n_people: int,
                  sim_length_bp: int, demo_model: str | None,
                  population: str, rec_rate: float, mu: float,
                  seed: int):
    """Run msprime once for a single contig of length ``sim_length_bp``.

    Pulled out so both the whole-chromosome path and the chunked path
    share the same simulation invocation. Returns the raw tree
    sequence; the caller is responsible for iterating variants and
    freeing the tree sequence when done.
    """
    import msprime
    import stdpopsim

    if demo_model is not None:
        species = stdpopsim.get_species("HomSap")
        # `right=sim_length_bp` slices the contig to a prefix region.
        # Each chunk simulates *as if* it were the first
        # sim_length_bp bp of `chrom`, with chunk-specific seeds
        # making chunks independent.
        contig = species.get_contig(chrom, right=sim_length_bp)
        model = species.get_demographic_model(demo_model)
        pop_names = [p.name for p in model.populations]
        if population not in pop_names:
            raise ValueError(
                f"population {population!r} not in {demo_model}: "
                f"{pop_names}"
            )
        engine = stdpopsim.get_engine("msprime")
        return engine.simulate(model, contig, {population: n_people},
                               seed=seed)
    ts = msprime.sim_ancestry(
        samples=n_people,
        population_size=10_000,
        sequence_length=sim_length_bp,
        recombination_rate=rec_rate,
        random_seed=seed,
    )
    return msprime.sim_mutations(
        ts, rate=mu, random_seed=seed,
        model=msprime.BinaryMutationModel(),
    )


def _chunk_overlap_bp(chunk_size_bp: int) -> int:
    """Per-chunk overlap margin in bp — clamped to the configured min/max."""
    overlap = int(chunk_size_bp * CHUNK_OVERLAP_FRACTION)
    return max(CHUNK_OVERLAP_MIN_BP, min(CHUNK_OVERLAP_MAX_BP, overlap))


def _simulate_chromosome_chunked(
    chrom: str, build: str, n_people: int,
    sim_length_bp: int, chunk_size_bp: int,
    demo_model: str | None, population: str,
    rec_rate: float, mu: float,
    rng: random.Random, titv_target: float,
) -> list:
    """Drive the per-chrom simulation as a sequence of independent
    sub-chunks. Each chunk's tree sequence is built, iterated, and
    freed before the next chunk's simulation starts.

    Determinism: a single ``chrom_seed`` is drawn from ``rng`` up
    front, and each chunk's per-msprime seed is derived from
    ``chrom_seed + chunk_index`` via a Knuth-style multiplicative mix.
    This way the same ``--seed`` produces the same chunks regardless
    of how many were already simulated in a previous resume.
    """
    overlap_bp = _chunk_overlap_bp(chunk_size_bp)
    chrom_seed = rng.randint(1, 2**31 - 1)

    sites: list = []
    chunk_index = 0
    cursor_bp = 0
    while cursor_bp < sim_length_bp:
        chunk_end_bp = min(cursor_bp + chunk_size_bp, sim_length_bp)
        is_last = (chunk_end_bp >= sim_length_bp)
        chunk_extent_bp = (chunk_end_bp - cursor_bp
                           + (overlap_bp if not is_last else 0))

        # Knuth-style multiplicative hash to mix chrom_seed and
        # chunk_index into a uniform 31-bit chunk seed. Avoids the
        # rng-state-dependence of `rng.randint` per chunk so a resumed
        # run sees the same chunk seed regardless of which chunks
        # have already been simulated.
        chunk_seed = (chrom_seed + chunk_index * 0x9E3779B9) & 0x7FFFFFFF
        if chunk_seed == 0:
            chunk_seed = 1  # msprime rejects seed=0
        chunk_rng = random.Random(chunk_seed)

        ts = _simulate_one(
            chrom, build, n_people, chunk_extent_bp,
            demo_model, population, rec_rate, mu, chunk_seed,
        )
        # ``chunk_size_bp`` is the chunk's logical size in chunk-local
        # coordinates — variants past it are in the trailing-overlap
        # region and belong to the next chunk's simulation, so we
        # drop them to keep the per-chrom site list dedup-free.
        keep_below = chunk_end_bp - cursor_bp
        chunk_sites = _tree_sequence_to_sites(
            ts, chrom, n_people, chunk_rng, titv_target,
            position_offset=cursor_bp,
            keep_below_bp=keep_below,
        )
        # Drop the tree sequence reference so its msprime working
        # memory can be reclaimed before the next chunk's sim.
        del ts
        sites.extend(chunk_sites)
        cursor_bp = chunk_end_bp
        chunk_index += 1
    return sites


def _tree_sequence_to_sites(ts, chrom: str, n_people: int,
                            rng: random.Random,
                            titv_target: float,
                            position_offset: int = 0,
                            keep_below_bp: int | None = None) -> list:
    """Convert an msprime TreeSequence to cohort site dicts.

    Phase 5c: emits sparse ``carriers`` rather than dense
    ``gts: list[str]``. Per-site memory drops from
    ``O(n_people)`` GT strings to ``O(non-zero entries)`` int
    tuples — the difference between hundreds of MB and a few MB
    per site at ``n=100 000``. See ``cohort_sites.py`` for the
    helper functions that round-trip between the two shapes.

    Phase 5f chunked-simulation hooks:

    - ``position_offset`` shifts every emitted variant's POS by the
      given bp count, so a chunk simulating "the first chunk_extent
      of chr1" can land its variants at the correct cohort
      coordinates without re-simulating the chromosome prefix.
    - ``keep_below_bp`` (when set) drops variants whose simulated
      position is at or above this bp threshold. Used by the chunked
      path to discard variants in the trailing-overlap region of
      each chunk (those positions belong to the next chunk's
      simulation, kept dedup-free by emitting them only once from
      whichever chunk owns them).
    """
    n_haplotypes = 2 * n_people
    sites = []
    used_positions: set = set()
    for var in ts.variants():
        # Drop trailing-overlap variants on the chunked path. We test
        # against the *simulated* position (pre-offset) so the
        # comparison is straightforward: if the simulated coordinate
        # is ≥ keep_below_bp, the variant belongs to the next chunk.
        if (keep_below_bp is not None
                and int(var.site.position) >= keep_below_bp):
            continue

        # msprime positions are float along [0, sequence_length). VCF
        # positions are 1-based integers. Advance past collisions so no
        # two sites land on the same POS. The chunk position_offset
        # shifts emitted POS into cohort coordinates.
        pos = int(var.site.position) + 1 + position_offset
        while pos in used_positions:
            pos += 1
        used_positions.add(pos)

        # var.genotypes: array of allele indices per haplotype slot.
        # For BinaryMutationModel this is {0, 1}; for stdpopsim's
        # default JC69MutationModel, recurrent mutations at the same
        # position can produce alleles 2+ (multi-allelic). The cohort
        # site dict declares a single ALT base, so emitting allele
        # index 2+ references an ALT that doesn't exist — bcftools
        # rejects it downstream with `Incorrect allele ("2") in
        # <SAMPLE>`. Multi-allelic would need per-alt accounting (one
        # ALT base per index, per-alt AC/AF, multi-alt VCF output);
        # for now we keep cohort sites strictly biallelic and skip
        # the rest.
        gts_arr = var.genotypes
        if int(gts_arr.max(initial=0)) > 1:
            continue
        n_alt_haplotypes = int((gts_arr > 0).sum())
        if n_alt_haplotypes == 0 or n_alt_haplotypes == n_haplotypes:
            # Fixed sites (no variation across cohort) are not variants.
            # Shouldn't happen with binary mutations in practice, but we
            # defend anyway so the output stays well-formed.
            continue

        ref = rng.choice(("A", "C", "G", "T"))
        alt = choose_alt(ref, rng, target=titv_target)
        assert alt is not None  # ref is always a standard base

        # Sparse carriers: only emit the non-zero entries. With binary
        # mutations every non-zero entry is allele index 1, so
        # ``carriers`` reduces to a list of haplotype indices paired
        # with the constant 1; we still store the tuple form for shape
        # parity with the multi-allelic legacy path.
        carriers = [
            (int(idx), int(allele))
            for idx, allele in enumerate(gts_arr) if allele > 0
        ]
        sites.append({
            "chrom": chrom,
            "pos": pos,
            "id": ".",
            "ref": ref,
            "alts": [alt],
            "afs": [n_alt_haplotypes / n_haplotypes],
            "acs": [n_alt_haplotypes],
            "n_haplotypes": n_haplotypes,
            "carriers": carriers,
        })
    return sites


def _simulate_chromosome_from_seed(chrom: str, build: str, n_people: int,
                                   length_mb: float,
                                   demo_model: str | None,
                                   population: str, rec_rate: float,
                                   mu: float, seed: int,
                                   titv_target: float,
                                   chunk_size_mb: float = 0.0) -> list:
    """Worker entry point — builds its own rng from the given seed.

    Module-level so it's picklable across the ProcessPoolExecutor task
    boundary. ``chunk_size_mb`` is forwarded to
    :func:`simulate_chromosome`; ``0`` (the default) means single-pass
    full-chromosome simulation, ``> 0`` means split into chunks.
    """
    rng = random.Random(seed)
    return simulate_chromosome(
        chrom, build, n_people, length_mb, demo_model, population,
        rec_rate, mu, rng, titv_target,
        chunk_size_mb=chunk_size_mb,
    )


def simulate_cohort_iter(chromosomes: list, build: str, n_people: int,
                         length_mb: float, demo_model: str | None,
                         population: str, rec_rate: float, mu: float,
                         rng: random.Random,
                         titv_target: float = DEFAULT_TARGET_TITV,
                         verbose: bool = False,
                         workers: int = 1,
                         chunk_size_mb: float = 0.0):
    """Yield ``(chrom, sites)`` pairs one chromosome at a time.

    The streaming counterpart of :func:`simulate_cohort` — used by the
    Phase 5b cohort-streaming path so the in-memory footprint at any
    moment is bounded by a single chromosome's sites rather than the
    whole cohort. The flat-list :func:`simulate_cohort` below stays as
    a thin wrapper for callers (admixture, tests) that need everything
    materialised.

    Phase 5e Phase A: simulation is now serial across chromosomes,
    regardless of ``workers``. The pre-5e parallel-chromosome path
    multiplied msprime tree-sequence RAM by the worker count and
    OOM-killed workstation-class hosts at ``n=3000+``; that's the
    failure mode 5e was scoped to fix. The cohort BCF write is what
    parallelises now (in the caller via
    :func:`bcf_writer.write_cohort_bcf_parallel`), not the simulation.
    The ``workers`` argument is kept on the signature for API
    compatibility but is no longer consulted by this generator.

    Determinism: per-chromosome seeds are still pre-derived from
    ``rng`` upfront, identical to :func:`simulate_cohort`, so the
    same ``--seed`` yields the same per-chrom site list across runs.
    Output is byte-identical for any ``workers`` value because the
    simulation path no longer branches on it.
    """
    del workers  # signature kept for API stability — see docstring
    seeds = [rng.randint(1, 2**31 - 1) for _ in chromosomes]

    for chrom, seed in zip(chromosomes, seeds):
        if verbose:
            print(f"  simulating chrom {chrom} (length {length_mb} "
                  f"Mb, model={demo_model or 'uniform'})...",
                  file=sys.stderr)
        sites = _simulate_chromosome_from_seed(
            chrom, build, n_people, length_mb, demo_model,
            population, rec_rate, mu, seed, titv_target,
            chunk_size_mb,
        )
        if verbose:
            print(f"    {len(sites)} variable sites on chrom "
                  f"{chrom}", file=sys.stderr)
        yield chrom, sites
        # Drop the generator-local binding so that, once the
        # consumer's own ``del sites`` runs, the per-chrom site
        # list (multi-GB at n≈3000 because the sparse-carriers
        # representation has heavy Python tuple overhead) is
        # actually collected before the next chromosome's
        # simulation starts. Otherwise the generator frame keeps
        # the previous chrom's list alive across the yield, and
        # peak working set is two chromosomes wide.
        sites = None


def simulate_cohort(chromosomes: list, build: str, n_people: int,
                    length_mb: float, demo_model: str | None,
                    population: str, rec_rate: float, mu: float,
                    rng: random.Random,
                    titv_target: float = DEFAULT_TARGET_TITV,
                    verbose: bool = False,
                    workers: int = 1,
                    chunk_size_mb: float = 0.0) -> list:
    """Simulate a cohort and return all sites as a single flat list.

    Thin wrapper over :func:`simulate_cohort_iter` for callers that
    need the full cohort materialised (admixture path, fixture builders
    in tests). Phase 5b's cohort-streaming flow uses
    :func:`simulate_cohort_iter` directly so peak RAM stays bounded by
    one chromosome's working set.
    """
    all_sites: list = []
    for _, sites in simulate_cohort_iter(
        chromosomes, build, n_people, length_mb, demo_model,
        population, rec_rate, mu, rng, titv_target, verbose, workers,
        chunk_size_mb,
    ):
        all_sites.extend(sites)
    return all_sites
