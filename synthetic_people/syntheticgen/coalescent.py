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


# ---------------------------------------------------------------------------
# Streaming-cohort primitives (PR 2 of the option-3 refactor)
# ---------------------------------------------------------------------------
#
# At WGS scale (chr1 = 249 Mb, ~1.5 M variants × n=3000) the in-memory
# sites list returned by ``_tree_sequence_to_sites`` reaches ~17 GB
# (memprof28). The streaming approach replaces "materialise everything,
# then mutate in place, then sort, then stream into Arrow" with:
#
#   pass 1: walk ts.variants() once, accept/reject under the same
#           filters, consume rng identically (ref + alt picks), build
#           a light ``sites_meta`` list of one tuple per accepted site.
#
#   plan:   the cli (PR 3) calls ``clinvar.plan_inject_clinvar`` /
#           ``dbsnp.plan_inject_rsids`` / ``cosmic.plan_inject_cosmic``
#           against ``sites_meta`` to get ``{index: overlay_record}``
#           plans. ``annotate_clinvar``'s match set is similarly
#           pre-computed from ``sites_meta`` (no carriers needed for
#           a (chrom, pos, ref, alt) key match).
#
#   pass 2: walk ts.variants() again, this time building full site
#           dicts (with carriers). The rng is restored to its
#           pre-pass-1 state via getstate/setstate so the same
#           per-variant rng.choice + choose_alt sequence happens —
#           guaranteeing byte-identical ref/alt picks. Pass 2 yields
#           sites in pos order: tree-derived sites are emitted at
#           tree pos (which is monotone); overlay-injected sites are
#           held in a small heap and drained when the tree-walk pos
#           crosses the overlay's pos.
#
# Memory bound: ~75 MB sites_meta + ~few-hundred MB peak heap at
# density=0.04 × 1.5 M sites = ~60 k injected records × ~few KB
# carriers each. Two orders of magnitude smaller than the 17 GB
# materialised path while preserving byte-identical output at every
# fixed seed (master rng + overlay rng).


def _tree_sequence_to_sites_meta(ts, chrom: str, n_people: int,
                                 rng: random.Random,
                                 titv_target: float,
                                 position_offset: int = 0,
                                 keep_below_bp: int | None = None
                                 ) -> list:
    """Pass 1 of the streaming-cohort walk: produce one
    ``(chrom, pos, ref, alt, n_alt_haplotypes, n_haplotypes)`` tuple
    per accepted variant *without* building carriers.

    Consumes rng in the same sequence as :func:`_tree_sequence_to_sites`
    — same filters, same ``rng.choice`` + ``choose_alt`` calls per
    accepted variant, same per-call position dedup advancement. The
    only difference is what's stored: a 6-tuple instead of a site
    dict with full ``carriers``.

    The tuple shape carries everything the overlay planners and the
    pass-2 walk need:

    - ``(chrom, pos)`` for the ``plan_inject_*`` index-picking
      planners and the ``annotate_clinvar`` match-key lookup.
    - ``(ref, alt)`` so pass 2 can avoid re-consuming rng for these
      picks (and so ``annotate_clinvar`` keys match the legacy
      ``(chrom, pos, ref, alt)`` form).
    - ``(n_alt_haplotypes, n_haplotypes)`` so AC/AN/AF can be re-
      derived in pass 2 without re-scanning ``var.genotypes``.

    At WGS chr1 × n=3000: ~1.5 M tuples × ~80 bytes ≈ 120 MB.
    """
    n_haplotypes = 2 * n_people
    meta: list = []
    used_positions: set = set()
    for var in ts.variants():
        if (keep_below_bp is not None
                and int(var.site.position) >= keep_below_bp):
            continue

        pos = int(var.site.position) + 1 + position_offset
        while pos in used_positions:
            pos += 1
        used_positions.add(pos)

        gts_arr = var.genotypes
        if int(gts_arr.max(initial=0)) > 1:
            continue
        n_alt_haplotypes = int((gts_arr > 0).sum())
        if n_alt_haplotypes == 0 or n_alt_haplotypes == n_haplotypes:
            continue

        ref = rng.choice(("A", "C", "G", "T"))
        alt = choose_alt(ref, rng, target=titv_target)
        assert alt is not None

        meta.append((chrom, pos, ref, alt,
                     n_alt_haplotypes, n_haplotypes))
    return meta


def _build_carriers_from_variant(var):
    """Extract the sparse ``[(hap_idx, allele_idx), ...]`` carriers
    list from a tskit :class:`Variant`. Same logic as the inline list
    comprehension in :func:`_tree_sequence_to_sites`, lifted to a
    helper so pass 2 of the streaming walk can call it without
    duplicating the per-haplotype loop.
    """
    gts_arr = var.genotypes
    return [(int(idx), int(allele))
            for idx, allele in enumerate(gts_arr) if allele > 0]


def _stream_cohort_pass2(ts, chrom: str, n_people: int,
                        rng: random.Random,
                        titv_target: float,
                        sites_meta: list,
                        inject_map: dict,
                        annotation_map: dict,
                        position_offset: int = 0,
                        keep_below_bp: int | None = None):
    """Pass 2 of the streaming-cohort walk: yield full site dicts in
    pos-sorted order, applying pre-computed overlays.

    Walks ``ts.variants()`` a second time, consumes rng identically to
    pass 1 (so per-variant ``rng.choice`` + ``choose_alt`` produce the
    same ref/alt picks the planners already decided against — caller is
    responsible for ``getstate``/``setstate`` around the two passes).
    For each accepted site:

    - If the post-filter index is in ``inject_map``: builds the
      overlay site dict (overlay's pos / ref / alts / id / overlay-
      specific INFO fields, paired with the tree variant's carriers
      and the meta's AC/AN/AF), and pushes it into a min-heap keyed
      on overlay pos.
    - Otherwise: builds the tree-derived site dict (with annotation
      applied if applicable) at the tree pos, and pushes that into
      the same heap.

    After each iteration we compute a **safe-yield threshold** =
    smallest pos any future site can have. Future tree sites have
    pos > current_pos (monotone walk + dedup); future overlay sites
    have pos >= ``overlay_positions_sorted[overlay_next_idx]``. The
    threshold is ``min(pos + 1, next-overlay-pos)``. Anything in the
    heap with ``pos < threshold`` is safe to yield.

    Result: a stream of site dicts in monotone-by-pos order,
    structurally identical to what the materialised path produces
    post-``sites.sort(key=lambda s: s["pos"])``.

    Memory model: buffer holds tree-derived sites that are "ahead of"
    an un-emitted overlay (e.g. tree at pos=120k waiting for overlay
    at pos=75k) plus overlay sites awaiting their turn. With overlay
    positions roughly uniform across the chromosome the buffer stays
    small (O(sqrt(N_inject))); with overlays pathologically clustered
    at one end it can grow to ~all tree sites pending the cluster,
    which is the bounded worst case the materialised path always
    paid. Practical density × position distributions keep it in the
    hundreds-of-MB range — far below the 17 GB materialised cost.
    """
    import heapq

    # Pre-compute the future-min overlay position correctly.
    #
    # Earlier draft used ``overlay_positions_sorted[overlay_next_idx]``,
    # but ``overlay_next_idx`` counts overlays in *tree-walk order* while
    # the array was sorted by *position* — so the look-up reported the
    # (k+1)-th smallest position, not "smallest pos among overlays not
    # yet encountered." That false bound let tree sites yield ahead of
    # a future overlay whose pos was smaller than the current tree pos
    # but larger than already-encountered overlay positions. Concretely:
    # if the first overlay we hit has pos 500 and a later overlay has
    # pos 100, the buggy code allowed tree sites at pos 200 to yield
    # before pos 100 was even discovered. PR #51 Copilot review caught
    # the bug; this PR fixes it.
    #
    # The fix: sort the overlays by *accepted_idx* (tree-walk order),
    # then build a suffix-min array over the resulting pos sequence.
    # At any step we know the smallest pos among un-encountered
    # overlays by reading ``suffix_min[overlay_next_ptr]``, where the
    # pointer advances each time we encounter an injected idx.
    overlay_by_idx_order = sorted(
        [(idx, rec["pos"]) for idx, rec in inject_map.items()],
        key=lambda kv: kv[0],
    )
    n_inj = len(overlay_by_idx_order)
    suffix_min_overlay_pos: list = [float("inf")] * (n_inj + 1)
    for i in range(n_inj - 1, -1, -1):
        suffix_min_overlay_pos[i] = min(
            overlay_by_idx_order[i][1], suffix_min_overlay_pos[i + 1],
        )
    overlay_next_ptr = 0

    buffer: list = []  # min-heap of (pos, counter, site_dict)
    counter = 0
    accepted_idx = -1
    used_positions: set = set()
    n_haplotypes_total_cohort = n_haplotypes_total(n_people)

    for var in ts.variants():
        if (keep_below_bp is not None
                and int(var.site.position) >= keep_below_bp):
            continue

        pos = int(var.site.position) + 1 + position_offset
        while pos in used_positions:
            pos += 1
        used_positions.add(pos)

        gts_arr = var.genotypes
        if int(gts_arr.max(initial=0)) > 1:
            continue
        n_alt_haplotypes = int((gts_arr > 0).sum())
        if (n_alt_haplotypes == 0
                or n_alt_haplotypes == n_haplotypes_total_cohort):
            continue

        ref = rng.choice(("A", "C", "G", "T"))
        alt = choose_alt(ref, rng, target=titv_target)
        assert alt is not None

        accepted_idx += 1
        _, _, meta_ref, meta_alt, ac, n_haplotypes = sites_meta[accepted_idx]
        assert ref == meta_ref and alt == meta_alt, (
            f"streaming pass 2 rng divergence at idx {accepted_idx}: "
            f"got ({ref!r},{alt!r}) expected ({meta_ref!r},{meta_alt!r})"
        )

        carriers = _build_carriers_from_variant(var)

        if accepted_idx in inject_map:
            overlay = inject_map[accepted_idx]
            site = {
                "chrom": chrom,
                "pos": overlay["pos"],
                "id": overlay.get("id", "."),
                "ref": overlay["ref"],
                "alts": list(overlay["alts"]),
                "afs": [ac / n_haplotypes],
                "acs": [ac],
                "n_haplotypes": n_haplotypes,
                "carriers": carriers,
            }
            for k in ("clnsig", "clndn", "cosmic_id", "cosmic_gene"):
                if k in overlay:
                    site[k] = overlay[k]
            heapq.heappush(buffer, (site["pos"], counter, site))
            counter += 1
            overlay_next_ptr += 1
        else:
            site = {
                "chrom": chrom,
                "pos": pos,
                "id": ".",
                "ref": ref,
                "alts": [alt],
                "afs": [ac / n_haplotypes],
                "acs": [ac],
                "n_haplotypes": n_haplotypes,
                "carriers": carriers,
            }
            if accepted_idx in annotation_map:
                ann = annotation_map[accepted_idx]
                site["clnsig"] = ann["clnsig"]
                site["clndn"] = ann["clndn"]
                if site.get("id") in (None, "", ".") and ann.get("id"):
                    site["id"] = ann["id"]
            heapq.heappush(buffer, (pos, counter, site))
            counter += 1

        # Safe-yield threshold: smallest pos any future site can have.
        # Future tree pos > current pos (monotone + dedup); future
        # overlay pos is given by the suffix-min over un-encountered
        # overlays (see overlay_by_idx_order construction above).
        smallest_future_pos = min(
            pos + 1, suffix_min_overlay_pos[overlay_next_ptr],
        )

        while buffer and buffer[0][0] < smallest_future_pos:
            yield heapq.heappop(buffer)[2]

    # Drain remaining buffer entries (all guaranteed to be in order
    # since no future site can land below their pos at end-of-walk).
    while buffer:
        yield heapq.heappop(buffer)[2]


def n_haplotypes_total(n_people: int) -> int:
    """Cohort haplotype count = 2 × n_people (diploid). Helper used by
    the streaming pass 2 so the monomorphic-site filter line stays
    readable."""
    return 2 * n_people


def _build_annotation_map(sites_meta: list,
                          clinvar_records: list | None) -> dict:
    """Pre-compute the ``annotate_clinvar`` index→record mapping from a
    light sites_meta list, without needing the full sites list.

    Matches on ``(chrom, pos, ref, alt)`` identically to
    :func:`syntheticgen.clinvar.annotate_clinvar`. Returns
    ``{accepted_idx: clinvar_record}`` for indices that match. The
    streaming pass 2 then attaches ``clnsig`` / ``clndn`` / ``id``
    fields to those sites as they pass through.
    """
    if not clinvar_records:
        return {}
    index: dict = {}
    for r in clinvar_records:
        index[(r["chrom"], r["pos"], r["ref"], r["alt"])] = r
    out: dict = {}
    for i, meta in enumerate(sites_meta):
        chrom, pos, ref, alt = meta[0], meta[1], meta[2], meta[3]
        rec = index.get((chrom, pos, ref, alt))
        if rec is not None:
            out[i] = rec
    return out


def _plan_cohort_overlays(
    sites_meta: list,
    annotation_map: dict,
    *,
    clinvar_records: list | None,
    clinvar_inject_density: float,
    rsid_pool: list | None,
    rsid_density: float,
    cosmic_records: list | None,
    cosmic_inject_density: float,
    overlay_rng: random.Random,
):
    """Run the chained-overlay planners (clinvar → rsids → cosmic)
    against ``sites_meta`` *while mirroring the materialised path's
    sort-between-injects semantics*. Returns ``(clinvar_plan,
    rsid_plan, cosmic_plan)`` all keyed by **tree-walk idx**
    (= ``accepted_idx``) so the streaming pass 2 can look injections
    up by the same index it advances.

    Why the sort-between-injects matters:

      The materialised path is::

          inject_clinvar(sites, ...)   # mutates pos in place, then sorts
          clinvar_reserved = {i: sites[i].clnsig}   # POST-SORT indices
          inject_rsids(sites, ..., reserve=clinvar_reserved)
          ...

      ``inject_rsids`` sees the *post-clinvar-sort* sites_meta. Its
      ``used_keys`` collision set contains the new clinvar positions
      (not the original tree positions at the injected indices); its
      ``candidate`` list shuffles different integers than a streaming
      pre-sort caller would see. PR #51's first draft ran all three
      planners against the original tree-walk-order ``sites_meta``,
      which diverged from the materialised path on chained overlays
      — Copilot caught it. This helper mirrors the sort-between-
      injects logic so byte-identical output is restored.

    Mechanism:

      A light "view" tracks each tree-walk-idx's CURRENT (chrom, pos)
      plus two boolean flags — ``has_clnsig`` and ``has_rs_id`` —
      that match exactly what the materialised cli's reserve sets
      key on (``s.get("clnsig")`` and ``s["id"].startswith("rs")``).
      After each plan we mutate the view's positions for injected
      indices, update flags, and stable-sort by (chrom, pos). The
      next planner sees the post-sort view; we map its output
      indices back to tree-walk idx for the streaming pass 2.
    """
    from .clinvar import plan_inject_clinvar
    from .dbsnp import plan_inject_rsids
    from .cosmic import plan_inject_cosmic

    # View entry: [tree_walk_idx, chrom, pos, has_clnsig, has_rs_id]
    view: list = []
    for i, m in enumerate(sites_meta):
        view.append([i, m[0], m[1],
                     i in annotation_map,  # clinvar match by annotate
                     False])

    # 1) Plan clinvar inject against the tree-walk-order view
    #    (sites_meta_chrpos at this step is the (chrom, pos) of the
    #    view, which is still tree-walk order since no sort has run).
    clinvar_plan_treewalk = plan_inject_clinvar(
        [(v[1], v[2]) for v in view],
        clinvar_records or [], clinvar_inject_density, overlay_rng,
    )
    # Apply mutations
    for tw_idx, rec in clinvar_plan_treewalk.items():
        v = view[tw_idx]
        v[2] = rec["pos"]
        v[3] = True  # has_clnsig
    # Sort by (chrom, pos). Stable.
    view.sort(key=lambda e: (e[1], e[2]))

    # 2) Reserve = indices (post-clinvar-sort) where has_clnsig is True.
    clinvar_reserved_post_sort = {i for i, v in enumerate(view) if v[3]}

    # 3) Plan rsids against the post-clinvar-sort view
    rsid_plan_post_sort = plan_inject_rsids(
        [(v[1], v[2]) for v in view],
        rsid_pool or [], rsid_density, overlay_rng,
        reserve_indices=clinvar_reserved_post_sort,
    )
    # Map back to tree-walk idx
    rsid_plan_treewalk: dict = {}
    for ps_idx, rec in rsid_plan_post_sort.items():
        tw_idx = view[ps_idx][0]
        rsid_plan_treewalk[tw_idx] = rec
        v = view[ps_idx]
        v[2] = rec["pos"]
        v[4] = True  # has_rs_id (inject_rsids sets id=rec["rsid"] which starts with "rs")
    view.sort(key=lambda e: (e[1], e[2]))

    # 4) Reserve = indices (post-rsid-sort) where has_clnsig OR has_rs_id.
    all_reserved_post_sort = {
        i for i, v in enumerate(view) if v[3] or v[4]
    }

    # 5) Plan cosmic against the post-rsid-sort view
    cosmic_plan_post_sort = plan_inject_cosmic(
        [(v[1], v[2]) for v in view],
        cosmic_records or [], cosmic_inject_density, overlay_rng,
        reserve_indices=all_reserved_post_sort,
    )
    cosmic_plan_treewalk: dict = {}
    for ps_idx, rec in cosmic_plan_post_sort.items():
        tw_idx = view[ps_idx][0]
        cosmic_plan_treewalk[tw_idx] = rec

    return clinvar_plan_treewalk, rsid_plan_treewalk, cosmic_plan_treewalk


def stream_cohort_sites(
    ts, chrom: str, n_people: int,
    rng: random.Random,
    titv_target: float = DEFAULT_TARGET_TITV,
    *,
    clinvar_records: list | None = None,
    clinvar_inject_density: float = 0.0,
    rsid_pool: list | None = None,
    rsid_density: float = 0.0,
    cosmic_records: list | None = None,
    cosmic_inject_density: float = 0.0,
    overlay_rng: random.Random | None = None,
    position_offset: int = 0,
    keep_below_bp: int | None = None,
):
    """Top-level streaming-cohort entry point.

    Replaces the materialised pattern::

        sites = _tree_sequence_to_sites(ts, ..., rng=rng)
        annotate_clinvar(sites, clinvar_records)
        inject_clinvar(sites, clinvar_records, density, overlay_rng)
        inject_rsids(sites, rsid_pool, density, overlay_rng, reserve_indices=...)
        inject_cosmic(sites, cosmic_records, density, overlay_rng, reserve_indices=...)
        sites.sort(key=lambda s: s["pos"])
        # ... downstream consumer iterates sites

    with::

        for site in stream_cohort_sites(ts, chrom, n_people, rng,
                                        clinvar_records=clinvar_records,
                                        clinvar_inject_density=density,
                                        ..., overlay_rng=overlay_rng):
            # ... downstream consumer

    Output is **byte-identical** at every fixed combination of seeds
    (``rng`` and ``overlay_rng``) to the materialised path. The rng
    consumes the same number of randoms (each variant consumes
    ``rng.choice + choose_alt``); the per-pass rng state is bracketed
    via ``getstate``/``setstate`` so two passes produce the same
    cumulative consumption as one.

    Memory bound: ``sites_meta`` (~120 MB at WGS chr1) + overlay heap
    (~few hundred MB peak) — two orders of magnitude smaller than the
    materialised path.

    Caller still owns the ``ts.dump`` lifecycle / overlay record
    sourcing / sort-on-completion logic — this function only replaces
    the in-memory site list with a streaming iterator.
    """
    state_before = rng.getstate()

    sites_meta = _tree_sequence_to_sites_meta(
        ts, chrom, n_people, rng, titv_target,
        position_offset=position_offset,
        keep_below_bp=keep_below_bp,
    )

    annotation_map = _build_annotation_map(sites_meta, clinvar_records)

    rng_for_plans = overlay_rng if overlay_rng is not None else rng

    clinvar_plan, rsid_plan, cosmic_plan = _plan_cohort_overlays(
        sites_meta, annotation_map,
        clinvar_records=clinvar_records,
        clinvar_inject_density=clinvar_inject_density,
        rsid_pool=rsid_pool, rsid_density=rsid_density,
        cosmic_records=cosmic_records,
        cosmic_inject_density=cosmic_inject_density,
        overlay_rng=rng_for_plans,
    )

    # Combine: each plan is keyed by tree-walk idx; the three sets of
    # keys are disjoint by construction (reserve_indices), so update
    # order doesn't matter for correctness.
    inject_map: dict = {}
    inject_map.update(clinvar_plan)
    inject_map.update(rsid_plan)
    inject_map.update(cosmic_plan)

    # Restore rng state for pass 2 — same per-variant rng calls reproduce
    # the same ref/alt picks the planners already locked in.
    rng.setstate(state_before)

    yield from _stream_cohort_pass2(
        ts, chrom, n_people, rng, titv_target,
        sites_meta, inject_map, annotation_map,
        position_offset=position_offset,
        keep_below_bp=keep_below_bp,
    )


def simulate_chromosome_ts(chrom: str, build: str, n_people: int,
                          length_mb: float, demo_model: str | None,
                          population: str, rec_rate: float, mu: float,
                          rng: random.Random,
                          chunk_size_mb: float = 0.0):
    """Sibling of :func:`simulate_chromosome` that returns the raw
    TreeSequence (instead of a materialised site list) so the caller
    can stream variants via :func:`stream_cohort_sites`.

    Consumes the master ``rng`` identically to :func:`simulate_chromosome`
    *up to the point where _tree_sequence_to_sites would start consuming
    rng for ref/alt picks*. That is: one ``rng.randint`` call for the
    msprime seed, then return. The caller is then expected to pass the
    same ``rng`` into :func:`stream_cohort_sites`, which will consume
    further (identically to the materialised path) for its own
    ref/alt picks. Net cumulative rng consumption matches the
    materialised path call-for-call.

    For PR 3 the chunked path (``chunk_size_mb > 0`` producing multiple
    chunks) is not yet supported in streaming mode; the caller should
    auto-pick the materialised path when chunking is requested.
    """
    _require_deps()
    if chrom not in BUILDS[build]["contigs"]:
        raise ValueError(f"unknown chromosome {chrom!r} for build {build}")
    chrom_length = BUILDS[build]["contigs"][chrom]
    sim_length_bp = int(chrom_length if length_mb <= 0
                        else min(chrom_length, length_mb * 1_000_000))
    chunk_size_bp = int(chunk_size_mb * 1_000_000) if chunk_size_mb > 0 else 0
    if chunk_size_bp > 0 and chunk_size_bp < sim_length_bp:
        raise NotImplementedError(
            "Streaming-cohort mode does not yet support "
            "--chr-chunk-mb > 0. Pass --chr-chunk-mb 0 or use "
            "--cohort-mode arrow / sites_list."
        )
    seed = rng.randint(1, 2**31 - 1)
    return _simulate_one(chrom, build, n_people, sim_length_bp,
                         demo_model, population, rec_rate, mu, seed)


def simulate_cohort_ts_iter(chromosomes: list, build: str, n_people: int,
                            length_mb: float, demo_model: str | None,
                            population: str, rec_rate: float, mu: float,
                            rng: random.Random,
                            verbose: bool = False,
                            chunk_size_mb: float = 0.0):
    """Yield ``(chrom, ts, walk_rng)`` triples one chromosome at a time.

    Streaming counterpart of :func:`simulate_cohort_iter`. Where
    ``simulate_cohort_iter`` materialises the per-chrom sites list,
    this generator yields the raw TreeSequence and a fresh rng
    instance the caller can pass to :func:`stream_cohort_sites`.

    Per-chrom rng derivation matches :func:`simulate_cohort_iter` so
    a given ``rng`` seed produces byte-identical TreeSequences (and
    therefore byte-identical cohort BCFs after streaming + writing)
    regardless of which iterator the cli routes through.
    """
    seeds = [rng.randint(1, 2**31 - 1) for _ in chromosomes]

    for chrom, seed in zip(chromosomes, seeds):
        if verbose:
            print(f"  simulating chrom {chrom} (length {length_mb} "
                  f"Mb, model={demo_model or 'uniform'})...",
                  file=sys.stderr)
        chrom_rng = random.Random(seed)
        ts = simulate_chromosome_ts(
            chrom, build, n_people, length_mb, demo_model,
            population, rec_rate, mu, chrom_rng,
            chunk_size_mb=chunk_size_mb,
        )
        if verbose:
            print(f"    chrom {chrom} TreeSequence ready "
                  f"({ts.num_sites} variant sites)", file=sys.stderr)
        yield chrom, ts, chrom_rng
        # Drop local refs so the per-chrom TreeSequence can be GC'd
        # before the next chromosome's simulation starts (mirrors the
        # cleanup in :func:`simulate_cohort_iter`).
        ts = None
        chrom_rng = None


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
