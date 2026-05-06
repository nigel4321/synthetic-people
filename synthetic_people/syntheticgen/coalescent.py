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

import multiprocessing as mp
import random
import sys
from concurrent.futures import ProcessPoolExecutor

from .builds import BUILDS
from .titv import DEFAULT_TARGET_TITV, choose_alt


DEFAULT_DEMO_MODEL = "OutOfAfrica_3G09"
DEFAULT_POPULATION = "CEU"
DEFAULT_REC_RATE = 1e-8        # uniform recombination rate (per bp per gen)
DEFAULT_MU = 1.29e-8           # human average mutation rate
DEFAULT_CHR_LENGTH_MB = 5.0    # simulated span per chromosome (0 = full)


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
                        titv_target: float = DEFAULT_TARGET_TITV) -> list:
    """Simulate one chromosome and return a list of cohort site dicts.

    With `demo_model` set, we drive the simulation through stdpopsim so
    the model's own demographic history and (optionally) recombination
    map are applied. With `demo_model=None` we fall back to a constant-
    size Ne=10_000 single-population coalescent driven directly by
    `msprime.sim_ancestry` — faster and dependency-lighter, but no LD
    map and no realistic demography.
    """
    _require_deps()
    import msprime
    import stdpopsim

    if chrom not in BUILDS[build]["contigs"]:
        raise ValueError(f"unknown chromosome {chrom!r} for build {build}")
    chrom_length = BUILDS[build]["contigs"][chrom]
    sim_length = int(chrom_length if length_mb <= 0
                     else min(chrom_length, length_mb * 1_000_000))

    # Seeds: use rng to generate a deterministic int seed for msprime so
    # outer-level --seed controls the coalescent path too.
    seed = rng.randint(1, 2**31 - 1)

    if demo_model is not None:
        species = stdpopsim.get_species("HomSap")
        # `right=sim_length` slices the chromosome to a prefix region.
        contig = species.get_contig(chrom, right=sim_length)
        model = species.get_demographic_model(demo_model)
        pop_names = [p.name for p in model.populations]
        if population not in pop_names:
            raise ValueError(
                f"population {population!r} not in {demo_model}: "
                f"{pop_names}"
            )
        engine = stdpopsim.get_engine("msprime")
        ts = engine.simulate(model, contig, {population: n_people},
                             seed=seed)
    else:
        ts = msprime.sim_ancestry(
            samples=n_people,
            population_size=10_000,
            sequence_length=sim_length,
            recombination_rate=rec_rate,
            random_seed=seed,
        )
        ts = msprime.sim_mutations(
            ts, rate=mu, random_seed=seed,
            model=msprime.BinaryMutationModel(),
        )

    return _tree_sequence_to_sites(ts, chrom, n_people, rng, titv_target)


def _tree_sequence_to_sites(ts, chrom: str, n_people: int,
                            rng: random.Random,
                            titv_target: float) -> list:
    """Convert an msprime TreeSequence to cohort site dicts."""
    n_haplotypes = 2 * n_people
    sites = []
    used_positions: set = set()
    for var in ts.variants():
        # msprime positions are float along [0, sequence_length). VCF
        # positions are 1-based integers. Advance past collisions so no
        # two sites land on the same POS.
        pos = int(var.site.position) + 1
        while pos in used_positions:
            pos += 1
        used_positions.add(pos)

        # var.genotypes: array of allele indices per haplotype slot.
        # For BinaryMutationModel this is {0, 1}; genotypes > 0 carry
        # the derived allele (ALT). Multi-allelic would need per-allele
        # handling; at M5 we stick to biallelic.
        gts_arr = var.genotypes
        n_alt_haplotypes = int((gts_arr > 0).sum())
        if n_alt_haplotypes == 0 or n_alt_haplotypes == n_haplotypes:
            # Fixed sites (no variation across cohort) are not variants.
            # Shouldn't happen with binary mutations in practice, but we
            # defend anyway so the output stays well-formed.
            continue

        ref = rng.choice(("A", "C", "G", "T"))
        alt = choose_alt(ref, rng, target=titv_target)
        assert alt is not None  # ref is always a standard base

        gts = [f"{int(gts_arr[2 * i])}|{int(gts_arr[2 * i + 1])}"
               for i in range(n_people)]
        sites.append({
            "chrom": chrom,
            "pos": pos,
            "id": ".",
            "ref": ref,
            "alts": [alt],
            "afs": [n_alt_haplotypes / n_haplotypes],
            "acs": [n_alt_haplotypes],
            "gts": gts,
        })
    return sites


def _simulate_chromosome_from_seed(chrom: str, build: str, n_people: int,
                                   length_mb: float,
                                   demo_model: str | None,
                                   population: str, rec_rate: float,
                                   mu: float, seed: int,
                                   titv_target: float) -> list:
    """Worker entry point — builds its own rng from the given seed.

    Module-level so it's picklable across the ProcessPoolExecutor task
    boundary.
    """
    rng = random.Random(seed)
    return simulate_chromosome(
        chrom, build, n_people, length_mb, demo_model, population,
        rec_rate, mu, rng, titv_target,
    )


def simulate_cohort_iter(chromosomes: list, build: str, n_people: int,
                         length_mb: float, demo_model: str | None,
                         population: str, rec_rate: float, mu: float,
                         rng: random.Random,
                         titv_target: float = DEFAULT_TARGET_TITV,
                         verbose: bool = False,
                         workers: int = 1):
    """Yield ``(chrom, sites)`` pairs one chromosome at a time.

    The streaming counterpart of :func:`simulate_cohort` — used by the
    Phase 5b cohort-streaming path so the in-memory footprint at any
    moment is bounded by a single chromosome's sites rather than the
    whole cohort. The flat-list :func:`simulate_cohort` below stays as
    a thin wrapper for callers (admixture, tests) that need everything
    materialised.

    Determinism: per-chromosome seeds are pre-derived from ``rng``
    *before* spawning workers, identical to :func:`simulate_cohort`,
    so the same ``--seed`` yields the same per-chrom site list whether
    the caller iterates lazily or collects upfront.

    With ``workers > 1`` and more than one chromosome, msprime runs
    each chromosome in its own fork-pool worker. Yields are emitted in
    submission order (not completion order) so a downstream BCF writer
    sees chromosomes in the requested order.
    """
    seeds = [rng.randint(1, 2**31 - 1) for _ in chromosomes]

    use_pool = workers > 1 and len(chromosomes) > 1
    if not use_pool:
        for chrom, seed in zip(chromosomes, seeds):
            if verbose:
                print(f"  simulating chrom {chrom} (length {length_mb} "
                      f"Mb, model={demo_model or 'uniform'})...",
                      file=sys.stderr)
            sites = _simulate_chromosome_from_seed(
                chrom, build, n_people, length_mb, demo_model,
                population, rec_rate, mu, seed, titv_target,
            )
            if verbose:
                print(f"    {len(sites)} variable sites on chrom "
                      f"{chrom}", file=sys.stderr)
            yield chrom, sites
        return

    if verbose:
        print(f"  simulating {len(chromosomes)} chromosomes "
              f"(length {length_mb} Mb, "
              f"model={demo_model or 'uniform'}) "
              f"in parallel with {workers} workers...",
              file=sys.stderr)

    # fork shares the parent's already-loaded msprime/stdpopsim modules
    # via copy-on-write — much faster startup than spawn, which would
    # re-import them per worker.
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as ex:
        futures = [
            (chrom, ex.submit(
                _simulate_chromosome_from_seed,
                chrom, build, n_people, length_mb, demo_model,
                population, rec_rate, mu, seed, titv_target,
            ))
            for chrom, seed in zip(chromosomes, seeds)
        ]
        for chrom, fut in futures:
            sites = fut.result()
            if verbose:
                print(f"    {len(sites)} variable sites on chrom "
                      f"{chrom}", file=sys.stderr)
            yield chrom, sites


def simulate_cohort(chromosomes: list, build: str, n_people: int,
                    length_mb: float, demo_model: str | None,
                    population: str, rec_rate: float, mu: float,
                    rng: random.Random,
                    titv_target: float = DEFAULT_TARGET_TITV,
                    verbose: bool = False,
                    workers: int = 1) -> list:
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
    ):
        all_sites.extend(sites)
    return all_sites
