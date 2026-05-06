"""UK-cohort admixture simulator with local ancestry truth tracking.

Builds a `demes`-defined demography in which three source demes (EUR,
SAS, AFR) contribute, via a single admixture pulse `PULSE_TIME`
generations ago, into a "UK" deme that the cohort is sampled from.
`msprime.sim_ancestry` runs with `record_migrations=True`, so for each
haplotype we can walk the tree back to the lineage node spanning the
pulse time and look up which source deme it migrated into. That gives
per-haplotype local ancestry segments — the truth set written to
`out/ancestry/person_<N>.bed`.

Mutations are drawn with `BinaryMutationModel` and REF/ALT bases are
synthesised through the M3+ Ti/Tv calibrator, mirroring the M5 path.
The output `sites` list has the same shape M4/M5 produce, so the
writer and CLI per-person loop are unchanged.
"""

from __future__ import annotations

import multiprocessing as mp
import random
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from .builds import BUILDS
from .titv import DEFAULT_TARGET_TITV, choose_alt


DEFAULT_EUR_FRAC = 0.60
DEFAULT_SAS_FRAC = 0.25
DEFAULT_AFR_FRAC = 0.15

# Pulse 20 generations back ≈ 600 years — covers post-Industrial-Revolution
# plus modern-era migration into the UK. Configurable via build_uk_demography.
PULSE_TIME = 20.0

SOURCE_POPS = ("EUR", "SAS", "AFR")


def _require_deps() -> None:
    try:
        import msprime  # noqa: F401
        import demes  # noqa: F401
        import tskit  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "Admixture path requires msprime + demes + tskit. Install via "
            "`pip install -r synthetic_people/requirements.txt`."
        ) from exc


def build_uk_demography(eur_frac: float = DEFAULT_EUR_FRAC,
                        sas_frac: float = DEFAULT_SAS_FRAC,
                        afr_frac: float = DEFAULT_AFR_FRAC,
                        pulse_time: float = PULSE_TIME,
                        uk_size: int = 50_000):
    """Return a `demes.Graph` for an EUR+SAS+AFR → UK admixture pulse.

    Population sizes mirror the Gutenkunst out-of-Africa parameterisation
    (ANC≈12.3k, OOA bottleneck≈2.1k, AFR≈12.3k, EUR/SAS present-day≈10k).
    SAS branches from OOA at the same time as EUR — a deliberate
    simplification: we are not trying to match a particular published
    EUR/SAS split estimate, only to give the two demes their own
    coalescent histories before the admixture pulse.
    """
    _require_deps()
    import demes

    total = eur_frac + sas_frac + afr_frac
    if abs(total - 1.0) > 1e-6:
        raise ValueError(
            f"ancestry proportions must sum to 1.0; got "
            f"EUR={eur_frac}, SAS={sas_frac}, AFR={afr_frac} "
            f"(sum={total:.6f})"
        )
    if min(eur_frac, sas_frac, afr_frac) < 0:
        raise ValueError("ancestry proportions must be non-negative")
    if pulse_time <= 0:
        raise ValueError("pulse_time must be > 0 generations")

    b = demes.Builder(time_units="generations")
    b.add_deme("ANC", epochs=[dict(start_size=12_300, end_time=5_920)])
    b.add_deme("AFR", ancestors=["ANC"],
               epochs=[dict(start_size=12_300, end_time=0)])
    b.add_deme("OOA", ancestors=["ANC"], start_time=5_920,
               epochs=[dict(start_size=2_100, end_time=2_040)])
    b.add_deme("EUR", ancestors=["OOA"],
               epochs=[dict(start_size=10_000, end_time=0)])
    b.add_deme("SAS", ancestors=["OOA"],
               epochs=[dict(start_size=10_000, end_time=0)])
    b.add_deme("UK", start_time=pulse_time,
               ancestors=list(SOURCE_POPS),
               proportions=[eur_frac, sas_frac, afr_frac],
               epochs=[dict(start_size=uk_size, end_time=0)])
    return b.resolve()


def simulate_chromosome(chrom: str, build: str, n_people: int,
                        length_mb: float, proportions: tuple,
                        rec_rate: float, mu: float,
                        rng: random.Random,
                        titv_target: float = DEFAULT_TARGET_TITV,
                        pulse_time: float = PULSE_TIME) -> tuple:
    """Simulate one chromosome of UK admixture.

    Returns ``(sites, person_segments)`` where ``person_segments[i]`` is
    a list of ``(start, end, hap1_pop, hap2_pop)`` tuples for person
    ``i`` on this chromosome, with positions in the chromosome's
    integer 0-based half-open coordinate space.
    """
    _require_deps()
    import msprime

    if chrom not in BUILDS[build]["contigs"]:
        raise ValueError(f"unknown chromosome {chrom!r} for build {build}")
    chrom_length = BUILDS[build]["contigs"][chrom]
    sim_length = int(chrom_length if length_mb <= 0
                     else min(chrom_length, length_mb * 1_000_000))

    eur, sas, afr = proportions
    graph = build_uk_demography(eur, sas, afr, pulse_time=pulse_time)
    demog = msprime.Demography.from_demes(graph)

    seed_anc = rng.randint(1, 2**31 - 1)
    seed_mut = rng.randint(1, 2**31 - 1)

    ts = msprime.sim_ancestry(
        samples={"UK": n_people},
        demography=demog,
        sequence_length=sim_length,
        recombination_rate=rec_rate,
        random_seed=seed_anc,
        record_migrations=True,
    )
    ts = msprime.sim_mutations(
        ts, rate=mu, random_seed=seed_mut,
        model=msprime.BinaryMutationModel(),
    )

    sites = _tree_sequence_to_sites(ts, chrom, n_people, rng, titv_target)
    person_segs = _local_ancestry(ts, n_people, pulse_time)
    return sites, person_segs


def _tree_sequence_to_sites(ts, chrom: str, n_people: int,
                            rng: random.Random,
                            titv_target: float) -> list:
    """Phase 5c: emits sparse carriers, matching coalescent.py.

    The admixture path uses BinaryMutationModel (only allele indices
    {0, 1}), so the multi-allelic skip in coalescent's converter is
    unnecessary here — but the rest of the structure stays in sync
    so a downstream consumer doesn't have to know which path
    produced a given site.
    """
    n_haplotypes = 2 * n_people
    sites: list = []
    used: set = set()
    for var in ts.variants():
        pos = int(var.site.position) + 1
        while pos in used:
            pos += 1
        used.add(pos)

        gts_arr = var.genotypes
        nalt = int((gts_arr > 0).sum())
        if nalt == 0 or nalt == n_haplotypes:
            continue

        ref = rng.choice(("A", "C", "G", "T"))
        alt = choose_alt(ref, rng, target=titv_target)
        assert alt is not None

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
            "afs": [nalt / n_haplotypes],
            "acs": [nalt],
            "n_haplotypes": n_haplotypes,
            "carriers": carriers,
        })
    return sites


def _local_ancestry(ts, n_people: int, pulse_time: float) -> list:
    """Return per-person list of ``(start, end, h1_pop, h2_pop)`` segments.

    Walks each haplotype-sample's lineage in the tree at every breakpoint
    until it finds the lineage node spanning ``pulse_time``; the t=pulse
    migration on that node tells us which source deme the lineage came
    from at this position. Adjacent same-ancestry segments are merged.
    Then haplotypes are paired into per-person joint intervals.
    """
    import tskit

    pop_names = [p.metadata["name"] for p in ts.populations()]
    mig_by_node = defaultdict(list)
    for m in ts.migrations():
        if m.time == pulse_time:
            mig_by_node[m.node].append((m.left, m.right, m.dest))

    samples = ts.samples()

    def lineage_pop(tree, sample_node: int, position: float) -> str:
        u = sample_node
        while True:
            p = tree.parent(u)
            if p == tskit.NULL:
                return pop_names[ts.node(u).population]
            if ts.node(p).time > pulse_time:
                break
            u = p
        for left, right, dest in mig_by_node.get(u, ()):
            if left <= position < right:
                return pop_names[dest]
        return pop_names[ts.node(u).population]

    hap_segs: list = [[] for _ in range(2 * n_people)]
    for tree in ts.trees():
        L = float(tree.interval.left)
        R = float(tree.interval.right)
        mid = (L + R) / 2
        for hap_i, sample_node in enumerate(samples):
            pop = lineage_pop(tree, sample_node, mid)
            segs = hap_segs[hap_i]
            if segs and segs[-1][2] == pop and segs[-1][1] == L:
                segs[-1] = (segs[-1][0], R, pop)
            else:
                segs.append((L, R, pop))

    person_segs: list = []
    for i in range(n_people):
        joint = _intersect_haplotype_segments(hap_segs[2 * i],
                                              hap_segs[2 * i + 1])
        person_segs.append(joint)
    return person_segs


def _intersect_haplotype_segments(h1: list, h2: list) -> list:
    """Merge two per-haplotype segment lists into joint per-person rows.

    Both inputs share the same total interval [0, L); output rows are
    ``(start, end, h1_pop, h2_pop)`` with integer coordinates suitable
    for BED output.
    """
    out: list = []
    i = j = 0
    while i < len(h1) and j < len(h2):
        a_s, a_e, a_p = h1[i]
        b_s, b_e, b_p = h2[j]
        s = max(a_s, b_s)
        e = min(a_e, b_e)
        if s < e:
            row = (int(s), int(e), a_p, b_p)
            if out and out[-1][2] == a_p and out[-1][3] == b_p \
                    and out[-1][1] == row[0]:
                out[-1] = (out[-1][0], row[1], a_p, b_p)
            else:
                out.append(row)
        if a_e <= b_e:
            i += 1
        else:
            j += 1
    return out


def _simulate_chromosome_from_seed(chrom: str, build: str, n_people: int,
                                   length_mb: float, proportions: tuple,
                                   rec_rate: float, mu: float, seed: int,
                                   titv_target: float,
                                   pulse_time: float) -> tuple:
    """Worker entry point — builds its own rng from the given seed."""
    rng = random.Random(seed)
    return simulate_chromosome(
        chrom, build, n_people, length_mb, proportions,
        rec_rate, mu, rng, titv_target, pulse_time=pulse_time,
    )


def simulate_cohort(chromosomes: list, build: str, n_people: int,
                    length_mb: float, proportions: tuple,
                    rec_rate: float, mu: float, rng: random.Random,
                    titv_target: float = DEFAULT_TARGET_TITV,
                    pulse_time: float = PULSE_TIME,
                    verbose: bool = False,
                    workers: int = 1) -> tuple:
    """Simulate the UK admixed cohort across one or more chromosomes.

    Returns ``(sites, ancestry)`` where ``ancestry[i]`` is a list of
    ``(chrom, start, end, h1_pop, h2_pop)`` rows for person ``i``,
    spanning all chromosomes in input order.

    With ``workers > 1`` and more than one chromosome, each chromosome
    runs in its own process. Per-chromosome seeds are pre-derived so
    output is deterministic regardless of the worker count.
    """
    eur, sas, afr = proportions
    if verbose:
        print(
            f"UK admixture proportions: EUR={eur:.2f} SAS={sas:.2f} "
            f"AFR={afr:.2f} (pulse {pulse_time:g} gens ago)",
            file=sys.stderr,
        )
    seeds = [rng.randint(1, 2**31 - 1) for _ in chromosomes]

    use_pool = workers > 1 and len(chromosomes) > 1
    all_sites: list = []
    ancestry: list = [[] for _ in range(n_people)]

    def _accumulate(chrom: str, sites: list, person_segs: list) -> None:
        all_sites.extend(sites)
        for i, segs in enumerate(person_segs):
            for L, R, h1, h2 in segs:
                ancestry[i].append((chrom, L, R, h1, h2))

    if not use_pool:
        for chrom, seed in zip(chromosomes, seeds):
            if verbose:
                print(f"  simulating chrom {chrom} (length {length_mb} "
                      f"Mb, UK admixture)...", file=sys.stderr)
            sites, person_segs = _simulate_chromosome_from_seed(
                chrom, build, n_people, length_mb, proportions,
                rec_rate, mu, seed, titv_target, pulse_time,
            )
            if verbose:
                print(f"    {len(sites)} variable sites on chrom "
                      f"{chrom}", file=sys.stderr)
            _accumulate(chrom, sites, person_segs)
        return all_sites, ancestry

    if verbose:
        print(f"  simulating {len(chromosomes)} chromosomes "
              f"in parallel with {workers} workers (UK admixture)...",
              file=sys.stderr)

    ctx = mp.get_context("fork")
    futures = []
    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as ex:
        for chrom, seed in zip(chromosomes, seeds):
            futures.append((chrom, ex.submit(
                _simulate_chromosome_from_seed,
                chrom, build, n_people, length_mb, proportions,
                rec_rate, mu, seed, titv_target, pulse_time,
            )))
        for chrom, fut in futures:
            sites, person_segs = fut.result()
            if verbose:
                print(f"    {len(sites)} variable sites on chrom "
                      f"{chrom}", file=sys.stderr)
            _accumulate(chrom, sites, person_segs)
    return all_sites, ancestry


def write_ancestry_bed(path: Path, segments: list) -> None:
    """Write one person's ancestry truth as BED.

    Columns: ``chrom\\tstart\\tend\\thap1_pop\\thap2_pop``. Coordinates
    follow BED convention (0-based half-open).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for chrom, start, end, h1, h2 in segments:
            f.write(f"{chrom}\t{start}\t{end}\t{h1}\t{h2}\n")


def ancestry_fractions(segments: list) -> dict:
    """Return ``{pop: fraction}`` averaged across both haplotypes by bp.

    Fractions are normalised against the total haplotype-bp covered by
    the segments (i.e. 2 × span per row). Source pops are always present
    in the output, with 0.0 if absent from the segments. Ancient demes
    (ANC, OOA) appear if any lineage has not yet found a SOURCE_POPS
    ancestor by `pulse_time` — should be rare with the default
    parameterisation and is reported faithfully when it happens.
    """
    totals: dict = defaultdict(float)
    total_bp = 0.0
    for _chrom, start, end, h1, h2 in segments:
        span = end - start
        totals[h1] += span
        totals[h2] += span
        total_bp += 2 * span
    keys = sorted(set(SOURCE_POPS) | set(totals))
    if total_bp == 0:
        return {p: 0.0 for p in keys}
    return {p: totals.get(p, 0.0) / total_bp for p in keys}
