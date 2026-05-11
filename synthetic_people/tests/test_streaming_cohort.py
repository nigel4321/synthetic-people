"""Tests for the streaming-cohort primitives (PR 2 of option 3).

The streaming path replaces the materialised
``_tree_sequence_to_sites`` + in-place ``inject_*`` + ``sites.sort``
pipeline with a two-pass walk over the TreeSequence that holds only
~120 MB of light metadata at WGS scale. This file verifies the key
invariant: **byte-identical output at every fixed seed combination**.

If the streaming path ever diverges, downstream cohort BCFs would
differ from the materialised path even at the same ``--seed`` —
silently breaking reproducibility. These tests lock that down at the
generator level before the cli wiring (PR 3) routes production runs
through it.

Tests gate on ``msprime`` being importable. Without it the
TreeSequence fixtures can't be built; tests skip cleanly.
"""

from __future__ import annotations

import importlib.util
import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

_HAVE_MSPRIME = importlib.util.find_spec("msprime") is not None


def _build_tree_sequence(*, n_people: int = 4, length_bp: int = 200_000,
                        seed: int = 42, mu: float = 1e-6,
                        rec_rate: float = 1e-7):
    """Build a small TreeSequence fixture deterministic in the seed.

    Small enough to run in a few seconds; large enough to produce a
    few hundred variants so the overlay-injection paths actually have
    something to land on.
    """
    import msprime
    ts = msprime.sim_ancestry(
        samples=n_people,
        sequence_length=length_bp,
        recombination_rate=rec_rate,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(
        ts, rate=mu, random_seed=seed,
        model=msprime.BinaryMutationModel(),
    )
    return ts


def _materialised_path(ts, chrom: str, n_people: int, *,
                      rng_seed: int,
                      overlay_rng_seed: int | None = None,
                      clinvar_records=None,
                      clinvar_inject_density: float = 0.0,
                      rsid_pool=None,
                      rsid_density: float = 0.0,
                      cosmic_records=None,
                      cosmic_inject_density: float = 0.0):
    """Reproduce the cli's current per-chrom flow exactly: build
    sites in RAM, run annotate + inject overlays in order, sort by
    (chrom, pos). The streaming output must match this list exactly."""
    from syntheticgen.coalescent import _tree_sequence_to_sites
    from syntheticgen.titv import DEFAULT_TARGET_TITV
    from syntheticgen.clinvar import annotate_clinvar, inject_clinvar
    from syntheticgen.dbsnp import inject_rsids
    from syntheticgen.cosmic import inject_cosmic

    rng = random.Random(rng_seed)
    sites = _tree_sequence_to_sites(
        ts, chrom, n_people, rng, DEFAULT_TARGET_TITV,
    )
    # Match the streamer's contract: overlay_rng defaults to the main
    # rng when no overlay seed is provided. Without this default, the
    # injectors below would receive ``None`` and raise on
    # ``None.shuffle``/``sample`` if any injection density were > 0.
    overlay_rng = (random.Random(overlay_rng_seed)
                   if overlay_rng_seed is not None else rng)
    if clinvar_records:
        annotate_clinvar(sites, clinvar_records)
        if clinvar_inject_density > 0:
            inject_clinvar(sites, clinvar_records,
                           clinvar_inject_density, overlay_rng)
    clinvar_reserved = {i for i, s in enumerate(sites)
                        if s.get("clnsig")}
    if rsid_pool and rsid_density > 0:
        inject_rsids(sites, rsid_pool, rsid_density, overlay_rng,
                     reserve_indices=clinvar_reserved)
    all_overlay_reserved = {
        i for i, s in enumerate(sites)
        if s.get("clnsig") or str(s.get("id", "")).startswith("rs")
    }
    if cosmic_records and cosmic_inject_density > 0:
        inject_cosmic(sites, cosmic_records, cosmic_inject_density,
                      overlay_rng,
                      reserve_indices=all_overlay_reserved)
    sites.sort(key=lambda s: s["pos"])
    return sites


def _streamed_path(ts, chrom: str, n_people: int, *,
                  rng_seed: int,
                  overlay_rng_seed: int | None = None,
                  clinvar_records=None,
                  clinvar_inject_density: float = 0.0,
                  rsid_pool=None,
                  rsid_density: float = 0.0,
                  cosmic_records=None,
                  cosmic_inject_density: float = 0.0):
    """Drain the streaming generator into a list for comparison."""
    from syntheticgen.coalescent import stream_cohort_sites
    from syntheticgen.titv import DEFAULT_TARGET_TITV

    rng = random.Random(rng_seed)
    # Mirror the materialised helper above: default overlay_rng to the
    # main rng when no overlay seed is provided. ``stream_cohort_sites``
    # already does this internally; passing ``None`` here is harmless
    # but explicit is clearer.
    overlay_rng = (random.Random(overlay_rng_seed)
                   if overlay_rng_seed is not None else rng)
    return list(stream_cohort_sites(
        ts, chrom, n_people, rng,
        DEFAULT_TARGET_TITV,
        clinvar_records=clinvar_records,
        clinvar_inject_density=clinvar_inject_density,
        rsid_pool=rsid_pool,
        rsid_density=rsid_density,
        cosmic_records=cosmic_records,
        cosmic_inject_density=cosmic_inject_density,
        overlay_rng=overlay_rng,
    ))


def _clinvar_recs_on_chrom(chrom: str, positions, ref="A", alt="C"):
    return [{"chrom": chrom, "pos": p, "ref": ref, "alt": alt,
             "id": f"VCV{i:04d}", "clnsig": "Pathogenic",
             "clndn": f"Disease_{i}"}
            for i, p in enumerate(positions)]


def _rsid_recs_on_chrom(chrom: str, positions, ref="A", alt="C"):
    return [{"chrom": chrom, "pos": p, "ref": ref, "alt": alt,
             "rsid": f"rs{i+10001}"}
            for i, p in enumerate(positions)]


def _cosmic_recs_on_chrom(chrom: str, positions, ref="A", alt="C"):
    return [{"chrom": chrom, "pos": p, "ref": ref, "alt": alt,
             "id": f"COSV{i:04d}", "gene": f"GENE{i % 5}"}
            for i, p in enumerate(positions)]


@unittest.skipUnless(_HAVE_MSPRIME, "msprime not installed")
class TreeSequenceToSitesMetaTest(unittest.TestCase):
    """Pass 1 produces ``(chrom, pos, ref, alt, ac, n_haplotypes)``
    tuples that match the corresponding fields from the materialised
    site dicts at the same seed."""

    def test_meta_matches_materialised_per_field(self):
        from syntheticgen.coalescent import (
            _tree_sequence_to_sites,
            _tree_sequence_to_sites_meta,
        )
        from syntheticgen.titv import DEFAULT_TARGET_TITV
        ts = _build_tree_sequence(n_people=4, length_bp=100_000)

        rng_a = random.Random(7)
        sites = _tree_sequence_to_sites(
            ts, "22", 4, rng_a, DEFAULT_TARGET_TITV,
        )

        rng_b = random.Random(7)
        meta = _tree_sequence_to_sites_meta(
            ts, "22", 4, rng_b, DEFAULT_TARGET_TITV,
        )

        self.assertEqual(len(sites), len(meta),
                         "meta and full path must accept the same variants")
        for s, m in zip(sites, meta):
            chrom, pos, ref, alt, ac, n_haps = m
            self.assertEqual(s["chrom"], chrom)
            self.assertEqual(s["pos"], pos)
            self.assertEqual(s["ref"], ref)
            self.assertEqual(s["alts"], [alt])
            self.assertEqual(s["acs"], [ac])
            self.assertEqual(s["n_haplotypes"], n_haps)

    def test_meta_leaves_rng_in_same_state_as_full_path(self):
        # Both paths must consume the same number of rng draws.
        from syntheticgen.coalescent import (
            _tree_sequence_to_sites,
            _tree_sequence_to_sites_meta,
        )
        from syntheticgen.titv import DEFAULT_TARGET_TITV
        ts = _build_tree_sequence(n_people=4, length_bp=100_000)

        rng_a = random.Random(13)
        _tree_sequence_to_sites(
            ts, "22", 4, rng_a, DEFAULT_TARGET_TITV,
        )
        state_a = rng_a.getstate()

        rng_b = random.Random(13)
        _tree_sequence_to_sites_meta(
            ts, "22", 4, rng_b, DEFAULT_TARGET_TITV,
        )
        state_b = rng_b.getstate()

        self.assertEqual(state_a, state_b,
                         "meta walk must consume rng identically to "
                         "the materialised walk for byte-identical "
                         "downstream reproducibility")


@unittest.skipUnless(_HAVE_MSPRIME, "msprime not installed")
class StreamCohortSitesParityTest(unittest.TestCase):
    """Top-level: streaming path produces byte-identical sites to the
    materialised path under every combination of overlay settings."""

    def setUp(self):
        self.ts = _build_tree_sequence(
            n_people=6, length_bp=200_000, seed=42,
        )
        self.chrom = "22"
        self.n_people = 6

    def _assert_byte_identical(self, materialised, streamed):
        self.assertEqual(
            len(materialised), len(streamed),
            f"length differs: materialised={len(materialised)}, "
            f"streamed={len(streamed)}",
        )
        for i, (m, s) in enumerate(zip(materialised, streamed)):
            self.assertEqual(
                m, s,
                f"site {i} diverges:\n  materialised={m}\n  streamed={s}",
            )

    def test_no_overlays(self):
        for seed in (0, 1, 42, 100):
            with self.subTest(seed=seed):
                m = _materialised_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed,
                )
                s = _streamed_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed,
                )
                self._assert_byte_identical(m, s)

    def test_clinvar_annotation_only(self):
        # Pick a few positions that the tree variants are likely to
        # land near; we don't need exact matches, just enough records
        # to exercise the annotation map code path.
        meta_recs = _clinvar_recs_on_chrom(
            self.chrom, range(1_000, 200_000, 2_000),
        )
        for seed in (0, 1, 42):
            with self.subTest(seed=seed):
                m = _materialised_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed, overlay_rng_seed=seed + 100,
                    clinvar_records=meta_recs,
                )
                s = _streamed_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed, overlay_rng_seed=seed + 100,
                    clinvar_records=meta_recs,
                )
                self._assert_byte_identical(m, s)

    def test_clinvar_inject(self):
        recs = _clinvar_recs_on_chrom(
            self.chrom, [50_000, 75_000, 100_000, 125_000, 150_000],
        )
        for seed in (0, 1, 42):
            with self.subTest(seed=seed):
                m = _materialised_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed, overlay_rng_seed=seed + 200,
                    clinvar_records=recs,
                    clinvar_inject_density=0.05,
                )
                s = _streamed_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed, overlay_rng_seed=seed + 200,
                    clinvar_records=recs,
                    clinvar_inject_density=0.05,
                )
                self._assert_byte_identical(m, s)

    def test_rsid_inject(self):
        rsids = _rsid_recs_on_chrom(
            self.chrom, [25_000, 60_000, 90_000, 140_000, 175_000],
        )
        for seed in (0, 1, 42):
            with self.subTest(seed=seed):
                m = _materialised_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed, overlay_rng_seed=seed + 300,
                    rsid_pool=rsids, rsid_density=0.05,
                )
                s = _streamed_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed, overlay_rng_seed=seed + 300,
                    rsid_pool=rsids, rsid_density=0.05,
                )
                self._assert_byte_identical(m, s)

    def test_cosmic_inject(self):
        cosmic = _cosmic_recs_on_chrom(
            self.chrom, [30_000, 65_000, 105_000, 155_000],
        )
        for seed in (0, 1, 42):
            with self.subTest(seed=seed):
                m = _materialised_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed, overlay_rng_seed=seed + 400,
                    cosmic_records=cosmic,
                    cosmic_inject_density=0.05,
                )
                s = _streamed_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed, overlay_rng_seed=seed + 400,
                    cosmic_records=cosmic,
                    cosmic_inject_density=0.05,
                )
                self._assert_byte_identical(m, s)

    def test_all_overlays_chained(self):
        # The full cli-equivalent chain: ClinVar annotate + inject,
        # rsID inject (with reserve), COSMIC inject (with reserve).
        clinvar = _clinvar_recs_on_chrom(
            self.chrom, range(1_000, 200_000, 2_500),
        )
        rsids = _rsid_recs_on_chrom(
            self.chrom, range(2_500, 200_000, 3_000),
        )
        cosmic = _cosmic_recs_on_chrom(
            self.chrom, range(4_000, 200_000, 5_000),
        )
        for seed in (0, 1, 42):
            with self.subTest(seed=seed):
                m = _materialised_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed, overlay_rng_seed=seed + 500,
                    clinvar_records=clinvar,
                    clinvar_inject_density=0.04,
                    rsid_pool=rsids, rsid_density=0.05,
                    cosmic_records=cosmic,
                    cosmic_inject_density=0.03,
                )
                s = _streamed_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed, overlay_rng_seed=seed + 500,
                    clinvar_records=clinvar,
                    clinvar_inject_density=0.04,
                    rsid_pool=rsids, rsid_density=0.05,
                    cosmic_records=cosmic,
                    cosmic_inject_density=0.03,
                )
                self._assert_byte_identical(m, s)


@unittest.skipUnless(_HAVE_MSPRIME, "msprime not installed")
class RngStateAfterStreamTest(unittest.TestCase):
    """Downstream of the cohort-write step, the master rng continues
    to be consumed (e.g. per-person VCF derivation if the cli routed
    through it directly). The streaming path must leave the rng in
    the same state as the materialised path so any downstream rng
    consumer produces the same draws."""

    def test_rng_state_after_stream_matches_materialised(self):
        from syntheticgen.coalescent import (
            _tree_sequence_to_sites,
            stream_cohort_sites,
        )
        from syntheticgen.titv import DEFAULT_TARGET_TITV
        ts = _build_tree_sequence(n_people=4, length_bp=100_000)

        # Materialised: consumes rng once.
        rng_a = random.Random(99)
        _tree_sequence_to_sites(
            ts, "22", 4, rng_a, DEFAULT_TARGET_TITV,
        )
        state_a = rng_a.getstate()

        # Streamed: must end in the same state.
        rng_b = random.Random(99)
        list(stream_cohort_sites(
            ts, "22", 4, rng_b, DEFAULT_TARGET_TITV,
        ))
        state_b = rng_b.getstate()

        self.assertEqual(state_a, state_b)


@unittest.skipUnless(_HAVE_MSPRIME, "msprime not installed")
class StreamPosMonotoneTest(unittest.TestCase):
    """The streamer's primary contract is that emitted sites are in
    pos-sorted order — that's what Arrow batches and bcftools merge
    require. Verify even when overlays scramble tree-walk order."""

    def test_pos_monotone_with_overlays(self):
        ts = _build_tree_sequence(n_people=4, length_bp=200_000)
        # ClinVar records spread across the entire range
        # — guarantees some injections at positions that differ
        # significantly from the tree-walk position they replace.
        clinvar = _clinvar_recs_on_chrom(
            "22", range(500, 200_000, 1_000),
        )
        from syntheticgen.coalescent import stream_cohort_sites
        from syntheticgen.titv import DEFAULT_TARGET_TITV
        rng = random.Random(1)
        overlay_rng = random.Random(2)
        sites = list(stream_cohort_sites(
            ts, "22", 4, rng, DEFAULT_TARGET_TITV,
            clinvar_records=clinvar,
            clinvar_inject_density=0.30,  # heavy injection
            overlay_rng=overlay_rng,
        ))
        positions = [s["pos"] for s in sites]
        self.assertEqual(positions, sorted(positions),
                         "streaming path must emit sites in pos-sorted "
                         "order even under heavy overlay injection")


@unittest.skipUnless(_HAVE_MSPRIME, "msprime not installed")
class StreamPosSortWithOutOfOrderOverlaysRegressionTest(unittest.TestCase):
    """Regression test for the PR #51 Copilot review finding: the
    safe-yield bound in ``_stream_cohort_pass2`` looked up
    ``overlay_positions_sorted[overlay_next_idx]`` where the array was
    sorted by *position* but the index counted overlays in tree-walk
    order. The look-up could report a future-min larger than the
    actual smallest un-encountered overlay pos, letting tree sites
    yield ahead of overlay sites and breaking pos-sortedness.

    The fix uses a suffix-min over overlays sorted by *accepted_idx*.
    This test constructs an inject_map that explicitly trips the
    pre-fix algorithm (high-pos overlay at early idx, low-pos overlay
    at later idx) and asserts the output is still pos-sorted."""

    def test_pos_sorted_when_first_overlay_high_later_overlay_low(self):
        from syntheticgen.coalescent import (
            _tree_sequence_to_sites_meta,
            _stream_cohort_pass2,
        )
        from syntheticgen.titv import DEFAULT_TARGET_TITV

        # Seed 4 with this fixture shape was the concrete repro
        # produced during the PR #51 review investigation.
        ts = _build_tree_sequence(
            n_people=4, length_bp=300_000, seed=4,
        )
        rng = random.Random(4)
        meta = _tree_sequence_to_sites_meta(
            ts, "22", 4, rng, DEFAULT_TARGET_TITV,
        )
        self.assertGreaterEqual(len(meta), 4,
            "fixture must produce at least 4 sites to exercise the "
            "high-then-low overlay-pos pattern")

        # Synthetic inject_map: idx 1 -> very high pos, idx 3 -> low pos.
        # The pre-fix algorithm would yield the tree site at meta[2]'s
        # pos before seeing the idx-3 overlay at a smaller pos.
        tree_positions = [m[1] for m in meta]
        high_pos = tree_positions[3] + 1000
        low_pos = max(1, tree_positions[1] - 5)
        inject_map = {
            1: {"pos": high_pos, "ref": "A", "alts": ["G"], "id": "."},
            3: {"pos": low_pos, "ref": "A", "alts": ["G"], "id": "."},
        }

        rng2 = random.Random(4)
        sites = list(_stream_cohort_pass2(
            ts, "22", 4, rng2, DEFAULT_TARGET_TITV,
            meta, inject_map, {},
        ))
        positions = [s["pos"] for s in sites]
        self.assertEqual(
            positions, sorted(positions),
            f"streaming output must remain pos-sorted under out-of-"
            f"order overlay positions; got {positions}",
        )


@unittest.skipUnless(_HAVE_MSPRIME, "msprime not installed")
class StreamCohortChainedOverlayParityAtScaleTest(unittest.TestCase):
    """Regression test for the PR #51 Copilot review finding on
    chained overlays.

    The materialised path sorts the sites list inside ``inject_*``
    after each mutation, so the next planner's ``used_keys`` and
    ``reserve_indices`` come from the post-sort view. PR #51's first
    draft of the streaming pipeline ran all planners against the
    original tree-walk-order ``sites_meta``, which diverged on
    chained overlays at large enough scale.

    The previous small-fixture parity tests (n=6, 200 kb, density
    0.04) failed to expose this — they typically produced only ~2
    injections per overlay, too few to trigger a chained divergence.
    These tests run at n=20, 1 Mb, density 0.10 across multiple seeds
    so the bug, if reintroduced, fails loudly."""

    def setUp(self):
        self.ts = _build_tree_sequence(
            n_people=20, length_bp=1_000_000, seed=42,
        )
        self.chrom = "22"
        self.n_people = 20
        self.clinvar = _clinvar_recs_on_chrom(
            self.chrom, range(50_000, 1_000_000, 7_500),
        )
        self.rsids = _rsid_recs_on_chrom(
            self.chrom, range(75_000, 1_000_000, 6_500),
        )
        self.cosmic = _cosmic_recs_on_chrom(
            self.chrom, range(100_000, 1_000_000, 8_000),
        )

    def _assert_byte_identical(self, materialised, streamed):
        self.assertEqual(
            len(materialised), len(streamed),
            f"length differs: m={len(materialised)} s={len(streamed)}",
        )
        for i, (m, s) in enumerate(zip(materialised, streamed)):
            self.assertEqual(
                m, s,
                f"site {i} diverges:\n  materialised={m}\n  streamed={s}",
            )

    def test_chained_overlays_at_density_10pct(self):
        for seed in range(5):
            with self.subTest(seed=seed):
                m = _materialised_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed, overlay_rng_seed=seed + 1000,
                    clinvar_records=self.clinvar,
                    clinvar_inject_density=0.10,
                    rsid_pool=self.rsids, rsid_density=0.10,
                    cosmic_records=self.cosmic,
                    cosmic_inject_density=0.10,
                )
                s = _streamed_path(
                    self.ts, self.chrom, self.n_people,
                    rng_seed=seed, overlay_rng_seed=seed + 1000,
                    clinvar_records=self.clinvar,
                    clinvar_inject_density=0.10,
                    rsid_pool=self.rsids, rsid_density=0.10,
                    cosmic_records=self.cosmic,
                    cosmic_inject_density=0.10,
                )
                self._assert_byte_identical(m, s)


if __name__ == "__main__":
    unittest.main()
