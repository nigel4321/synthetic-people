"""Cohort-streaming resume contract (Phase 5b2).

Long cohort runs (n=100k × 22 chromosomes × full-chr length) take
hours; an OOM kill or SIGINT mid-stream should not require redoing
hours of completed work. This module persists a tiny ``cohort.meta.json``
file alongside the cohort BCF directory recording:

* the parameter set the run was started with (seed, n, build, chroms,
  demo model, population, chr_length, rec rate, mu) — used to detect
  resume mismatches;
* the sample IDs and per-person seeds drawn at the start of the run
  (so the resumed run uses the exact same values rather than re-
  drawing from a master rng that's now in a different state);
* per-chromosome overlay seeds (so each chromosome's overlay rng is
  independent of run order and resumes deterministically);
* the list of chromosomes whose cohort BCF has finished writing
  (so resume skips them).

On startup the streamed pipeline calls :func:`load_or_create_meta`,
which returns either a loaded ``Resume`` object (when the existing
meta.json matches the new run's params) or a freshly-derived one
written to disk. Each completed chromosome subsequently calls
:func:`mark_chromosome_completed` to update the on-disk record.

``--no-resume`` forces a fresh start: any existing meta.json plus
the entire cohort/ directory get wiped before deriving new seeds,
so a partial run from previous parameters can't accidentally bleed
into a re-run.
"""

from __future__ import annotations

import json
import random
import shutil
from dataclasses import dataclass, field
from pathlib import Path


# Schema version 2 (2026-05-15): added ``sexes`` for M13.1 per-person
# sex assignment. Schema-version 1 meta files surface as
# ResumeMismatch — the user runs once with --no-resume to migrate.
_SCHEMA_VERSION = 2


def _params_for(args, chromosomes: list) -> dict:
    """Extract the resume-identity-relevant params from argparse args.

    Two runs are considered the same cohort if these match exactly.
    Anything that affects the simulated content (seed, demography,
    chromosome span) is in here; cosmetic args (output_dir,
    --workers, --mode, overlay densities) are deliberately not —
    overlays are seeded per-chromosome so changing density mid-resume
    would be detected by an overlay-density mismatch separately if we
    cared, but for now we treat density as runtime-mutable.
    """
    return {
        "seed": args.seed,
        "n_people": args.n,
        "build": args.build,
        "chromosomes": list(chromosomes),
        "chr_length_mb": args.chr_length_mb,
        "demo_model": args.demo_model,
        "population": args.population,
        "rec_rate": args.rec_rate,
        "mu": args.mu,
        # M13.1: male_fraction is part of the cohort identity — changing
        # it between runs would silently reuse the old sex assignments
        # from the persisted meta. Surfacing it here triggers
        # ResumeMismatch when it changes, forcing the user to consciously
        # choose --no-resume.
        "male_fraction": args.male_fraction,
    }


def _draw_sexes(seed, n: int, male_fraction: float) -> list[str]:
    """Draw ``n`` per-person sex assignments without consuming the
    master rng.

    M13.1 contract: M13.1 records per-person sex but does not change
    any simulation behaviour at the same ``--seed``. To honour that,
    sex draws must NOT advance the master rng — doing so would shift
    every downstream rng consumer (overlay seeds, per-chrom seeds,
    error model, etc.) and a fixed-seed run would no longer reproduce
    pre-M13.1 output.

    We seed a dedicated ``random.Random`` deterministically from the
    master seed (and a fixed salt so the sex rng never collides with
    any other derived rng in this codebase). The master rng is
    untouched.

    Returns a list of ``"m"`` / ``"f"`` strings, deterministic given
    ``(seed, n, male_fraction)``. When ``seed`` is None we fall back
    to ``0`` so the sex draws stay reproducible across runs that omit
    ``--seed`` (matching the existing "no-seed → constant rng state"
    behaviour for those debug-only runs).
    """
    salt = 0x5E_5E_5E_5E  # arbitrary fixed salt — "SE" for "sex"
    sex_rng = random.Random((seed or 0) ^ salt)
    return [
        "m" if sex_rng.random() < male_fraction else "f"
        for _ in range(n)
    ]


@dataclass
class Resume:
    """Loaded or freshly-derived resume state for a streamed run."""

    meta_path: Path
    params: dict
    samples: list
    person_seeds: list
    overlay_seeds: dict       # chrom -> int
    # M13.1: per-person sex assignment, parallel-indexed to ``samples``.
    # Persisted so a resumed run sees the same sexes the original run
    # drew, same contract as person_seeds.
    sexes: list = field(default_factory=list)
    completed_chromosomes: list = field(default_factory=list)

    def is_chromosome_done(self, chrom: str) -> bool:
        return chrom in self.completed_chromosomes

    def mark_chromosome_done(self, chrom: str) -> None:
        if chrom in self.completed_chromosomes:
            return
        self.completed_chromosomes.append(chrom)
        self._save()

    def _save(self) -> None:
        payload = {
            "schema_version": _SCHEMA_VERSION,
            "params": self.params,
            "samples": self.samples,
            "person_seeds": self.person_seeds,
            "sexes": self.sexes,
            "overlay_seeds": self.overlay_seeds,
            "completed_chromosomes": self.completed_chromosomes,
        }
        # Atomic-ish write: stage to a sibling tmp file then rename.
        # On POSIX a same-directory rename is atomic, so the on-disk
        # meta.json is always either the previous version or the new
        # one — never a torn partial write.
        tmp = self.meta_path.with_suffix(self.meta_path.suffix + ".tmp")
        with open(tmp, "w") as fh:
            json.dump(payload, fh, indent=2)
        tmp.replace(self.meta_path)


class ResumeMismatch(Exception):
    """Existing cohort.meta.json was for a different cohort.

    Surfaces in the CLI as a clear error with a hint about
    ``--no-resume``. Caller is expected to convert into ``sys.exit``
    with a user-facing message rather than letting the exception
    propagate.
    """


def load_or_create_meta(args, chromosomes: list, cohort_dir: Path,
                        rng: random.Random,
                        force_fresh: bool = False) -> Resume:
    """Load an existing resume record or create a new one.

    When an existing ``cohort.meta.json`` is present and its params
    match the current run, the loaded record is returned (with
    completed_chromosomes carried through so the streamed loop can
    skip them). Mismatched params raise :exc:`ResumeMismatch`. When
    ``force_fresh`` is True any existing meta.json plus the cohort/
    directory get wiped first; a freshly-derived record is then
    written to disk.

    Sample IDs, per-person seeds, and per-chromosome overlay seeds are
    drawn from ``rng`` only on the create path, so a resumed run
    sees the exact same values as the original run did.
    """
    meta_path = cohort_dir / "cohort.meta.json"
    cohort_dir.mkdir(parents=True, exist_ok=True)

    if force_fresh:
        # --no-resume: nuke any prior state. We delete everything
        # under cohort/ rather than just the meta file so a stale
        # cohort.chr*.bcf from a previous (mismatched-param) run
        # can't accidentally be reused.
        for path in cohort_dir.iterdir():
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                shutil.rmtree(path)
    elif meta_path.is_file():
        with open(meta_path) as fh:
            payload = json.load(fh)
        if payload.get("schema_version") != _SCHEMA_VERSION:
            raise ResumeMismatch(
                f"cohort.meta.json schema version "
                f"{payload.get('schema_version')!r} is not supported "
                f"(expected {_SCHEMA_VERSION}). Use --no-resume to "
                f"start fresh."
            )
        params = _params_for(args, chromosomes)
        if payload.get("params") != params:
            mismatched = [
                k for k in params
                if payload.get("params", {}).get(k) != params[k]
            ]
            raise ResumeMismatch(
                f"cohort.meta.json was written with different "
                f"parameters: {mismatched} mismatch. Use a different "
                f"--output-dir, or pass --no-resume to wipe the "
                f"existing cohort/ directory and start fresh."
            )
        # Everything checks out; reuse the prior state.
        return Resume(
            meta_path=meta_path,
            params=payload["params"],
            samples=payload["samples"],
            person_seeds=payload["person_seeds"],
            sexes=payload["sexes"],
            overlay_seeds=payload["overlay_seeds"],
            completed_chromosomes=payload.get(
                "completed_chromosomes", []),
        )

    # Fresh start: derive new state from the master rng.
    from .background import draw_sample_ids   # local import; avoids cycle
    samples = draw_sample_ids(args.n, rng)
    person_seeds = [rng.randint(1, 2**31 - 1) for _ in range(args.n)]
    # M13.1: sexes are drawn from a SEPARATE rng (see ``_draw_sexes``).
    # Drawing from ``rng`` here would advance the master rng state and
    # shift every downstream consumer, breaking the "no simulation
    # change at fixed seed" contract.
    sexes = _draw_sexes(args.seed, args.n, args.male_fraction)
    overlay_seeds = {
        chrom: rng.randint(1, 2**31 - 1)
        for chrom in chromosomes
    }
    resume = Resume(
        meta_path=meta_path,
        params=_params_for(args, chromosomes),
        samples=samples,
        person_seeds=person_seeds,
        sexes=sexes,
        overlay_seeds=overlay_seeds,
        completed_chromosomes=[],
    )
    resume._save()
    return resume
