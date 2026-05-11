"""COSMIC overlay (M7, optional / `--somatic`).

COSMIC requires a registration-gated download (no public anonymous
URL), so this module never auto-fetches. The user supplies a local
path to a COSMIC VCF (e.g. `Cosmic_GenomeScreensMutant_v100_GRCh38.vcf`
or any rebuild that follows the same INFO schema), and we overlay or
inject records onto cohort sites the same way ClinVar does.

If `--somatic` is set without a path, we emit a clear instruction and
exit non-zero — the spec mandates the optional gating, but silently
skipping is too quiet for a flag the user explicitly opted into.

The relevant COSMIC INFO tags vary slightly by release; we read what
exists and skip what doesn't:
- `GENE` / `GENE_NAME` — affected gene symbol
- `LEGACY_ID` / `COSMIC_ID` — COSV/COSM identifier (also in the ID col)
- `CDS` / `AA` — coding-DNA / protein consequence string

Fields land in `INFO/COSMIC_GENE` and `INFO/COSMIC_ID` on the output
records; the writer/header pick those up via the M7 INFO declarations.
"""

from __future__ import annotations

import random
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import Any


DEFAULT_COSMIC_INJECT_DENSITY = 0.005


def load_cosmic_records(cosmic_vcf: Path,
                        chromosomes: list[str],
                        max_per_chrom: int | None = 5000,
                        ) -> list[dict]:
    """Stream COSMIC entries restricted to `chromosomes`.

    Returns one dict per (chrom, pos, ref, alt) row. Robust to missing
    optional INFO tags — older COSMIC releases don't all carry every
    field. Multi-allelic ALTs are split.
    """
    out: list[dict] = []
    for chrom in chromosomes:
        # Permissive query — COSMIC schema varies; we extract whatever
        # is there and fall back to "." for missing tags.
        fmt = ("%CHROM\t%POS\t%ID\t%REF\t%ALT\t"
               "%INFO/GENE\t%INFO/LEGACY_ID\t%INFO/CDS\t%INFO/AA\n")
        cmd = ["bcftools", "query", "-r", chrom, "-f", fmt, str(cosmic_vcf)]
        n_kept = 0
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True,
                              stderr=subprocess.DEVNULL) as proc:
            assert proc.stdout is not None
            for line in proc.stdout:
                if max_per_chrom and n_kept >= max_per_chrom:
                    proc.terminate()
                    break
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                c, pos, vid, ref, alt, gene, lid, cds, aa = parts[:9]
                if len(ref) > 50:
                    continue
                for a in alt.split(","):
                    if a.startswith("<") or len(a) > 50:
                        continue
                    out.append({
                        "chrom": c,
                        "pos": int(pos),
                        "id": vid if vid and vid != "." else (
                            lid if lid and lid != "." else "."),
                        "ref": ref,
                        "alt": a,
                        "gene": gene if gene and gene != "." else "",
                        "cds": cds if cds and cds != "." else "",
                        "aa": aa if aa and aa != "." else "",
                    })
                n_kept += 1
    return out


def plan_inject_cosmic(
    sites_meta: Sequence[tuple[str, int]],
    cosmic_records: list[dict],
    density: float,
    rng: random.Random,
    reserve_indices: set[int] | None = None,
) -> dict[int, dict[str, Any]]:
    """Decide which sites to replace with COSMIC records.

    Twin of :func:`syntheticgen.clinvar.plan_inject_clinvar` for the
    somatic-overlay path. ``sites_meta`` is a sequence of ``(chrom,
    pos)`` tuples; ``reserve_indices`` excludes sites already claimed
    by an earlier overlay.

    Returns ``{site_index: overlay_record}`` where each
    ``overlay_record`` carries ``pos`` / ``ref`` / ``alts`` and
    optionally ``id`` / ``cosmic_gene`` / ``cosmic_id`` when the
    source record has non-``"."`` ``id`` / non-empty ``gene``. The
    conditional-set semantics mirror :func:`inject_cosmic` exactly so
    the apply path stays trivial.

    RNG consumption order matches :func:`inject_cosmic` exactly.
    """
    if density <= 0 or not sites_meta or not cosmic_records:
        return {}
    chrom_set = {meta[0] for meta in sites_meta}
    pool = [r for r in cosmic_records if r["chrom"] in chrom_set]
    if not pool:
        return {}

    reserve_indices = reserve_indices or set()
    candidate = [i for i in range(len(sites_meta)) if i not in reserve_indices]
    if not candidate:
        return {}

    n_target = max(1, int(round(density * len(sites_meta))))
    n_target = min(n_target, len(candidate), len(pool))

    rng.shuffle(candidate)
    pool_choices = rng.sample(pool, n_target)

    used_keys: set = {(meta[0], meta[1]) for meta in sites_meta}
    plan: dict = {}
    cursor = 0
    for rec in pool_choices:
        key = (rec["chrom"], rec["pos"])
        if key in used_keys:
            continue
        target_i = None
        while cursor < len(candidate):
            i = candidate[cursor]
            cursor += 1
            if sites_meta[i][0] == rec["chrom"]:
                target_i = i
                break
        if target_i is None:
            break
        old_key = sites_meta[target_i]
        used_keys.discard(old_key)
        used_keys.add(key)
        overlay: dict = {
            "pos": rec["pos"],
            "ref": rec["ref"],
            "alts": [rec["alt"]],
        }
        if rec["id"] and rec["id"] != ".":
            overlay["id"] = rec["id"]
            overlay["cosmic_id"] = rec["id"]
        if rec["gene"]:
            overlay["cosmic_gene"] = rec["gene"]
        plan[target_i] = overlay
    return plan


def inject_cosmic(sites: list[dict],
                  cosmic_records: list[dict],
                  density: float,
                  rng: random.Random,
                  reserve_indices: set[int] | None = None) -> int:
    """Replace `density` × len(sites) cohort sites with COSMIC records.

    Mirrors `clinvar.inject_clinvar` but writes COSMIC_GENE / COSMIC_ID
    (and the ID column when present). Returns the number of injections.

    Implementation note: delegates the rng-consuming planning to
    :func:`plan_inject_cosmic` and then applies the resulting
    ``{index: overlay}`` plan in place.
    """
    sites_meta = [(s["chrom"], s["pos"]) for s in sites]
    plan = plan_inject_cosmic(
        sites_meta, cosmic_records, density, rng,
        reserve_indices=reserve_indices,
    )
    for idx, rec in plan.items():
        site = sites[idx]
        site["pos"] = rec["pos"]
        site["ref"] = rec["ref"]
        site["alts"] = rec["alts"]
        if "id" in rec:
            site["id"] = rec["id"]
        if "cosmic_id" in rec:
            site["cosmic_id"] = rec["cosmic_id"]
        if "cosmic_gene" in rec:
            site["cosmic_gene"] = rec["cosmic_gene"]
    sites.sort(key=lambda s: (s["chrom"], s["pos"]))
    return len(plan)
