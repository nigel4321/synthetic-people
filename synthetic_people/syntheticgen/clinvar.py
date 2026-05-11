"""ClinVar fetch + candidate loading + cohort overlay (M7).

Two roles:

* The M1 path: load a small set of pathogenic candidates and use one as
  the per-person "highlighted" variant.
* The M7 path: build a (chrom, pos, ref, alt)-keyed index of ClinVar
  records and either annotate coalescent-produced sites that happen to
  collide with one, or inject a fraction of ClinVar records into the
  cohort so CLNSIG/CLNDN appear at realistic positions. Coalescent sims
  cover positions 1..sim_length while ClinVar sits at real chromosome
  coordinates, so collision-only annotation almost never fires —
  injection is the practical mechanism that makes ClinVar visible in
  the output.

ClinVar's INFO/RS field carries dbSNP rs numbers, so the same cached
file doubles as a rsID source — see `dbsnp.py`.
"""

from __future__ import annotations

import random
import subprocess
import sys
import urllib.request
from pathlib import Path

from .builds import BUILDS


DEFAULT_SIG_FILTER = {
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic",
}

# Default fraction of cohort sites to overwrite with injected ClinVar
# records. Roughly the per-genome ClinVar-known-pathogenic density is
# very low; this knob is exposed mainly so the validation suite can see
# CLNSIG-bearing records in every batch.
DEFAULT_CLINVAR_INJECT_DENSITY = 0.01


def _sanitize_info_value(v: str) -> str:
    """Replace characters that would break VCF INFO parsing.

    ClinVar already uses underscores for spaces; this is defensive only.
    """
    return v.replace(";", ",").replace("=", "_").replace(" ", "_")


def _download(url: str, dest: Path) -> None:
    """Stream-download `url` to `dest.part` then rename to `dest`.

    Writing to a .part file means a partial download never masquerades as
    a completed one on the next run.
    """
    tmp = dest.with_name(dest.name + ".part")
    print(f"  downloading {url}", file=sys.stderr)
    with urllib.request.urlopen(url) as resp:
        total = int(resp.getheader("Content-Length") or 0)
        downloaded = 0
        last_report = 0
        with open(tmp, "wb") as fh:
            while True:
                chunk = resp.read(1 << 20)  # 1 MiB
                if not chunk:
                    break
                fh.write(chunk)
                downloaded += len(chunk)
                if total and downloaded - last_report >= 10 * (1 << 20):
                    pct = downloaded * 100 / total
                    print(f"    {downloaded/1e6:7.1f} / {total/1e6:.1f} MB "
                          f"({pct:.0f}%)", file=sys.stderr)
                    last_report = downloaded
    tmp.rename(dest)


def fetch_clinvar(cache_dir: Path, build: str) -> Path:
    """Ensure the ClinVar VCF + tabix index are cached. Returns VCF path."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    vcf_url = BUILDS[build]["clinvar_url"]
    tbi_url = vcf_url + ".tbi"
    vcf_path = cache_dir / f"clinvar_{build}.vcf.gz"
    tbi_path = vcf_path.with_suffix(vcf_path.suffix + ".tbi")
    if not vcf_path.exists():
        _download(vcf_url, vcf_path)
    if not tbi_path.exists():
        _download(tbi_url, tbi_path)
    return vcf_path


def load_highlighted_candidates(clinvar_vcf: Path,
                                sig_filter: set[str]) -> list[dict]:
    """Stream ClinVar, keeping records whose CLNSIG matches the filter.

    Returns a list of dicts with chrom/pos/id/ref/alt/clnsig/clndn. Skips
    multi-allelic sites and long indels to keep synthetic output simple.
    """
    cmd = [
        "bcftools", "query",
        "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNDN\n",
        str(clinvar_vcf),
    ]
    out: list[dict] = []
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True,
                          stderr=subprocess.DEVNULL) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            chrom, pos, vid, ref, alt, clnsig, clndn = parts[:7]
            if not clnsig or clnsig == ".":
                continue
            # CLNSIG can be pipe- or comma-joined for multi-condition records.
            sigs = {s.strip() for s in clnsig.replace("|", ",").split(",")}
            if not sigs & sig_filter:
                continue
            # Highlighted variants stay single-alt: the "one clinically-
            # highlighted variant per person" concept is inherently one alt.
            #
            # Empty / "." ALT records are also dropped — injecting one
            # while keeping the cohort GT block would land "1|0" against
            # ALT=".", which is invalid VCF and crashes downstream
            # bcftools stats with "Requested allele outside valid range".
            if "," in alt or not alt or alt == "." or alt.startswith("<") or \
                    len(ref) > 50 or len(alt) > 50:
                continue
            out.append({
                "chrom": chrom,
                "pos": int(pos),
                "id": vid if vid else ".",
                "ref": ref,
                "alts": [alt],
                "afs": [None],  # ClinVar doesn't supply AF; filled by writer.
                "clnsig": _sanitize_info_value(clnsig),
                "clndn": _sanitize_info_value(clndn),
            })
    return out


def load_clinvar_index(clinvar_vcf: Path,
                       chromosomes: list[str],
                       sig_filter: set[str] | None = None,
                       max_per_chrom: int | None = None
                       ) -> list[dict]:
    """Stream ClinVar restricted to `chromosomes` into a flat record list.

    Each entry has chrom/pos/ref/alt/clnsig/clndn/id (the ClinVar VCV id)
    and rsid (from INFO/RS, "" if absent). Multi-allelic sites are kept
    by splitting on commas at parse time. Long indels (>50 bp) and
    symbolic ALTs are dropped — same filter as the highlighted-candidate
    loader, so injected records remain "writeable" through the standard
    record path.

    `sig_filter` defaults to the pathogenic-set the highlighted path
    uses; pass an empty set to skip the CLNSIG filter and load all
    annotated records.
    """
    if sig_filter is None:
        sig_filter = DEFAULT_SIG_FILTER
    out: list[dict] = []
    for chrom in chromosomes:
        cmd = [
            "bcftools", "query",
            "-r", chrom,
            "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t"
                  "%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/RS\n",
            str(clinvar_vcf),
        ]
        n_kept = 0
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True,
                              stderr=subprocess.DEVNULL) as proc:
            assert proc.stdout is not None
            for line in proc.stdout:
                if max_per_chrom and n_kept >= max_per_chrom:
                    proc.terminate()
                    break
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                c, pos, vid, ref, alt, clnsig, clndn, rs = parts[:8]
                if not clnsig or clnsig == ".":
                    continue
                if sig_filter:
                    sigs = {s.strip() for s in
                            clnsig.replace("|", ",").split(",")}
                    if not sigs & sig_filter:
                        continue
                if len(ref) > 50:
                    continue
                # Multi-allelic sites: emit one row per alt so each is
                # individually injectable as a biallelic record.
                for a in alt.split(","):
                    # See `load_candidates` for the rationale on the
                    # empty / "." ALT skip.
                    if not a or a == "." or a.startswith("<") or len(a) > 50:
                        continue
                    out.append({
                        "chrom": c,
                        "pos": int(pos),
                        "id": vid if vid else ".",
                        "ref": ref,
                        "alt": a,
                        "clnsig": _sanitize_info_value(clnsig),
                        "clndn": _sanitize_info_value(clndn),
                        "rsid": rs if rs and rs != "." else "",
                    })
                n_kept += 1
    return out


def annotate_clinvar(sites: list[dict],
                     clinvar_records: list[dict]) -> int:
    """Overlay CLNSIG/CLNDN/id onto cohort sites that match a ClinVar entry.

    Matches on (chrom, pos, ref, alt[0]); biallelic-only by design (the
    cohort path is biallelic from M5 onwards). Mutates `sites` in place
    and returns the count of annotated sites.
    """
    index: dict = {}
    for r in clinvar_records:
        index[(r["chrom"], r["pos"], r["ref"], r["alt"])] = r
    n = 0
    for s in sites:
        key = (s["chrom"], s["pos"], s["ref"], s["alts"][0])
        rec = index.get(key)
        if rec is None:
            continue
        s["clnsig"] = rec["clnsig"]
        s["clndn"] = rec["clndn"]
        if s.get("id") in (None, "", ".") and rec.get("id"):
            s["id"] = rec["id"]
        n += 1
    return n


def plan_inject_clinvar(
    sites_meta: list,
    clinvar_records: list[dict],
    density: float,
    rng: random.Random,
) -> dict:
    """Decide which sites to replace with ClinVar records.

    ``sites_meta`` is a sequence of ``(chrom, pos)`` tuples giving the
    light-weight view of the cohort sites — exactly enough information
    for the planner to (a) pick injection indices, (b) avoid position
    collisions with existing sites. The full site dicts (with carriers,
    AF, etc.) are not needed here.

    Returns a ``{site_index: overlay_record}`` dict where each
    ``overlay_record`` carries the fields ``inject_clinvar`` would
    write into the site at that index: ``pos`` / ``ref`` / ``alts``
    (list) / ``id`` / ``clnsig`` / ``clndn``.

    RNG consumption order matches ``inject_clinvar`` exactly:

      1. ``rng.shuffle`` of ``range(len(sites_meta))``.
      2. ``rng.sample`` of ``n_target`` records from the chrom-filtered
         ClinVar pool.

    Phase 5d streaming-cohort follow-up: this function is the
    rng-consuming half of the overlay logic, extracted so it can be
    called from a streaming context (where the full sites list never
    exists in memory) and from the legacy in-place
    :func:`inject_clinvar` (where it still does). Both paths give
    byte-identical output at every fixed seed.
    """
    if density <= 0 or not sites_meta or not clinvar_records:
        return {}
    chrom_set = {meta[0] for meta in sites_meta}
    pool = [r for r in clinvar_records if r["chrom"] in chrom_set]
    if not pool:
        return {}
    n_target = max(1, int(round(density * len(sites_meta))))
    n_target = min(n_target, len(sites_meta), len(pool))

    site_indices = list(range(len(sites_meta)))
    rng.shuffle(site_indices)
    pool_choices = rng.sample(pool, n_target)

    used_keys: set = {(meta[0], meta[1]) for meta in sites_meta}
    plan: dict = {}
    cursor = 0
    for rec in pool_choices:
        key = (rec["chrom"], rec["pos"])
        # Skip ClinVar records whose coordinate is already occupied —
        # keeps positions unique without re-deriving them.
        if key in used_keys:
            continue
        target_i = None
        while cursor < len(site_indices):
            i = site_indices[cursor]
            cursor += 1
            if sites_meta[i][0] == rec["chrom"]:
                target_i = i
                break
        if target_i is None:
            break
        old_key = sites_meta[target_i]
        used_keys.discard(old_key)
        used_keys.add(key)
        plan[target_i] = {
            "pos": rec["pos"],
            "ref": rec["ref"],
            "alts": [rec["alt"]],
            "id": rec["id"],
            "clnsig": rec["clnsig"],
            "clndn": rec["clndn"],
        }
    return plan


def inject_clinvar(sites: list[dict],
                   clinvar_records: list[dict],
                   density: float,
                   rng: random.Random) -> int:
    """Replace `density` × len(sites) cohort sites with ClinVar records.

    The site's GT block (the carrier of LD structure) is preserved; only
    coordinates, REF/ALT, ID and CLNSIG/CLNDN are overwritten. Records
    are picked without replacement from `clinvar_records`, restricted to
    chromosomes that actually appear in `sites`. After injection sites
    remain biallelic and writeable through the standard record path.
    Returns the number of injections performed.

    Sites are sorted by (chrom, pos) on exit, so the cohort site list
    stays monotone for downstream consumers.

    Implementation note: delegates the rng-consuming planning to
    :func:`plan_inject_clinvar` and then applies the resulting
    ``{index: overlay}`` plan in place. Byte-identical at every fixed
    seed to the pre-refactor implementation.
    """
    sites_meta = [(s["chrom"], s["pos"]) for s in sites]
    plan = plan_inject_clinvar(sites_meta, clinvar_records, density, rng)
    for idx, rec in plan.items():
        site = sites[idx]
        site["pos"] = rec["pos"]
        site["ref"] = rec["ref"]
        site["alts"] = rec["alts"]
        site["id"] = rec["id"]
        site["clnsig"] = rec["clnsig"]
        site["clndn"] = rec["clndn"]
    sites.sort(key=lambda s: (s["chrom"], s["pos"]))
    return len(plan)
