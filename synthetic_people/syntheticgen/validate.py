"""Validation suite analytics (M10).

Pure-stats primitives that summarise a generated cohort against the
acceptance criteria in `SYHTHETIC_PROJECT.md` §6:

* **Variant statistics** — Ti/Tv, Het/Hom per sample, AF distribution,
  singleton fraction, indel length histogram, SV count summary.
* **LD decay curve** — r² between SNP pairs versus physical distance,
  binned on a log scale to match the "biological expectation" plot.
* **PCA** — cohort scatter via scikit-learn, and (optionally) projection
  against a 1000 Genomes reference loaded from a glob of phase-3 VCFs.

Heavy IO sticks to streaming `bcftools query` so very large cohorts
remain tractable. Plotting is gated behind the matplotlib import so a
caller running on a headless box without matplotlib still gets the
stats and structured artefacts written to disk.
"""

from __future__ import annotations

import math
import subprocess
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

from .titv import is_transition, titv_ratio


# Default LD-decay binning: log-spaced across 100 bp – 500 kb, matching
# the canonical "r² vs distance" plot for human WGS. 12 bins is plenty
# to see the monotone decay without over-fragmenting at low counts.
DEFAULT_LD_BINS_KB = (
    (0.1, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 5.0),
    (5.0, 10.0), (10.0, 20.0), (20.0, 50.0), (50.0, 100.0),
    (100.0, 200.0), (200.0, 500.0),
)
# Cap the number of SNP pairs sampled per bin to keep wall time bounded
# even on a 50k-record cohort.
DEFAULT_LD_PAIRS_PER_BIN = 5_000


# ---------------------------------------------------------------------------
# VCF iteration
# ---------------------------------------------------------------------------


@dataclass
class Record:
    """One VCF record — already split per-allele if multi-allelic."""
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: str
    dp: int
    gq: int
    ad_ref: int
    ad_alt: int
    info: dict


def _parse_info(info_str: str) -> dict:
    out: dict = {}
    if info_str in ("", "."):
        return out
    for kv in info_str.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        else:
            out[kv] = True
    return out


def iter_records(vcf_path: Path) -> Iterator[Record]:
    """Stream a single-sample VCF as `Record` objects, splitting
    multi-allelics into one record per ALT allele.

    Uses `bcftools query` to keep memory flat regardless of file size.
    """
    fmt = (
        "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO\\t[%GT\\t%DP\\t%GQ\\t%AD]\\n"
    )
    proc = subprocess.run(
        ["bcftools", "query", "-f", fmt, str(vcf_path)],
        check=True, capture_output=True, text=True,
    )
    for line in proc.stdout.splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        chrom, pos, ref, alt_field, info_str, gt, dp, gq, ad = parts[:9]
        info = _parse_info(info_str)
        try:
            pos_i = int(pos)
        except ValueError:
            continue
        try:
            dp_i = int(dp) if dp not in (".", "") else 0
        except ValueError:
            dp_i = 0
        try:
            gq_i = int(gq) if gq not in (".", "") else 0
        except ValueError:
            gq_i = 0
        ad_parts = ad.split(",") if ad and ad != "." else []
        ad_ref = int(ad_parts[0]) if ad_parts and ad_parts[0].isdigit() \
            else 0
        # Multi-allelics: emit one Record per ALT, with its corresponding
        # AD entry (ad_parts[i+1] for ALT i).
        alts = alt_field.split(",")
        for i, alt in enumerate(alts):
            ad_alt = (
                int(ad_parts[i + 1])
                if len(ad_parts) > i + 1 and ad_parts[i + 1].isdigit()
                else 0
            )
            yield Record(
                chrom=chrom, pos=pos_i, ref=ref, alt=alt, gt=gt,
                dp=dp_i, gq=gq_i, ad_ref=ad_ref, ad_alt=ad_alt,
                info=info,
            )


# ---------------------------------------------------------------------------
# Per-sample stats
# ---------------------------------------------------------------------------


@dataclass
class SampleStats:
    """Aggregate per-VCF (i.e. per-person) summary."""
    name: str
    n_records: int = 0
    n_snv: int = 0
    n_indel: int = 0
    n_sv: int = 0
    n_ti: int = 0
    n_tv: int = 0
    n_het: int = 0
    n_hom_alt: int = 0
    n_hom_ref: int = 0
    n_dropout: int = 0
    af_values: list = field(default_factory=list)
    indel_lengths: list = field(default_factory=list)
    sv_by_type: dict = field(default_factory=lambda: defaultdict(int))
    singletons: int = 0  # records with realised AC == 1 (per VCF)
    # Tier 1 validation additions (2026-05-12):
    # Realised overlay-density counters — number of records carrying
    # each overlay marker. Lets the validator compare what landed in
    # the output against what the manifest requested
    # (catches density-target drift in the overlay pipeline).
    n_with_rs: int = 0          # INFO/RS non-empty (rsID overlay)
    n_with_clnsig: int = 0      # INFO/CLNSIG non-empty (ClinVar overlay)
    n_with_cosmic_id: int = 0   # INFO/COSMIC_ID non-empty (COSMIC overlay)
    # Per-chromosome breakouts — each inner dict tracks the metrics
    # we'd care about for chrom-specific regression detection (e.g.
    # an X-only bug after M13 lands would be invisible to cohort-wide
    # aggregates but obvious in the per-chrom Ti/Tv).
    by_chrom: dict = field(default_factory=lambda: defaultdict(_chrom_bucket))


def _chrom_bucket() -> dict:
    """Default counters for a per-chromosome bucket inside
    ``SampleStats.by_chrom``."""
    return {
        "n_records": 0,
        "n_snv": 0,
        "n_indel": 0,
        "n_sv": 0,
        "n_ti": 0,
        "n_tv": 0,
        "n_het": 0,
        "n_hom_alt": 0,
        "n_hom_ref": 0,
        "n_dropout": 0,
    }


def _classify_record(rec: Record) -> str:
    """Return 'snv', 'indel', or 'sv'."""
    if rec.alt.startswith("<") and rec.alt.endswith(">"):
        return "sv"
    if len(rec.ref) == 1 and len(rec.alt) == 1:
        return "snv"
    return "indel"


def _is_dropout(gt: str) -> bool:
    return "." in gt and gt != "0|0"  # ./. or .|. style


def _gt_dosage(gt: str) -> int:
    """Sum of alt-allele indices in a phased GT. -1 for missing."""
    if "." in gt:
        return -1
    parts = gt.replace("/", "|").split("|")
    if len(parts) != 2:
        return -1
    try:
        a, b = int(parts[0]), int(parts[1])
    except ValueError:
        return -1
    # For multi-allelic dosage we want "this allele's dosage", but the
    # caller iterates per-ALT records so 1 == "this ALT" by construction.
    return (1 if a > 0 else 0) + (1 if b > 0 else 0)


def summarise_vcf(vcf_path: Path, name: str | None = None) -> SampleStats:
    """Walk a VCF and return aggregate counts."""
    if name is None:
        name = vcf_path.stem.replace(".vcf", "")
    s = SampleStats(name=name)
    for rec in iter_records(vcf_path):
        s.n_records += 1
        chrom_bucket = s.by_chrom[rec.chrom]
        chrom_bucket["n_records"] += 1
        kind = _classify_record(rec)
        if kind == "snv":
            s.n_snv += 1
            chrom_bucket["n_snv"] += 1
            if is_transition(rec.ref, rec.alt):
                s.n_ti += 1
                chrom_bucket["n_ti"] += 1
            else:
                s.n_tv += 1
                chrom_bucket["n_tv"] += 1
        elif kind == "indel":
            s.n_indel += 1
            chrom_bucket["n_indel"] += 1
            s.indel_lengths.append(len(rec.alt) - len(rec.ref))
        else:
            s.n_sv += 1
            chrom_bucket["n_sv"] += 1
            svtype = rec.info.get("SVTYPE", "OTHER")
            s.sv_by_type[svtype] += 1

        # Dosage / het-hom / dropout
        if _is_dropout(rec.gt):
            s.n_dropout += 1
            chrom_bucket["n_dropout"] += 1
        else:
            d = _gt_dosage(rec.gt)
            if d == 0:
                s.n_hom_ref += 1
                chrom_bucket["n_hom_ref"] += 1
            elif d == 1:
                s.n_het += 1
                chrom_bucket["n_het"] += 1
            elif d == 2:
                s.n_hom_alt += 1
                chrom_bucket["n_hom_alt"] += 1

        # AF (if present in INFO)
        af_str = rec.info.get("AF")
        if af_str:
            try:
                # Single-sample: AC<=2 / AN=2 → 0.5 or 1.0; we store the
                # per-record AF for histogramming.
                s.af_values.append(float(af_str.split(",")[0]))
            except ValueError:
                pass

        # Singleton: AC == 1 (single-sample VCF).
        ac_str = rec.info.get("AC")
        if ac_str == "1":
            s.singletons += 1

        # Overlay-density counters (Tier 1 validation): a record
        # "carries" an overlay marker if the corresponding INFO field
        # is set and not "." / empty. The cli's overlay paths populate
        # ``INFO/RS`` (dbsnp), ``INFO/CLNSIG`` (clinvar), and
        # ``INFO/COSMIC_ID`` (cosmic). True flag-style fields (no "=")
        # are parsed by ``_parse_info`` as ``True`` — also counts as
        # carrying the marker.
        rs_val = rec.info.get("RS")
        if rs_val is True or (rs_val and rs_val != "."):
            s.n_with_rs += 1
        clnsig_val = rec.info.get("CLNSIG")
        if clnsig_val is True or (clnsig_val and clnsig_val != "."):
            s.n_with_clnsig += 1
        cosmic_val = rec.info.get("COSMIC_ID")
        if cosmic_val is True or (cosmic_val and cosmic_val != "."):
            s.n_with_cosmic_id += 1
    return s


def titv_from_stats(samples: Iterable[SampleStats]) -> float:
    ti = sum(s.n_ti for s in samples)
    tv = sum(s.n_tv for s in samples)
    if tv == 0:
        return float("inf") if ti else 0.0
    return ti / tv


def het_hom_ratio(s: SampleStats) -> float:
    """Het / Hom-alt. Returns inf if no hom-alt; 0 if no het."""
    if s.n_hom_alt == 0:
        return float("inf") if s.n_het else 0.0
    return s.n_het / s.n_hom_alt


def af_histogram(samples: Iterable[SampleStats], n_bins: int = 20):
    """Linear-binned AF histogram over [0, 1]. Returns (edges, counts)."""
    import numpy as np  # local — numpy is a hard dep at M10

    values: list = []
    for s in samples:
        values.extend(s.af_values)
    if not values:
        edges = [i / n_bins for i in range(n_bins + 1)]
        return edges, [0] * n_bins
    arr = np.asarray(values)
    counts, edges = np.histogram(arr, bins=n_bins, range=(0.0, 1.0))
    return list(edges), counts.tolist()


def aggregate_indel_lengths(samples: Iterable[SampleStats]) -> dict:
    """Return {length_bp: count} aggregated across the cohort.

    Lengths are stored as `len(ALT) - len(REF)`, so insertions are
    positive and deletions negative.
    """
    out: dict = defaultdict(int)
    for s in samples:
        for L in s.indel_lengths:
            out[L] += 1
    return dict(out)


def aggregate_sv_summary(samples: Iterable[SampleStats]) -> dict:
    out: dict = defaultdict(int)
    for s in samples:
        for k, v in s.sv_by_type.items():
            out[k] += v
    return dict(out)


def cohort_chrom_stats(samples: Iterable[SampleStats]) -> dict:
    """Aggregate per-chromosome counters across the cohort.

    Returns ``{chrom: {n_records, n_snv, n_indel, n_sv, n_ti, n_tv,
    n_het, n_hom_alt, n_hom_ref, n_dropout, titv}}``. Per-chrom Ti/Tv
    is computed from the per-chrom ti/tv counts (NaN if tv=0).

    Surfaces chrom-specific regressions that cohort-wide aggregates
    hide — e.g. a Y-chromosome ploidy bug after M13 lands would
    barely move the cohort-wide het-rate but would jump out as a
    chrX-only het excess in this table.
    """
    out: dict = defaultdict(_chrom_bucket)
    for s in samples:
        for chrom, bucket in s.by_chrom.items():
            for k, v in bucket.items():
                out[chrom][k] += v
    # Tack a per-chrom titv on each row for direct readability.
    for chrom, bucket in out.items():
        tv = bucket["n_tv"]
        ti = bucket["n_ti"]
        bucket["titv"] = (ti / tv) if tv > 0 else (
            float("inf") if ti else 0.0
        )
    # Sort by canonical chromosome order: 1-22, X, Y, MT, other.
    return dict(sorted(out.items(), key=_chrom_sort_key))


def cohort_overlay_density(samples: Iterable[SampleStats]) -> dict:
    """Compute the realised overlay-density fractions across the cohort.

    Returns a dict with one entry per overlay channel
    (``rsid``, ``clinvar``, ``cosmic``) plus an aggregate
    ``n_records`` count. Each entry is ``{"n": int, "fraction":
    float}`` — the count of records carrying that overlay marker
    and the fraction of all records.

    Used by ``validate_batch.py`` to compare the realised density
    against the manifest's requested density (catches drift in the
    overlay pipeline before it shows up as missing CLNSIG /
    missing rsID counts at use time).
    """
    total = sum(s.n_records for s in samples)
    rs = sum(s.n_with_rs for s in samples)
    cln = sum(s.n_with_clnsig for s in samples)
    cos = sum(s.n_with_cosmic_id for s in samples)
    def _frac(n: int) -> float:
        return (n / total) if total > 0 else 0.0
    return {
        "n_records": total,
        "rsid": {"n": rs, "fraction": _frac(rs)},
        "clinvar": {"n": cln, "fraction": _frac(cln)},
        "cosmic": {"n": cos, "fraction": _frac(cos)},
    }


# Canonical chromosome ordering: numeric 1-22 first, then X / Y / MT,
# then everything else alphabetical. Mirrors VCF convention so
# downstream tools see expected order.
_SPECIAL_ORDER = {"X": 23, "Y": 24, "MT": 25, "M": 25}


def _chrom_sort_key(item) -> tuple:
    chrom, _ = item
    # Strip leading "chr" if present so "chr22" and "22" sort together.
    label = chrom.removeprefix("chr") if chrom.startswith("chr") else chrom
    if label.isdigit():
        return (0, int(label))
    if label in _SPECIAL_ORDER:
        return (0, _SPECIAL_ORDER[label])
    return (1, label)


def check_ref_against_fasta(
    vcf_path: Path, fasta_path: Path,
) -> dict:
    """Run ``bcftools norm --check-ref w`` against a FASTA and report
    pass/fail + mismatch count for one VCF.

    Returns ``{"path": str, "passed": bool, "mismatches": int,
    "errored": bool, "stderr_tail": str}``.

    A "pass" means every record's ``REF`` matches the reference base
    at that position. The current synthetic-people output uses
    fabricated REF (``rng.choice("ACGT")``, see
    ``coalescent.py:440``) so today this check will fail on every
    record — the gate is in place so when M12 wires in the real
    FASTA, a passing run is empirical evidence the wiring works.

    Skip cleanly (returns ``{"errored": True, ...}``) when bcftools
    or the FASTA is unavailable — the caller can surface that
    distinct from a real mismatch.
    """
    try:
        proc = subprocess.run(
            ["bcftools", "norm", "--check-ref", "w",
             "-f", str(fasta_path), str(vcf_path), "-Ou"],
            capture_output=True, check=False, timeout=600,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        return {
            "path": str(vcf_path),
            "passed": False,
            "mismatches": 0,
            "errored": True,
            "stderr_tail": f"<could not run bcftools: {exc}>",
        }
    stderr = proc.stderr.decode(errors="replace")
    # bcftools --check-ref w writes a warning line per mismatch:
    #   "REF_MISMATCH\t<chrom>\t<pos>\t<ref-in-vcf>\t<ref-in-fasta>"
    # Count those; the line prefix is the stable contract.
    mismatches = sum(
        1 for ln in stderr.splitlines() if ln.startswith("REF_MISMATCH")
    )
    return {
        "path": str(vcf_path),
        "passed": (proc.returncode == 0 and mismatches == 0),
        "mismatches": mismatches,
        "errored": proc.returncode != 0,
        "stderr_tail": stderr[-1500:] if stderr else "",
    }


# ---------------------------------------------------------------------------
# Genotype matrix + LD decay
# ---------------------------------------------------------------------------


def build_genotype_matrix(vcf_paths, max_records: int | None = None):
    """Construct a `(n_samples, n_variants)` integer dosage matrix from a
    list of single-sample VCFs.

    Returns `(matrix, positions, chroms)` where `positions[i]` is the
    1-based POS for column `i` and `chroms[i]` is the chromosome.

    Sites are taken as the **union** across all input VCFs; missing
    genotypes are coded `-1`. Multi-allelic sites are reduced to
    "alt-vs-rest" by treating any non-zero allele as alt (consistent
    with the per-allele record splitting in `iter_records`, but here we
    collapse back so each (chrom, pos) is one column).
    """
    import numpy as np

    sites: dict = {}  # (chrom, pos) → column index
    cols_chrom: list = []
    cols_pos: list = []
    rows: list = []
    n_samples = len(list(vcf_paths))

    # We need to resolve the iterator once and keep it; redo a pass by
    # taking the original list:
    paths = list(vcf_paths)
    rows = [dict() for _ in paths]

    for sample_idx, p in enumerate(paths):
        for rec in iter_records(p):
            if max_records is not None and \
                    len(rows[sample_idx]) >= max_records:
                break
            key = (rec.chrom, rec.pos)
            if key not in sites:
                sites[key] = len(cols_chrom)
                cols_chrom.append(rec.chrom)
                cols_pos.append(rec.pos)
            d = _gt_dosage(rec.gt)
            # If multiple ALTs land at the same (chrom, pos), the last
            # one wins — but we apply max() so any alt-supporting call
            # keeps a non-zero dosage.
            col = sites[key]
            existing = rows[sample_idx].get(col, -1)
            if existing < 0 or (d >= 0 and d > existing):
                rows[sample_idx][col] = d

    n_variants = len(cols_chrom)
    mat = np.full((len(paths), n_variants), -1, dtype=np.int8)
    for i, row in enumerate(rows):
        for col, d in row.items():
            mat[i, col] = d
    return mat, cols_pos, cols_chrom


def _r2_pair(g1, g2):
    """r² between two dosage vectors (shape (n_samples,)).

    Missing genotypes (`-1`) are excluded pairwise. Returns NaN if
    either vector has zero variance after filtering.
    """
    import numpy as np

    a = np.asarray(g1, dtype=float)
    b = np.asarray(g2, dtype=float)
    mask = (a >= 0) & (b >= 0)
    if mask.sum() < 4:
        return float("nan")
    a = a[mask]
    b = b[mask]
    va = a.var()
    vb = b.var()
    if va == 0 or vb == 0:
        return float("nan")
    cov = ((a - a.mean()) * (b - b.mean())).mean()
    r = cov / math.sqrt(va * vb)
    return r * r


def ld_decay(matrix, positions, chroms,
             distance_bins_kb=DEFAULT_LD_BINS_KB,
             pairs_per_bin: int = DEFAULT_LD_PAIRS_PER_BIN,
             rng=None):
    """Compute mean r² in each distance bin.

    Returns a list of dicts: `[{"low_kb": .., "high_kb": .., "n_pairs":
    .., "mean_r2": ..}, ...]`. Bins with no pairs report `mean_r2=NaN`.

    The pairing strategy is straightforward: for each bin, walk the
    sorted positions per chromosome and emit every pair whose
    bp-distance falls in the bin, capped at `pairs_per_bin`. A small
    PRNG (`rng=random.Random`) is used to subsample when the bin would
    otherwise overflow — making the curve reproducible.
    """
    import numpy as np

    if rng is None:
        import random
        rng = random.Random(0)

    # Group columns by chromosome
    by_chrom: dict = defaultdict(list)
    for idx, (c, p) in enumerate(zip(chroms, positions)):
        by_chrom[c].append((p, idx))
    for c in by_chrom:
        by_chrom[c].sort()

    out = []
    for lo_kb, hi_kb in distance_bins_kb:
        lo_bp = int(lo_kb * 1000)
        hi_bp = int(hi_kb * 1000)
        pair_indices: list = []
        for c, positions_indices in by_chrom.items():
            n = len(positions_indices)
            for i in range(n):
                p_i, col_i = positions_indices[i]
                # Walk forward until past hi_bp
                for j in range(i + 1, n):
                    p_j, col_j = positions_indices[j]
                    d = p_j - p_i
                    if d < lo_bp:
                        continue
                    if d >= hi_bp:
                        break
                    pair_indices.append((col_i, col_j))

        if len(pair_indices) > pairs_per_bin:
            pair_indices = rng.sample(pair_indices, pairs_per_bin)

        r2s = []
        for col_i, col_j in pair_indices:
            r2 = _r2_pair(matrix[:, col_i], matrix[:, col_j])
            if not math.isnan(r2):
                r2s.append(r2)

        out.append({
            "low_kb": lo_kb,
            "high_kb": hi_kb,
            "n_pairs": len(r2s),
            "mean_r2": float(np.mean(r2s)) if r2s else float("nan"),
        })
    return out


# ---------------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------------


def cohort_pca(matrix, n_components: int = 2):
    """Run PCA on a `(n_samples, n_variants)` dosage matrix.

    Sites with constant or near-constant calls across the cohort are
    dropped before the fit (zero-variance columns make PCA singular).
    Missing genotypes (`-1`) are imputed with the column mean for the
    fit only — typical practice for sample-wise PCA on sparse calls.

    Returns `(transformed, explained_variance_ratio, kept_variant_ids)`.
    """
    import numpy as np
    from sklearn.decomposition import PCA

    arr = np.asarray(matrix, dtype=float)
    # Build a mean-imputed copy for the fit
    masked = arr.copy()
    masked[masked < 0] = np.nan
    # Suppress numpy's "Mean of empty slice" — all-NaN columns are
    # expected (sites only present in samples that were dropped) and
    # we explicitly handle them by falling back to 0 below, after
    # which the column is constant and pruned by the variance filter.
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore", category=RuntimeWarning)
        col_mean = np.nanmean(masked, axis=0)
    col_mean = np.where(np.isnan(col_mean), 0.0, col_mean)
    inds = np.where(np.isnan(masked))
    masked[inds] = np.take(col_mean, inds[1])

    var = masked.var(axis=0)
    keep = var > 1e-9
    masked = masked[:, keep]
    if masked.shape[1] < n_components:
        return None, None, None

    pca = PCA(n_components=n_components)
    transformed = pca.fit_transform(masked)
    return transformed, pca.explained_variance_ratio_.tolist(), \
        keep.nonzero()[0].tolist()
