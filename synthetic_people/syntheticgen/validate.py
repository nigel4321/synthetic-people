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
import collections
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
    # Tier 2 validation additions (2026-05-13):
    # Per-Mb variant-density bins, keyed by chrom → bin_mb → count.
    # Surfaces gene-density variation in per-region density plots
    # (currently flat under uniform μ; should show structure
    # post-M14 once context-aware μ lands).
    density_bin_counts: dict = field(
        default_factory=lambda: defaultdict(lambda: defaultdict(int))
    )
    # Sampled DP / GQ / AD-ref-fraction values for distribution
    # sanity. Bounded by ``_QUALITY_SAMPLE_CAP`` per VCF so memory
    # stays flat at WGS scale. AD ref-fraction is recorded only at
    # heterozygous calls (where it should hit the empirical 0.475
    # Illumina+BWA-MEM ref-bias target).
    dp_samples: list = field(default_factory=list)
    gq_samples: list = field(default_factory=list)
    ad_het_ref_frac_samples: list = field(default_factory=list)
    # M13.2 sex-chromosome gate counters (2026-05-17). Populated
    # during ``summarise_vcf`` when ``build`` is supplied so the
    # PAR / non-PAR split is correct. ``cohort_sex_chrom_gates``
    # aggregates these against ``manifest.sex[]`` to produce the
    # three pass/fail gates that M13.3 will turn green.
    n_y_records: int = 0            # total chrY records (any pos)
    n_y_non_par_records: int = 0    # subset of chrY outside PAR1/PAR2
    n_y_non_par_het: int = 0        # of n_y_non_par_records, heterozygous
    n_mt_records: int = 0
    n_mt_het: int = 0               # heterozygous MT calls — always wrong


# Cap how many records contribute to the DP/GQ/AD-ratio sample
# arrays per VCF. 50K is enough to render meaningful histograms
# without bloating memory at WGS scale (~5M records per VCF at
# WGS chr22; full WGS ~50M).
_QUALITY_SAMPLE_CAP = 50_000

# Width of per-region density bins. 1 Mb is the natural unit for
# variant-count-per-Mb plots — fine-grained enough to show
# gene-dense vs gene-poor regions but coarse enough that bin
# counts stay meaningful at moderate cohort sizes.
_DENSITY_BIN_MB = 1


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


def summarise_vcf(
    vcf_path: Path, name: str | None = None,
    build: str | None = None,
) -> SampleStats:
    """Walk a VCF and return aggregate counts.

    ``build`` is required for the M13.2 sex-chromosome gate counters
    (``n_y_non_par_*`` and ``n_mt_*``) because PAR boundaries are
    build-specific. When ``build`` is None those counters are still
    populated where build-independent (e.g. ``n_y_records`` /
    ``n_mt_records``) but PAR / non-PAR splits collapse — every chrY
    record gets classified as non-PAR. Callers running the
    sex-chromosome gates must pass ``build``.
    """
    if name is None:
        name = vcf_path.stem.replace(".vcf", "")
    s = SampleStats(name=name)
    quality_samples_taken = 0
    # Local import to avoid a cycle through builds.py if validate.py
    # ever gets imported from inside the cli's import chain.
    from .builds import is_in_par
    for rec in iter_records(vcf_path):
        s.n_records += 1
        chrom_bucket = s.by_chrom[rec.chrom]
        chrom_bucket["n_records"] += 1

        # Tier 2: per-Mb density bin (1-based POS → 0-indexed Mb bin).
        s.density_bin_counts[rec.chrom][
            (rec.pos - 1) // (_DENSITY_BIN_MB * 1_000_000)
        ] += 1

        # Tier 2: sample DP/GQ/AD for distribution sanity, capped at
        # ``_QUALITY_SAMPLE_CAP`` to bound per-VCF memory. First-N
        # sampling (not reservoir) is fine because the cli writes
        # records in genomic order, not in any biased order that
        # would distort the early-N quality distribution.
        if quality_samples_taken < _QUALITY_SAMPLE_CAP:
            s.dp_samples.append(rec.dp)
            s.gq_samples.append(rec.gq)
            ad_total = rec.ad_ref + rec.ad_alt
            if _gt_dosage(rec.gt) == 1 and ad_total > 0:
                # Ref-allele fraction at heterozygous calls — should
                # cluster near the empirical 0.475 Illumina+BWA-MEM
                # ref-bias from ``quality.py``.
                s.ad_het_ref_frac_samples.append(
                    rec.ad_ref / ad_total
                )
            quality_samples_taken += 1
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
        is_het = False
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
                is_het = True
            elif d == 2:
                s.n_hom_alt += 1
                chrom_bucket["n_hom_alt"] += 1

        # M13.2 sex-chromosome gate counters. Normalise the chrom
        # spelling (FASTA-side ``chr22`` vs cli's ``22``) before
        # comparing — the cli emits unprefixed names today but VCFs
        # from external tools may use either convention.
        chrom_norm = rec.chrom[3:] if rec.chrom.startswith("chr") else rec.chrom
        if chrom_norm == "Y":
            s.n_y_records += 1
            in_par = (
                is_in_par("Y", rec.pos, build) if build is not None
                else False  # without a build, treat all Y as non-PAR
            )
            if not in_par:
                s.n_y_non_par_records += 1
                if is_het:
                    s.n_y_non_par_het += 1
        elif chrom_norm == "MT":
            s.n_mt_records += 1
            if is_het:
                s.n_mt_het += 1

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
    is the ti/tv ratio, with ``float('inf')`` when ``tv == 0 and
    ti > 0`` and ``0.0`` when both are zero. Sanitised to ``None``
    for JSON serialisation by ``validate_batch._jsonable``.

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


def cohort_sex_chrom_gates(
    samples: Iterable[SampleStats],
    sample_sex: dict[str, str] | None,
) -> dict:
    """M13.2: aggregate the three sex-chromosome validation gates.

    Three checks, all currently expected to FAIL on M13.1-era output
    (since the simulator still treats every chromosome as diploid).
    M13.3 will turn them green; this function is the empirical gate
    that proves M13.3 actually wired the ploidy lookup through.

    Args:
      samples: per-person ``SampleStats`` from ``summarise_vcf``.
        Each one must have been produced with ``build=...`` passed
        through, otherwise the PAR / non-PAR split is meaningless
        and the function reports ``status="skipped"``.
      sample_sex: ``{sample_name: "m"|"f"}`` from the manifest's
        top-level ``sex`` list (parallel-indexed to ``samples``).
        Pass ``None`` for legacy pre-M13.1 batches; the function
        reports ``status="skipped"`` because the sex labels are
        what determines which gate applies to which sample.

    Returns a dict with three sub-results, each shaped like::

        {
          "status": "pass" | "fail" | "skipped",
          "summary": "<one-line human-readable>",
          ... gate-specific counters ...
        }

    Plus a top-level ``"sample_sex_known": bool`` to make the skip
    reason discoverable in the JSON.
    """
    if sample_sex is None:
        return {
            "sample_sex_known": False,
            "y_het_in_males": {
                "status": "skipped",
                "summary": "no manifest.sex (pre-M13.1 batch)",
            },
            "female_y_absence": {
                "status": "skipped",
                "summary": "no manifest.sex (pre-M13.1 batch)",
            },
            "mt_no_heterozygous": {
                "status": "skipped",
                "summary": "no manifest.sex needed but reporting "
                           "skipped for symmetry with the other "
                           "two gates",
            },
        }

    # Gate 1: Y heterozygosity ≈ 0 in males. Aggregate non-PAR Y
    # heterozygous-call counts across every male sample. Today's
    # diploid-everywhere output produces ~50 % het rate by chance
    # at every Y position; M13.3 emits haploid GT and this drops
    # to 0.
    male_y_non_par = 0
    male_y_non_par_het = 0
    # Gate 2: Female chrY absence. Aggregate Y record counts across
    # every female. Today's output includes Y for both sexes; M13.3
    # will drop Y entirely for females.
    female_y_records = 0
    n_females_seen = 0
    # Gate 3: MT no-heterozygous. Heterozygous MT calls are wrong
    # regardless of sex (MT is haploid and clonally inherited).
    mt_records = 0
    mt_het = 0

    for s in samples:
        sex = sample_sex.get(s.name)
        if sex == "m":
            male_y_non_par += s.n_y_non_par_records
            male_y_non_par_het += s.n_y_non_par_het
        elif sex == "f":
            n_females_seen += 1
            female_y_records += s.n_y_records
        mt_records += s.n_mt_records
        mt_het += s.n_mt_het

    def _gate_y_het() -> dict:
        if male_y_non_par == 0:
            return {
                "status": "skipped",
                "summary": "no non-PAR chrY records in male VCFs "
                           "(chromosome not simulated, or empty)",
                "non_par_records": 0,
                "heterozygous_records": 0,
            }
        het_frac = male_y_non_par_het / male_y_non_par
        # Tolerance: post-M13.3 the rate should be exactly 0. Any
        # non-trivial fraction is a regression. We use 1 % to absorb
        # any future floating-point or edge-case noise without
        # accepting silent regressions.
        passed = male_y_non_par_het == 0
        return {
            "status": "pass" if passed else "fail",
            "summary": (
                f"{male_y_non_par_het} of {male_y_non_par} "
                f"non-PAR chrY records in male VCFs are het "
                f"({het_frac:.1%}); expected 0"
            ),
            "non_par_records": male_y_non_par,
            "heterozygous_records": male_y_non_par_het,
            "het_fraction": het_frac,
        }

    def _gate_female_y() -> dict:
        if n_females_seen == 0:
            return {
                "status": "skipped",
                "summary": "no female samples in cohort",
                "y_records": 0,
            }
        passed = female_y_records == 0
        return {
            "status": "pass" if passed else "fail",
            "summary": (
                f"{female_y_records} chrY records across "
                f"{n_females_seen} female VCFs; expected 0"
            ),
            "y_records": female_y_records,
            "n_females": n_females_seen,
        }

    def _gate_mt() -> dict:
        if mt_records == 0:
            return {
                "status": "skipped",
                "summary": "MT not simulated in this cohort",
                "records": 0,
                "heterozygous": 0,
            }
        het_frac = mt_het / mt_records
        passed = mt_het == 0
        return {
            "status": "pass" if passed else "fail",
            "summary": (
                f"{mt_het} of {mt_records} MT records are het "
                f"({het_frac:.1%}); expected 0 (MT is haploid)"
            ),
            "records": mt_records,
            "heterozygous": mt_het,
            "het_fraction": het_frac,
        }

    return {
        "sample_sex_known": True,
        "y_het_in_males": _gate_y_het(),
        "female_y_absence": _gate_female_y(),
        "mt_no_heterozygous": _gate_mt(),
    }


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

    Iterates ``samples`` once — necessary because the ``Iterable``
    signature lets callers pass a generator. Four separate ``sum()``
    calls would exhaust the generator after the first and zero-out
    every subsequent count.
    """
    total = rs = cln = cos = 0
    for s in samples:
        total += s.n_records
        rs += s.n_with_rs
        cln += s.n_with_clnsig
        cos += s.n_with_cosmic_id
    def _frac(n: int) -> float:
        return (n / total) if total > 0 else 0.0
    return {
        "n_records": total,
        "rsid": {"n": rs, "fraction": _frac(rs)},
        "clinvar": {"n": cln, "fraction": _frac(cln)},
        "cosmic": {"n": cos, "fraction": _frac(cos)},
    }


def cohort_per_region_density(
    samples: Iterable[SampleStats],
) -> dict:
    """Aggregate per-Mb variant density across the cohort (Tier 2 #6).

    Returns ``{chrom: [{"start_mb": int, "end_mb": int, "count":
    int}, ...]}`` sorted by chromosome (canonical 1-22, X, Y, MT
    order) and by bin start within each chrom. Counts are summed
    across samples; for a cohort of n persons each variant is
    counted n times (once per per-person VCF).

    Today's uniform-μ simulation produces flat density (modulo
    statistical noise); post-M14 with context-aware μ the density
    should show real gene-density structure. The validator
    surfaces both the JSON table and a per-chrom line plot.
    """
    aggregated: dict = defaultdict(lambda: defaultdict(int))
    for s in samples:
        for chrom, bins in s.density_bin_counts.items():
            for bin_idx, n in bins.items():
                aggregated[chrom][bin_idx] += n
    out: dict = {}
    for chrom, bins in sorted(
        aggregated.items(), key=_chrom_sort_key,
    ):
        rows = []
        for bin_idx in sorted(bins.keys()):
            start_mb = bin_idx * _DENSITY_BIN_MB
            end_mb = start_mb + _DENSITY_BIN_MB
            rows.append({
                "start_mb": start_mb,
                "end_mb": end_mb,
                "count": bins[bin_idx],
            })
        out[chrom] = rows
    return out


def cohort_quality_metrics(
    samples: Iterable[SampleStats],
) -> dict:
    """Aggregate DP / GQ / AD-ref-fraction sanity stats (Tier 2 #7).

    Returns summary statistics (count, mean, std, percentiles)
    rather than the raw samples — keeps ``summary.json`` compact
    while still exposing whether the empirical distributions
    match ``quality.py``'s claimed model:

    - **DP** ~ Poisson(λ ≈ ``DEFAULT_DP_MEAN``) with per-sample
      Gaussian jitter, so cohort mean should land near the
      configured DP mean.
    - **GQ** is recomputed from AD vs GT; clean hom-ref / hom-alt
      calls hit the cap, hets are slightly lower; cohort median
      should be high.
    - **AD-ref-fraction at hets** is the ref-allele share of
      reads at a heterozygous call. ``quality.py`` uses
      ``HET_ALT_FRAC = 0.475`` (the **alt** share, slightly below
      0.5 due to reference-allele alignment bias). The ref share
      we measure here is therefore ``1 - HET_ALT_FRAC = 0.525``.
      Drift from that band signals either the ``HET_ALT_FRAC``
      constant changed or the AD computation regressed.

    Targets are imported from ``syntheticgen.quality`` rather
    than duplicated as magic numbers — the validator and the
    simulator stay in lockstep automatically when those
    constants are tuned.

    Percentile/median calculations use ``statistics.median`` and
    ``numpy.percentile`` so they're true linear-interpolated
    statistics, not index picks. The old form
    (``values_sorted[int(0.10 * n)]``) misreported p90 as the
    max for small n and gave the upper-middle for even-n median.
    """
    import statistics
    import numpy as np
    from .quality import DEFAULT_DP_MEAN, HET_ALT_FRAC

    dp_all: list = []
    gq_all: list = []
    ad_all: list = []
    for s in samples:
        dp_all.extend(s.dp_samples)
        gq_all.extend(s.gq_samples)
        ad_all.extend(s.ad_het_ref_frac_samples)

    def _summarise(values: list, target: float | None = None) -> dict:
        if not values:
            return {
                "n": 0, "mean": None, "median": None,
                "stdev": None, "p10": None, "p90": None,
                "target": target,
            }
        n = len(values)
        arr = np.asarray(values)
        return {
            "n": n,
            "mean": statistics.fmean(values),
            "median": float(statistics.median(values)),
            "stdev": (statistics.stdev(values)
                      if n >= 2 else 0.0),
            "p10": float(np.percentile(arr, 10)),
            "p90": float(np.percentile(arr, 90)),
            "target": target,
        }

    return {
        # Targets pulled directly from ``syntheticgen/quality.py``
        # constants — drift surfaces the moment those change.
        "dp": _summarise(dp_all, target=DEFAULT_DP_MEAN),
        "gq": _summarise(gq_all, target=None),
        # Ref share = 1 - HET_ALT_FRAC, since we record the ref
        # fraction at hets (not the alt fraction).
        "ad_het_ref_fraction": _summarise(
            ad_all, target=1.0 - HET_ALT_FRAC,
        ),
    }


def cohort_f_statistic(matrix) -> dict:
    """Per-sample inbreeding coefficient F from a dosage matrix
    (Tier 2 #8).

    ``matrix`` is the ``(n_samples, n_variants)`` int dosage
    matrix from :func:`build_genotype_matrix` — 0 = hom-ref,
    1 = het, 2 = hom-alt, -1 = missing. Computes per-variant
    cohort allele frequency p, then per-sample:

    - observed_het = count of dosage==1 sites for this sample
    - expected_het = Σᵥ 2 · p_v · (1 - p_v) across the same sites
    - F = 1 − observed_het / expected_het

    Real outbred cohorts have F ≈ 0 (slightly negative is typical
    under finite-sample sampling variance). Strongly positive F
    (> ~0.05) indicates inbreeding or hidden structure; strongly
    negative F is a hint that the simulator emits too many hets
    relative to the expected HWE distribution.

    Returns ``{"per_sample": [...], "cohort_mean": float,
    "cohort_median": float}``. Per-sample entries include
    ``sample_idx``, ``observed_het``, ``expected_het``, ``f``.
    """
    import numpy as np
    import statistics

    arr = np.asarray(matrix)
    n_samples, n_variants = arr.shape
    if n_samples == 0 or n_variants == 0:
        return {"per_sample": [], "cohort_mean": None,
                "cohort_median": None}

    # Per-variant cohort AF, ignoring missing (-1) entries.
    # arr is int dosage; mask out missing.
    valid = arr >= 0
    n_called = valid.sum(axis=0)
    # Avoid divide-by-zero: variants where no sample has a call.
    safe_n = np.where(n_called > 0, n_called, 1)
    total_dosage = np.where(valid, arr, 0).sum(axis=0)
    p = total_dosage / (2.0 * safe_n)
    # Expected het per variant: 2 p (1 − p), zeroed for
    # zero-call variants.
    expected_het_per_variant = np.where(
        n_called > 0, 2.0 * p * (1.0 - p), 0.0,
    )

    out_per_sample = []
    fs: list = []
    for i in range(n_samples):
        sample_valid = valid[i]
        observed_het = int(((arr[i] == 1) & sample_valid).sum())
        # Sum expected_het only over the variants this sample
        # actually has a call at; otherwise we'd over-count.
        expected_het = float(
            (expected_het_per_variant * sample_valid).sum()
        )
        f = (1.0 - observed_het / expected_het
             if expected_het > 0 else None)
        out_per_sample.append({
            "sample_idx": i,
            "observed_het": observed_het,
            "expected_het": expected_het,
            "f": f,
        })
        if f is not None:
            fs.append(f)
    return {
        "per_sample": out_per_sample,
        "cohort_mean": (statistics.fmean(fs) if fs else None),
        "cohort_median": (
            statistics.median(fs) if fs else None
        ),
    }


def cohort_ancestry_tracts(
    ancestry_bed_paths: list,
) -> dict:
    """Tract-length distribution from admixture-mode ancestry BEDs
    (Tier 2 #9).

    Each per-person ancestry BED carries 5 columns:
    ``chrom\tstart\tend\thap1_pop\thap2_pop`` (one segment per row,
    0-based half-open). A "tract" is a maximal run of consecutive
    rows where the haplotype's ancestry is constant; tract length
    is ``end - start`` in bp.

    Returns ``{"by_population": {pop: {"n": int, "mean_bp": float,
    "median_bp": int}}, "mean_bp_across_pops": float}``.

    Reads the BEDs the M6 admixture path already writes, so this
    check needs no extra cli machinery — just runs when ancestry
    BEDs exist under the batch dir. Useful as a sanity check that
    the realised pulse-time matches the configured one (mean
    tract length L_morgans ≈ 1 / ((1 − α) · T) under a clean
    single-pulse admixture).
    """
    import statistics

    by_pop: dict = defaultdict(list)

    for bed_path in ancestry_bed_paths:
        try:
            text = Path(bed_path).read_text()
        except OSError:
            continue
        # Track tracts for each haplotype independently. State:
        # (chrom, pop) per haplotype; emit tract when either
        # changes.
        for hap_col in (3, 4):  # 0-indexed columns: hap1=3, hap2=4
            current_chrom = None
            current_pop = None
            current_start = 0
            current_end = 0
            for raw_line in text.splitlines():
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue
                chrom, start_s, end_s = parts[0], parts[1], parts[2]
                try:
                    start = int(start_s)
                    end = int(end_s)
                except ValueError:
                    continue
                pop = parts[hap_col]
                if (chrom != current_chrom
                        or pop != current_pop):
                    if current_pop is not None:
                        by_pop[current_pop].append(
                            current_end - current_start,
                        )
                    current_chrom = chrom
                    current_pop = pop
                    current_start = start
                    current_end = end
                else:
                    current_end = end
            # Flush the final tract.
            if current_pop is not None:
                by_pop[current_pop].append(
                    current_end - current_start,
                )

    out = {}
    all_lengths: list = []
    for pop, lengths in sorted(by_pop.items()):
        if not lengths:
            continue
        out[pop] = {
            "n": len(lengths),
            "mean_bp": statistics.fmean(lengths),
            # ``statistics.median`` does the right thing for even-n
            # (interpolates between the two middle values) — the
            # earlier ``lengths_sorted[len // 2]`` picked the upper
            # middle, which made the reported median diverge from
            # ``cohort_f_statistic``'s use of ``statistics.median``.
            "median_bp": int(statistics.median(lengths)),
        }
        all_lengths.extend(lengths)
    return {
        "by_population": out,
        "mean_bp_across_pops": (
            statistics.fmean(all_lengths) if all_lengths else None
        ),
    }


# Canonical chromosome ordering: numeric 1-22 first, then X / Y / MT,
# then everything else alphabetical. Mirrors VCF convention so
# downstream tools see expected order.
_SPECIAL_ORDER = {"X": 23, "Y": 24, "MT": 25, "M": 25}


def _chrom_sort_key(item) -> tuple:
    chrom, _ = item
    # Strip leading "chr" if present so "chr22" and "22" sort together
    # by their semantic rank. ``chrom`` (the raw key) is the tertiary
    # tie-breaker so a mixed-prefix dict (``{"22": ..., "chr22": ...}``)
    # has a deterministic order rather than depending on insertion
    # order alone.
    label = chrom.removeprefix("chr") if chrom.startswith("chr") else chrom
    if label.isdigit():
        return (0, int(label), chrom)
    if label in _SPECIAL_ORDER:
        return (0, _SPECIAL_ORDER[label], chrom)
    return (1, label, chrom)


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
    # ``-Ou`` makes bcftools write the normalised BCF stream to
    # stdout in uncompressed form — on a real WGS VCF that's
    # multiple GB of binary, so we always discard stdout. stderr
    # is streamed line-by-line rather than captured wholesale:
    # ``--check-ref w`` emits one ``REF_MISMATCH`` warning per
    # mismatched record, so on today's fabricated-REF synthetic
    # output the stderr volume can reach hundreds of MB per chrom
    # (one line per record × ~30 bytes × ~1M records). Counting
    # incrementally and retaining only ``_REF_CHECK_STDERR_TAIL``
    # bytes for diagnostics keeps memory bounded.
    try:
        proc = subprocess.Popen(
            ["bcftools", "norm", "--check-ref", "w",
             "-f", str(fasta_path), str(vcf_path), "-Ou"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError as exc:
        return {
            "path": str(vcf_path),
            "passed": False,
            "mismatches": 0,
            "errored": True,
            "stderr_tail": f"<could not run bcftools: {exc}>",
        }

    mismatches = 0
    tail = collections.deque(maxlen=_REF_CHECK_STDERR_TAIL_LINES)
    try:
        # bcftools writes one warning line per mismatch:
        #   "REF_MISMATCH\t<chrom>\t<pos>\t<ref-in-vcf>\t<ref-in-fasta>"
        # ``stderr`` is a binary stream; iterating yields one line at
        # a time so we never materialise the full stderr.
        assert proc.stderr is not None
        for raw_line in proc.stderr:
            line = raw_line.decode(errors="replace")
            if line.startswith("REF_MISMATCH"):
                mismatches += 1
            tail.append(line)
        rc = proc.wait(timeout=600)
    except subprocess.TimeoutExpired as exc:
        proc.kill()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            pass
        return {
            "path": str(vcf_path),
            "passed": False,
            "mismatches": mismatches,
            "errored": True,
            "stderr_tail": f"<timeout: {exc}>",
        }
    stderr_tail = "".join(tail)
    return {
        "path": str(vcf_path),
        "passed": (rc == 0 and mismatches == 0),
        "mismatches": mismatches,
        "errored": rc != 0,
        "stderr_tail": stderr_tail[-1500:],
    }


# Cap how many stderr lines we retain from bcftools for diagnostics.
# bcftools emits one short line per mismatched record, so 200 lines
# is plenty of context if the run fails while bounding the per-VCF
# memory cost to a few KB regardless of cohort size or chrom length.
_REF_CHECK_STDERR_TAIL_LINES = 200


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
