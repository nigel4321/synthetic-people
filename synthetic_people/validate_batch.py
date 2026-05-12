#!/usr/bin/env python3
"""Validation suite (M10) — run against a generated batch.

Walks every `person_*.vcf.gz` under the batch directory, computes the
acceptance-criteria stats from `SYHTHETIC_PROJECT.md` §6, and writes
plots + a Markdown report under `<batch>/validation/`.

Usage::

    .venv/bin/python validate_batch.py path/to/out/

Anything matplotlib-dependent is skipped with a warning if matplotlib
isn't available; the JSON / Markdown artefacts still land.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

# Make sibling package importable when invoked as a plain script.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from syntheticgen.validate import (  # noqa: E402
    aggregate_indel_lengths,
    aggregate_sv_summary,
    af_histogram,
    build_genotype_matrix,
    check_ref_against_fasta,
    cohort_chrom_stats,
    cohort_overlay_density,
    cohort_pca,
    het_hom_ratio,
    ld_decay,
    summarise_vcf,
    titv_from_stats,
)


def _try_plots():
    try:
        from syntheticgen import plots
        return plots
    except ImportError as e:
        print(f"[warn] plotting disabled — matplotlib unavailable: {e}",
              file=sys.stderr)
        return None


def _default_pca_labels(n_samples: int, manifest: dict | None) -> list:
    """Per-sample label used to colour PCA dots.

    If the batch is admixture mode and per-person ancestry fractions are
    in the manifest, label each person by their dominant component
    (so EUR-leaning people show in one colour, SAS in another, etc.).
    Otherwise everyone gets a single 'cohort' label.
    """
    if manifest and manifest.get("mode") == "admixture-uk":
        out = []
        for entry in manifest.get("people", []):
            fracs = entry.get("ancestry_fractions") or {}
            if fracs:
                dom = max(fracs, key=fracs.get)
                out.append(f"{dom}-dominant")
            else:
                out.append("unlabelled")
        if len(out) == n_samples:
            return out
    return ["cohort"] * n_samples


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Validate a synthetic batch (M10).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("batch_dir", type=Path,
                   help="Directory containing person_*.vcf.gz")
    p.add_argument("--out-dir", type=Path, default=None,
                   help="Where validation/ artefacts land. "
                        "Defaults to <batch_dir>/validation/")
    p.add_argument("--af-bins", type=int, default=20,
                   help="Number of AF histogram bins")
    p.add_argument("--pca-components", type=int, default=2,
                   help="Number of PCA components to compute")
    p.add_argument("--ld-pairs-per-bin", type=int, default=5_000,
                   help="Cap on SNP pairs sampled per LD-decay bin")
    p.add_argument("--seed", type=int, default=0,
                   help="RNG seed for LD-decay sub-sampling "
                        "(deterministic plots)")
    p.add_argument("--reference-fasta", type=Path, default=None,
                   help="Reference FASTA (e.g. GRCh38 primary "
                        "assembly) for the Tier 1 REF-matches-"
                        "reference gate. When provided, every "
                        "per-person VCF is checked via `bcftools "
                        "norm --check-ref w` and mismatches are "
                        "reported in the summary + Markdown "
                        "report. Today's synthetic output uses "
                        "fabricated REF (`rng.choice('ACGT')`) so "
                        "this check will fail on every record "
                        "until M12 wires in the real FASTA — "
                        "the gate is in place so a passing run "
                        "after M12 is empirical evidence the "
                        "wiring works. Skips cleanly when omitted.")
    args = p.parse_args(argv)

    if not args.batch_dir.is_dir():
        print(f"batch_dir not a directory: {args.batch_dir}",
              file=sys.stderr)
        return 1

    person_vcfs = sorted(args.batch_dir.glob("person_*.vcf.gz"))
    if not person_vcfs:
        print(f"no person_*.vcf.gz under {args.batch_dir}",
              file=sys.stderr)
        return 1

    out_dir = args.out_dir or args.batch_dir / "validation"
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = args.batch_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) \
        if manifest_path.exists() else None

    print(f"Validating {len(person_vcfs)} VCFs from {args.batch_dir}",
          file=sys.stderr)

    # --- per-sample summaries ---
    samples = []
    for v in person_vcfs:
        s = summarise_vcf(v)
        samples.append(s)
        print(f"  {s.name}: records={s.n_records} snv={s.n_snv} "
              f"indel={s.n_indel} sv={s.n_sv} het={s.n_het} "
              f"hom_alt={s.n_hom_alt} dropouts={s.n_dropout}",
              file=sys.stderr)

    titv = titv_from_stats(samples)
    cohort_het = sum(s.n_het for s in samples)
    cohort_hom_alt = sum(s.n_hom_alt for s in samples)
    cohort_n_records = sum(s.n_records for s in samples)
    cohort_n_singletons = sum(s.singletons for s in samples)
    af_edges, af_counts = af_histogram(samples, n_bins=args.af_bins)
    indel_hist = aggregate_indel_lengths(samples)
    sv_summary = aggregate_sv_summary(samples)
    # Tier 1 validation additions: per-chrom + overlay-density.
    chrom_stats = cohort_chrom_stats(samples)
    overlay_density = cohort_overlay_density(samples)
    overlay_requested = _overlay_requested_from_manifest(manifest)

    # Tier 1: REF-matches-FASTA gate — only when the caller passes
    # --reference-fasta. Today's synthetic REF is fabricated so
    # every record mismatches; the gate exists so a post-M12 run
    # against the real FASTA proves the wiring works.
    ref_check = None
    if args.reference_fasta is not None:
        if not args.reference_fasta.is_file():
            print(
                f"[warn] --reference-fasta not found: "
                f"{args.reference_fasta}; skipping REF check.",
                file=sys.stderr,
            )
        else:
            ref_check = []
            for v in person_vcfs:
                r = check_ref_against_fasta(v, args.reference_fasta)
                ref_check.append(r)
                status = (
                    "PASS" if r["passed"]
                    else ("ERROR" if r["errored"]
                          else f"FAIL ({r['mismatches']} mismatches)")
                )
                print(f"  REF check {v.name}: {status}",
                      file=sys.stderr)

    print(f"Cohort Ti/Tv = {titv:.3f}", file=sys.stderr)
    print(f"Cohort Het/Hom-alt = {cohort_het}/{cohort_hom_alt} = "
          f"{cohort_het / cohort_hom_alt if cohort_hom_alt else float('inf'):.3f}",
          file=sys.stderr)

    # --- LD decay ---
    print("Building genotype matrix for LD decay + PCA...",
          file=sys.stderr)
    matrix, positions, chroms = build_genotype_matrix(person_vcfs)
    print(f"  matrix shape: {matrix.shape}", file=sys.stderr)

    rng = random.Random(args.seed)
    ld = ld_decay(matrix, positions, chroms,
                  pairs_per_bin=args.ld_pairs_per_bin, rng=rng)
    for b in ld:
        print(f"  LD bin {b['low_kb']:>5.1f}–{b['high_kb']:<5.1f} kb: "
              f"n={b['n_pairs']:>5d}  mean_r²={b['mean_r2']:.4f}",
              file=sys.stderr)

    # --- PCA ---
    print("Running cohort PCA...", file=sys.stderr)
    transformed, explained, kept = cohort_pca(
        matrix, n_components=args.pca_components)
    if transformed is not None:
        print(f"  PCA explained variance ratio: {explained}",
              file=sys.stderr)
    else:
        print("  PCA skipped — insufficient variant columns",
              file=sys.stderr)

    pca_labels = _default_pca_labels(len(person_vcfs), manifest)

    # --- artefacts: JSON ---
    summary_json = {
        "batch_dir": str(args.batch_dir),
        "n_samples": len(person_vcfs),
        "n_records_total": cohort_n_records,
        "titv": titv,
        "het": cohort_het,
        "hom_alt": cohort_hom_alt,
        "het_hom_ratio": (cohort_het / cohort_hom_alt
                          if cohort_hom_alt else None),
        "singletons": cohort_n_singletons,
        "af_histogram": {"edges": af_edges, "counts": af_counts},
        "indel_length_distribution": {
            str(k): v for k, v in sorted(indel_hist.items())
        },
        "sv_summary": sv_summary,
        # Tier 1 additions: per-chrom breakouts, overlay density,
        # REF-vs-FASTA gate result.
        "chrom_stats": chrom_stats,
        "overlay_density": {
            "realised": overlay_density,
            "requested": overlay_requested,
        },
        "ref_check": ref_check,
        "ld_decay": ld,
        "pca": {
            "components": args.pca_components,
            "explained_variance_ratio": explained,
            "n_kept_variants": (len(kept) if kept is not None else 0),
            "transformed": (transformed.tolist()
                            if transformed is not None else None),
            "labels": pca_labels,
        },
        "per_sample": [
            {
                "name": s.name, "n_records": s.n_records,
                "n_snv": s.n_snv, "n_indel": s.n_indel, "n_sv": s.n_sv,
                "n_ti": s.n_ti, "n_tv": s.n_tv,
                "n_het": s.n_het, "n_hom_alt": s.n_hom_alt,
                "n_hom_ref": s.n_hom_ref, "n_dropout": s.n_dropout,
                "het_hom_ratio": (het_hom_ratio(s)
                                  if s.n_hom_alt else None),
                "singletons": s.singletons,
                "sv_by_type": dict(s.sv_by_type),
            }
            for s in samples
        ],
    }
    json_path = out_dir / "summary.json"
    json_path.write_text(json.dumps(summary_json, indent=2,
                                    default=str))
    print(f"Summary JSON → {json_path}", file=sys.stderr)

    # --- artefacts: plots ---
    plots = _try_plots()
    plot_paths: dict = {}
    if plots is not None:
        plot_paths["ld_decay"] = plots.plot_ld_decay(
            ld, out_dir / "ld_decay.png")
        plot_paths["af_histogram"] = plots.plot_af_histogram(
            af_edges, af_counts, out_dir / "af_histogram.png")
        plot_paths["indel_lengths"] = plots.plot_indel_lengths(
            indel_hist, out_dir / "indel_lengths.png")
        plot_paths["pca"] = plots.plot_pca(
            transformed, pca_labels, out_dir / "pca.png",
            explained=explained)
        for name, path in plot_paths.items():
            print(f"  plot {name} → {path}", file=sys.stderr)

    # --- artefacts: Markdown report ---
    report = _build_markdown_report(
        args.batch_dir, summary_json, plot_paths, manifest)
    report_path = out_dir / "report.md"
    report_path.write_text(report)
    print(f"Markdown report → {report_path}", file=sys.stderr)
    print("Validation complete.", file=sys.stderr)
    return 0


def _overlay_requested_from_manifest(manifest: dict | None) -> dict:
    """Pull the requested overlay densities from the cli manifest.

    Returns ``{"rsid": float, "clinvar": float, "cosmic": float}``
    with ``None`` for any field the manifest doesn't carry.
    Tier 1 validation pairs these with the realised counts from
    :func:`syntheticgen.validate.cohort_overlay_density` so the
    Markdown report can flag drift between requested and realised.
    """
    out = {"rsid": None, "clinvar": None, "cosmic": None}
    if not manifest:
        return out
    ov = manifest.get("overlays") or {}
    if "rsid_density" in ov:
        out["rsid"] = ov["rsid_density"]
    if "clinvar_inject_density" in ov:
        out["clinvar"] = ov["clinvar_inject_density"]
    if "cosmic_inject_density" in ov:
        out["cosmic"] = ov["cosmic_inject_density"]
    return out


# Tier 1: Ti/Tv tolerance band tightened from the old [1.7, 2.6] (which
# accepted any biologically-plausible noise) to [2.0, 2.2], matching
# the calibrator's actual WGS target and what real human WGS reports
# (typically 2.0–2.1 ± 0.05). A failing run is a useful signal that
# the Ti/Tv calibrator missed its target — either due to a small-n
# noisy cohort or a real regression.
TITV_BAND_LOW = 2.0
TITV_BAND_HIGH = 2.2


def _build_markdown_report(batch_dir: Path, summary: dict,
                           plots: dict, manifest: dict | None) -> str:
    titv = summary["titv"]
    titv_marker = (
        "OK"
        if TITV_BAND_LOW <= titv <= TITV_BAND_HIGH
        else f"outside [{TITV_BAND_LOW}, {TITV_BAND_HIGH}] "
             "(WGS calibrator target)"
    )

    lines = [
        f"# Validation report — {batch_dir.name}",
        "",
        f"- Samples: **{summary['n_samples']}**",
        f"- Records (cohort total): **{summary['n_records_total']}**",
        f"- Ti/Tv: **{titv:.3f}** ({titv_marker})",
        f"- Het / Hom-alt: **{summary['het']} / {summary['hom_alt']}**"
        + (f" = **{summary['het_hom_ratio']:.3f}**"
           if summary["het_hom_ratio"] is not None else ""),
        f"- Singletons: **{summary['singletons']}**",
    ]
    if manifest:
        lines.append(f"- Mode: `{manifest.get('mode', 'unknown')}`")
        lines.append(f"- Build: `{manifest.get('build', 'unknown')}`")
        lines.append(f"- Chromosomes: `{manifest.get('chromosomes', [])}`")
        if manifest.get("mode") == "admixture-uk":
            ap = manifest.get("ancestry_proportions") or {}
            lines.append(
                "- Requested ancestry: "
                + ", ".join(f"{k}={v:.2f}" for k, v in ap.items())
            )
    lines.append("")

    if plots.get("ld_decay"):
        lines += ["## LD decay", "",
                  f"![LD decay]({plots['ld_decay'].name})", ""]
    lines.append("Distance bin | Pairs | Mean r²")
    lines.append("---|---|---")
    for b in summary["ld_decay"]:
        r2_str = (f"{b['mean_r2']:.4f}"
                  if b["mean_r2"] == b["mean_r2"] else "—")
        lines.append(f"{b['low_kb']:.1f}–{b['high_kb']:.1f} kb | "
                     f"{b['n_pairs']} | {r2_str}")
    lines.append("")

    if plots.get("pca"):
        lines += ["## PCA", "",
                  f"![PCA]({plots['pca'].name})", ""]
    pca = summary["pca"]
    if pca.get("explained_variance_ratio"):
        evr = ", ".join(f"PC{i+1}={v:.3f}"
                        for i, v in enumerate(
                            pca["explained_variance_ratio"]))
        lines.append(f"Explained variance ratio: {evr}")
        lines.append("")

    if plots.get("af_histogram"):
        lines += ["## Allele frequency distribution", "",
                  f"![AF]({plots['af_histogram'].name})", ""]

    if plots.get("indel_lengths"):
        lines += ["## Indel length distribution", "",
                  f"![Indels]({plots['indel_lengths'].name})", ""]

    sv = summary.get("sv_summary") or {}
    if sv:
        lines += ["## Structural variants", ""]
        for t, n in sorted(sv.items()):
            lines.append(f"- {t}: **{n}**")
        lines.append("")

    # Tier 1: per-chromosome breakouts — surfaces chrom-specific
    # regressions that cohort-wide aggregates hide.
    chrom_stats = summary.get("chrom_stats") or {}
    if chrom_stats:
        lines += ["## Per-chromosome stats", "",
                  "Chrom | Records | SNV | Indel | SV | Ti | Tv | "
                  "Ti/Tv | Het | Hom-alt",
                  "---|---|---|---|---|---|---|---|---|---"]
        for chrom, b in chrom_stats.items():
            titv_val = b.get("titv", float("nan"))
            titv_str = (
                f"{titv_val:.3f}"
                if isinstance(titv_val, (int, float))
                and titv_val == titv_val
                and titv_val != float("inf")
                else "—"
            )
            lines.append(
                f"{chrom} | {b['n_records']} | {b['n_snv']} | "
                f"{b['n_indel']} | {b['n_sv']} | {b['n_ti']} | "
                f"{b['n_tv']} | {titv_str} | {b['n_het']} | "
                f"{b['n_hom_alt']}"
            )
        lines.append("")

    # Tier 1: realised vs requested overlay density. Catches drift
    # where (e.g.) --rsid-density 0.2 silently produced 0.05.
    od = summary.get("overlay_density") or {}
    realised = od.get("realised") or {}
    requested = od.get("requested") or {}
    if realised.get("n_records"):
        lines += ["## Overlay density (realised vs requested)", "",
                  "Channel | Records | Realised | Requested | Δ",
                  "---|---|---|---|---"]
        for channel in ("rsid", "clinvar", "cosmic"):
            r = realised.get(channel) or {}
            n = r.get("n", 0)
            frac = r.get("fraction", 0.0)
            req = requested.get(channel)
            req_str = f"{req:.4f}" if req is not None else "—"
            delta = (
                f"{frac - req:+.4f}" if req is not None else "—"
            )
            lines.append(
                f"{channel} | {n} | {frac:.4f} | {req_str} | {delta}"
            )
        lines.append("")

    # Tier 1: REF-matches-GRCh38 gate result (when --reference-fasta
    # is provided). Today's synthetic REF is fabricated so this
    # section will list a failure per VCF until M12 lands.
    rc = summary.get("ref_check")
    if rc is not None:
        lines += ["## REF-vs-FASTA check", "",
                  "VCF | Status | Mismatches",
                  "---|---|---"]
        for r in rc:
            name = Path(r["path"]).name
            if r.get("errored"):
                status = "ERROR"
            elif r.get("passed"):
                status = "PASS"
            else:
                status = "FAIL"
            lines.append(
                f"{name} | {status} | {r.get('mismatches', 0)}"
            )
        lines.append("")
        lines.append(
            "*Today's synthetic output uses fabricated REF "
            "(`rng.choice('ACGT')`); this section will report "
            "failures until M12 wires in the real FASTA.*"
        )
        lines.append("")

    lines += ["## Per-sample stats", "",
              "Sample | Records | SNV | Indel | SV | Het | Hom-alt | "
              "Het/Hom | Dropouts",
              "---|---|---|---|---|---|---|---|---"]
    for s in summary["per_sample"]:
        ratio = (f"{s['het_hom_ratio']:.3f}"
                 if s.get("het_hom_ratio") is not None else "—")
        lines.append(
            f"{s['name']} | {s['n_records']} | {s['n_snv']} | "
            f"{s['n_indel']} | {s['n_sv']} | {s['n_het']} | "
            f"{s['n_hom_alt']} | {ratio} | {s['n_dropout']}"
        )
    lines.append("")
    return "\n".join(lines)


if __name__ == "__main__":
    sys.exit(main())
