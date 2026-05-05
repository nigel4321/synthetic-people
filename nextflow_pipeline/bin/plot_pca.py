#!/usr/bin/env python3
"""Per-VCF cohort PCA plot.

Reads genotypes via ``bcftools query``, builds a samples × variants
dosage matrix (0 / 1 / 2 / NaN-for-missing), mean-imputes missing
calls, drops zero-variance columns, fits a 2-component PCA, and
writes a scatter PNG plus a small JSON sidecar with the per-sample
coordinates and the variance explained by each component.

Designed for the variant-scan pipeline's per-file fan-out: one
invocation per input VCF, so the result captures within-cohort
sample structure for that file. Cross-file (cohort-of-cohorts) PCA
makes assumptions the pipeline can't verify (sample-set overlap,
shared variants) and is intentionally out of scope.

Skip rather than fail when the input is too small to be meaningful
— ``--min-samples`` and ``--min-variants`` gate; a skip writes a
single-line JSON ``{"skipped": "<reason>"}`` and exits 0 so the
Nextflow process is happy. A genuine error (bcftools failure,
missing input, etc.) exits non-zero so the pipeline fails fast.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


# Mapping from a phased / unphased GT string to an alt-dosage [0, 1, 2].
# Anything else (multi-allelic alt indices, ``./.``, partial ``./0``)
# becomes NaN and is mean-imputed downstream. This is the same shape
# ``synthetic_people/syntheticgen/validate.py`` uses, kept simple here
# so the script is dependency-light.
_HOM_REF = {"0|0", "0/0"}
_HET = {"0|1", "1|0", "0/1", "1/0"}
_HOM_ALT = {"1|1", "1/1"}


def _gt_to_dosage(gt: str) -> float:
    if gt in _HOM_REF:
        return 0.0
    if gt in _HET:
        return 1.0
    if gt in _HOM_ALT:
        return 2.0
    # Multi-allelic, missing, partial — treat as missing.
    return float("nan")


def _query_genotypes(vcf: str, max_variants: int):
    """Stream ``bcftools query`` output into ``(samples, dosage_rows)``.

    Returns ``samples`` as a list of strings and ``dosage_rows`` as a list
    of lists of floats. Each row is one variant; columns line up with
    ``samples``. Stops after ``max_variants`` rows are collected — the
    rest of the bcftools output is discarded by closing our read end of
    the pipe, which lets bcftools exit on its next write (SIGPIPE).

    A non-zero bcftools exit code only means real trouble when we read
    every record bcftools intended to emit. If we voluntarily stopped
    early, the exit status is irrelevant — bcftools was killed by our
    pipe close, not by anything wrong with the input.
    """
    samples_proc = subprocess.run(
        ["bcftools", "query", "-l", vcf],
        capture_output=True, text=True, check=False,
    )
    if samples_proc.returncode != 0:
        raise RuntimeError(
            f"bcftools query -l failed (exit {samples_proc.returncode}): "
            f"{samples_proc.stderr.strip()[:200]}"
        )
    samples = [s for s in samples_proc.stdout.splitlines() if s.strip()]
    if not samples:
        return samples, []

    # `bcftools query -f '[%GT\t]\n'` emits one tab-separated row per
    # variant, columns in sample order. Trailing tab is fine — split
    # discards the empty trailing field below.
    proc = subprocess.Popen(
        ["bcftools", "query", "-f", "[%GT\t]\n", vcf],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    assert proc.stdout is not None
    assert proc.stderr is not None
    rows: list[list[float]] = []
    n_cols = len(samples)
    truncated = False
    try:
        for line in proc.stdout:
            line = line.rstrip("\n")
            if not line:
                continue
            cells = line.split("\t")
            # Strip the trailing empty cell from the format string's
            # final \t.
            if cells and cells[-1] == "":
                cells = cells[:-1]
            if len(cells) != n_cols:
                # Defensive: skip malformed rows rather than mis-align
                # the matrix. bcftools should emit one cell per sample
                # but multi-allelic + custom format strings can
                # occasionally surprise us; better to drop the row than
                # to corrupt PCA.
                continue
            rows.append([_gt_to_dosage(g) for g in cells])
            if len(rows) >= max_variants:
                truncated = True
                break
    finally:
        # Closing our read end signals EOF to bcftools's writer. It
        # exits on its next write — usually with SIGPIPE (negative
        # returncode 13 on POSIX) or 141. We only treat that as an
        # error when we did NOT voluntarily stop reading.
        proc.stdout.close()
        stderr_text = proc.stderr.read()
        proc.stderr.close()
        proc.wait()

    if not truncated and proc.returncode not in (0, None):
        raise RuntimeError(
            f"bcftools query failed (exit {proc.returncode}): "
            f"{stderr_text.strip()[:500] or '(no stderr captured)'}"
        )
    return samples, rows


def _write_skip(out_json: str, reason: str) -> None:
    with open(out_json, "w") as fh:
        json.dump({"skipped": reason}, fh, indent=2)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--vcf", required=True)
    p.add_argument("--name", required=True,
                   help="Used in plot title and embedded JSON")
    p.add_argument("--out-png", required=True)
    p.add_argument("--out-json", required=True,
                   help="Per-sample PC coords + variance explained")
    p.add_argument("--max-variants", type=int, default=10000,
                   help="Cap variants used for PCA (random prefix). "
                        "Larger values are slower and rarely change the "
                        "result for typical cohort sizes.")
    p.add_argument("--min-samples", type=int, default=3,
                   help="Skip when fewer than this many samples")
    p.add_argument("--min-variants", type=int, default=10,
                   help="Skip when fewer than this many usable variants "
                        "after zero-variance pruning")
    args = p.parse_args()

    if not Path(args.vcf).is_file():
        print(f"[plot_pca] missing input: {args.vcf}", file=sys.stderr)
        return 1

    # Cheap gates first — these don't need numpy / sklearn / matplotlib
    # so a "too small" skip is honest even on hosts without those deps.
    try:
        samples, rows = _query_genotypes(args.vcf, args.max_variants)
    except RuntimeError as exc:
        # Surface a clean, single-line error message rather than a Python
        # traceback — Nextflow logs are noisy enough already, and the
        # captured bcftools stderr in the message is the actionable bit.
        print(f"[plot_pca] {exc}", file=sys.stderr)
        return 1
    if len(samples) < args.min_samples:
        _write_skip(
            args.out_json,
            f"too few samples ({len(samples)} < {args.min_samples})",
        )
        Path(args.out_png).write_bytes(b"")
        return 0
    if len(rows) < args.min_variants:
        _write_skip(
            args.out_json,
            f"too few variants ({len(rows)} < {args.min_variants})",
        )
        Path(args.out_png).write_bytes(b"")
        return 0

    # Heavy deps are import-guarded so a host without them skips the
    # whole stage cleanly rather than crashing the pipeline. This is
    # the same pattern synthetic_people uses for its validate suite.
    try:
        import numpy as np
        from sklearn.decomposition import PCA
        import matplotlib
        matplotlib.use("Agg")  # headless-safe default
        import matplotlib.pyplot as plt
    except ImportError as exc:
        _write_skip(args.out_json, f"missing dep: {exc.name}")
        # Touch the PNG so the Nextflow output declaration is satisfied.
        Path(args.out_png).write_bytes(b"")
        print(f"[plot_pca] skip: {exc.name} not importable", file=sys.stderr)
        return 0

    # Build the samples × variants matrix (transpose: rows came in as
    # variant-major from bcftools query).
    matrix = np.array(rows, dtype=float).T  # shape (n_samples, n_variants)

    # Mean-impute missing calls column-by-column. A column that is
    # entirely missing collapses to NaN in the mean and gets dropped
    # below alongside zero-variance columns.
    col_means = np.nanmean(matrix, axis=0)
    nan_mask = np.isnan(matrix)
    matrix[nan_mask] = np.take(col_means, np.where(nan_mask)[1])

    # Drop columns that didn't contribute (all-missing or invariant).
    col_var = np.nanvar(matrix, axis=0)
    keep = np.isfinite(col_var) & (col_var > 0)
    matrix = matrix[:, keep]
    if matrix.shape[1] < args.min_variants:
        _write_skip(
            args.out_json,
            f"too few variable variants after pruning "
            f"({matrix.shape[1]} < {args.min_variants})",
        )
        Path(args.out_png).write_bytes(b"")
        return 0

    n_components = min(2, matrix.shape[0] - 1, matrix.shape[1])
    if n_components < 2:
        _write_skip(
            args.out_json,
            f"insufficient rank for 2-component PCA "
            f"(samples={matrix.shape[0]}, variants={matrix.shape[1]})",
        )
        Path(args.out_png).write_bytes(b"")
        return 0

    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(matrix)
    var_pct = [float(v * 100) for v in pca.explained_variance_ratio_]

    # Plot
    fig, ax = plt.subplots(figsize=(6.0, 5.0), dpi=120)
    ax.scatter(coords[:, 0], coords[:, 1], s=24, alpha=0.85,
               edgecolors="black", linewidths=0.4)
    ax.set_xlabel(f"PC1 ({var_pct[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_pct[1]:.1f}%)")
    ax.set_title(f"{args.name} — cohort PCA "
                 f"(n={matrix.shape[0]}, "
                 f"{matrix.shape[1]} variants)")
    ax.axhline(0, color="grey", linewidth=0.5, alpha=0.5)
    ax.axvline(0, color="grey", linewidth=0.5, alpha=0.5)
    fig.tight_layout()
    fig.savefig(args.out_png)
    plt.close(fig)

    summary = {
        "name": args.name,
        "n_samples": int(matrix.shape[0]),
        "n_variants_used": int(matrix.shape[1]),
        "explained_variance_pct": var_pct,
        "samples": [
            {"sample": s, "pc1": float(c[0]), "pc2": float(c[1])}
            for s, c in zip(samples, coords)
        ],
    }
    with open(args.out_json, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[plot_pca] {args.name}: PC1={var_pct[0]:.1f}% "
          f"PC2={var_pct[1]:.1f}% "
          f"({matrix.shape[0]} samples × {matrix.shape[1]} variants)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
