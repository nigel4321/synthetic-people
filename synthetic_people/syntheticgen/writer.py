"""Write a single-sample bgzipped + tabix-indexed VCF."""

from __future__ import annotations

import io
import random
import subprocess
from pathlib import Path

from .background import alt_dosages
from .builds import BUILDS, ploidy_for
from .errors import (
    DEFAULT_DROPOUT_RATE,
    DEFAULT_GT_ERROR_RATE,
    maybe_dropout,
    maybe_flip_gt,
    new_error_stats,
)
from .header import build_header
from .quality import (
    DEFAULT_DP_MEAN,
    DEFAULT_DP_SAMPLE_JITTER_SD,
    draw_site_quality,
    gq_from_ad,
    sample_lambda,
)
from .truth import TruthBedWriter, classify_golden


def _contig_sort_key(chrom: str, pos: int,
                     contig_order: dict) -> tuple:
    """Order chromosomes using the reference dict; unknowns sort last."""
    return (contig_order.get(chrom, len(contig_order)), pos)


def write_person_vcf(out_path: Path, person: dict, build: str,
                     rng: random.Random,
                     dp_mean: float = DEFAULT_DP_MEAN,
                     error_rate: float = 0.0,
                     dropout_rate: float = 0.0,
                     stats: dict | None = None,
                     truth_writer: TruthBedWriter | None = None,
                     sex: str | None = None) -> Path:
    """Write a single-sample bgzipped+indexed VCF.

    Each record carries simulated DP/GQ/AD alongside GT. The per-sample
    mean depth is drawn once (Gaussian around `dp_mean`) so different
    samples in a cohort show plausibly different coverage profiles.

    M9 sequencing-error parameters:
      * `error_rate` — per-call probability of a GT flip applied
        *after* AD has been drawn from the truth, so the recomputed GQ
        naturally drops where the call disagrees with the reads.
      * `dropout_rate` — per-call probability of a coverage dropout
        (DP=0, AD all-zero, GQ=0, GT=`./.`).
      * `stats` — optional dict mutated in place with running counts of
        `flipped`, `dropped`, `total_calls`. The caller seeds it with
        `errors.new_error_stats()`.

    M11 truth-set tracking:
      * `truth_writer` — optional `TruthBedWriter` that receives a row
        per golden record (highlighted / ClinVar / COSMIC / SV / rsID)
        and a row per noise event (FLIP / DROPOUT). The caller is
        responsible for `.close()`-ing the writer to flush.

    M13.3 haploid emission (2026-05-17):
      * `sex` — ``"m"`` / ``"f"`` / ``None``. When set, per-record
        ploidy is looked up via :func:`builds.ploidy_for` and the
        emitted GT is collapsed accordingly:

          - **ploidy 0** (chrY in females): the record is dropped
            entirely. The truth-BED writer also skips it — a record
            that didn't reach the VCF can't be a missed call against
            it.
          - **ploidy 1** (chrX non-PAR in males, chrY non-PAR in
            males, MT in everyone): GT collapses to a single-allele
            field (e.g. ``"1"`` rather than ``"1|1"``). The first
            haplotype is kept; this is biologically natural for MT
            (clonal) and chrY (haploid) and a deterministic choice
            for chrX-in-males. AN drops to 1, AD stays a count over
            n_alleles. The DP/GQ noise model is unchanged.
          - **ploidy 2** (autosomes, chrX in females, PAR positions
            in males): today's diploid behaviour.

        When `sex` is None the function preserves pre-M13.3 behaviour
        and emits diploid GTs for every record — back-compat for any
        caller that hasn't been updated to plumb per-person sex
        through yet.
    """
    contigs = BUILDS[build]["contigs"]
    contig_order = {c: i for i, c in enumerate(contigs)}

    records: list[tuple[dict, str, bool]] = []
    records.append((person["highlighted"], person["highlighted"]["gt"], True))
    for bg in person["background"]:
        records.append((bg, bg["gt"], False))
    records.sort(key=lambda r: _contig_sort_key(
        r[0]["chrom"], r[0]["pos"], contig_order))

    str_out = str(out_path)
    if not str_out.endswith(".vcf.gz"):
        raise ValueError("out_path must end in .vcf.gz")

    sample_lam = sample_lambda(dp_mean, DEFAULT_DP_SAMPLE_JITTER_SD, rng)

    # Pipe records straight into `bgzip -c` rather than writing a plain
    # `.vcf` and then shelling out to bgzip. Saves one disk pass per
    # person and removes a fork/exec.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "wb") as out_fh:
        proc = subprocess.Popen(
            ["bgzip", "-c"],
            stdin=subprocess.PIPE,
            stdout=out_fh,
        )
        # Closing `proc.stdin` signals EOF to bgzip. Wrap it in a
        # TextIOWrapper so the existing string-based writes don't have
        # to learn to encode bytes themselves.
        try:
            fh = io.TextIOWrapper(proc.stdin, encoding="utf-8",
                                  write_through=True,
                                  line_buffering=False)
            fh.write(build_header(build, person["sample_id"]))
            for variant, gt, is_hi in records:
                # M13.3: per-record ploidy filtering. When ``sex`` is
                # supplied, drop chrY records for females and collapse
                # haploid records to a single-allele GT.
                if sex is not None:
                    ploidy = ploidy_for(
                        variant["chrom"], sex, build, variant["pos"],
                    )
                    if ploidy == 0:
                        # chrY in females — record doesn't exist
                        # biologically; drop before AD draw + truth
                        # writer to avoid wasted work.
                        continue
                    if ploidy == 1:
                        # Collapse diploid GT to single-allele. Pick
                        # the first haplotype; for MT this is fine
                        # (clonal), for chrY-in-males this is the
                        # paternal haplotype, for chrX non-PAR in
                        # males the choice is arbitrary but
                        # deterministic.
                        first = gt.replace("/", "|").split("|", 1)[0]
                        gt = first
                else:
                    ploidy = 2

                alts = variant["alts"]
                n_alts = len(alts)
                n_alleles = n_alts + 1  # +1 for REF

                # AC is per-alt (Number=A); AN is total allele count
                # (1 for haploid emission, 2 for diploid).
                per_alt_dosage = alt_dosages(gt, n_alts)
                ac_str = ",".join(str(d) for d in per_alt_dosage)
                af_str = ",".join(f"{d/ploidy:.1f}"
                                  for d in per_alt_dosage)
                info_parts = [f"AC={ac_str}", f"AN={ploidy}",
                              f"AF={af_str}"]
                if is_hi:
                    info_parts.append("HIGHLIGHT")
                # Annotation overlays (M7) are picked up regardless of
                # whether this is the highlighted record.
                if variant.get("clnsig") and variant["clnsig"] != ".":
                    info_parts.append(f"CLNSIG={variant['clnsig']}")
                if variant.get("clndn") and variant["clndn"] != ".":
                    info_parts.append(f"CLNDN={variant['clndn']}")
                if variant.get("cosmic_id"):
                    info_parts.append(
                        f"COSMIC_ID={variant['cosmic_id']}")
                if variant.get("cosmic_gene"):
                    info_parts.append(
                        f"COSMIC_GENE={variant['cosmic_gene']}")
                # Structural variants (M8) declare SVTYPE / SVLEN /
                # END / CIPOS. The ALT is a symbolic <DEL>/<DUP>/<INV>
                # already in variant["alts"], so ",".join(alts) writes
                # the right thing — only INFO needs extra fields.
                if variant.get("svtype"):
                    info_parts.append(f"SVTYPE={variant['svtype']}")
                    if variant.get("svlen") is not None:
                        info_parts.append(f"SVLEN={variant['svlen']}")
                    if variant.get("end") is not None:
                        info_parts.append(f"END={variant['end']}")
                    if variant.get("cipos"):
                        lo, hi = variant["cipos"]
                        info_parts.append(f"CIPOS={lo},{hi}")

                # Truth-state DP/AD/GQ first; M9 noise is layered on
                # top so AD reflects what the reads would have looked
                # like under the *true* genotype and any GT flip has
                # to live with low GQ.
                dp, ad, gq = draw_site_quality(gt, n_alleles,
                                               sample_lam, rng)
                called_gt = gt
                event_kind: str | None = None
                if stats is not None:
                    stats["total_calls"] = \
                        stats.get("total_calls", 0) + 1
                if dropout_rate > 0 and maybe_dropout(rng, dropout_rate):
                    called_gt = "./."
                    dp = 0
                    ad = (0,) * n_alleles
                    gq = 0
                    event_kind = "DROPOUT"
                    if stats is not None:
                        stats["dropped"] = stats.get("dropped", 0) + 1
                elif error_rate > 0:
                    new_gt, flipped = maybe_flip_gt(gt, rng, error_rate)
                    if flipped:
                        called_gt = new_gt
                        gq = gq_from_ad(called_gt, ad)
                        event_kind = "FLIP"
                        if stats is not None:
                            stats["flipped"] = \
                                stats.get("flipped", 0) + 1
                if truth_writer is not None:
                    cat = classify_golden(variant, is_hi)
                    if cat is not None:
                        truth_writer.add_golden(variant, cat, gt)
                    if event_kind is not None:
                        truth_writer.add_noise(variant, event_kind,
                                               gt, called_gt)
                ad_str = ",".join(str(x) for x in ad)
                sample_field = f"{called_gt}:{dp}:{gq}:{ad_str}"
                fh.write("\t".join([
                    variant["chrom"], str(variant["pos"]),
                    variant.get("id") or ".",
                    variant["ref"], ",".join(alts),
                    "100", "PASS", ";".join(info_parts),
                    "GT:DP:GQ:AD", sample_field,
                ]) + "\n")
            fh.flush()
        finally:
            # detach() so closing the wrapper doesn't double-close
            # proc.stdin (we still need to close it explicitly so
            # bgzip sees EOF and exits).
            try:
                fh.detach()
            except Exception:
                pass
            if proc.stdin is not None and not proc.stdin.closed:
                proc.stdin.close()
        rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"bgzip exited with code {rc} writing "
                           f"{out_path}")
    subprocess.run(["tabix", "-p", "vcf", "-f", str_out],
                   check=True, capture_output=True)
    return out_path
