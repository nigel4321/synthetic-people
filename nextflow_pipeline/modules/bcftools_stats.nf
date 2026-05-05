// Per-VCF `bcftools stats` for MultiQC's native bcftools_stats parser.
//
// The output text file is consumed by MultiQC's built-in bcftools-stats
// module, which gives us per-sample Ti/Tv, indel-length distribution,
// SNP / indel / multi-allelic counts, and substitution-by-type bar charts
// in the cohort report — none of which the pipeline computes today.
//
// `-s -` activates per-sample stats blocks; without it bcftools emits
// only the aggregate SN section and MultiQC produces a single-sample
// report regardless of cohort size.
//
// Resilience: bcftools stats trips on records that bcftools considers
// malformed under its variant-type computation — symbolic ALT records
// (`<DEL>` / `<DUP>` / `<INV>`) and records where GT references an
// allele index that isn't present in ALT (e.g. ALT="." with GT="1|0",
// which synthetic_people occasionally emits). When the full-file pass
// fails we retry against a filtered subset that drops both classes,
// producing complete per-sample / SNV / indel stats over the records
// bcftools is happy with. The pipeline never aborts on a single
// problematic record.
process BCFTOOLS_STATS {
    tag "$name"

    input:
    tuple val(name), path(vcf), path(tbi)

    output:
    path "${name}.bcftools_stats.txt", emit: stats

    script:
    """
    set +e
    bcftools stats -s - "${vcf}" > "${name}.bcftools_stats.txt" 2> stats.err
    rc=\$?
    set -e

    if [ \$rc -ne 0 ]; then
        echo "[bcftools_stats] full-file pass exited \$rc; retrying with symbolic-ALT and missing-ALT records filtered out" >&2
        cat stats.err >&2 || true
        bcftools view -e 'ALT="." || ALT~"<"' "${vcf}" \\
            | bcftools stats -s - > "${name}.bcftools_stats.txt"
    fi
    """
}
