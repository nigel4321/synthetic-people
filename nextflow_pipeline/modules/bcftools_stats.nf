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
process BCFTOOLS_STATS {
    tag "$name"

    input:
    tuple val(name), path(vcf), path(tbi)

    output:
    path "${name}.bcftools_stats.txt", emit: stats

    script:
    """
    bcftools stats -s - "${vcf}" > "${name}.bcftools_stats.txt"
    """
}
