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
// Resilience: bcftools stats trips on records where GT references an
// allele index that isn't present in ALT — e.g. ALT="." paired with
// GT="1|0", which synthetic_people occasionally emits when the
// rsID-injection step rewrites ALT to "." but leaves the cohort GT
// block intact. Such records are malformed by VCF spec (the GT index
// must reference a real ALT) and contribute zero information to any
// stats: no substitution type, no indel length, no allele balance.
// Dropping them via filter is loss-free.
//
// Symbolic-ALT structural variants (`<DEL>` / `<DUP>` / `<INV>`) are
// kept — bcftools stats counts them under "others" without complaint
// (verified separately) and that count is the only thing the
// bcftools-stats MultiQC panel reports for SVs.
//
// Strategy: try the full-file pass first (cheapest, complete output
// when the input is clean), and on non-zero exit retry against
// `bcftools view -e 'ALT="."'`. The pipeline never aborts on a single
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
        echo "[bcftools_stats] full-file pass exited \$rc; retrying with missing-ALT records filtered out" >&2
        cat stats.err >&2 || true
        bcftools view -e 'ALT="."' "${vcf}" \\
            | bcftools stats -s - > "${name}.bcftools_stats.txt"
    fi
    """
}
