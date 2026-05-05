// Per-VCF cohort PCA: scatter PNG + JSON sidecar with per-sample PC
// coordinates and the explained-variance ratio.
//
// Three outputs per process:
//   - ``${name}.pca.png``       — the publishable scatter, lands in
//                                  ``results/pca/`` via the workflow's
//                                  publishDir on the workflow side.
//   - ``${name}.pca.json``      — machine-readable summary (per-sample
//                                  PC1/PC2 + variance explained, or a
//                                  ``{"skipped": "..."}`` marker when
//                                  the input is too small).
//   - ``${name}_pca_mqc.png``   — same image renamed with MultiQC's
//                                  custom-content suffix so the cohort
//                                  HTML report embeds it next to the
//                                  bcftools-stats panels.
//
// The script self-skips on too-few-samples / too-few-variants / missing
// heavy deps, writing a tombstone JSON and a zero-byte PNG so the
// channel always has the declared paths. This way one small input
// doesn't abort the whole pipeline.
process PCA_PLOT {
    tag "$name"

    // Publish only the human-facing pair to results/pca/. The
    // ``_pca_mqc.png`` copy is a consumable for MULTIQC and stays
    // inside the work dir; the glob below excludes it. publishDir
    // patterns are evaluated outside the input scope, so the glob
    // can't reference the process's `name` variable directly.
    publishDir "${params.outdir}/pca", mode: 'copy', pattern: "*.pca.{png,json}"

    input:
    tuple val(name), path(vcf), path(tbi)

    output:
    path "${name}.pca.png",      emit: png
    path "${name}.pca.json",     emit: json
    path "${name}_pca_mqc.png",  emit: mqc

    script:
    """
    plot_pca.py \\
        --vcf           "${vcf}" \\
        --name          "${name}" \\
        --out-png       "${name}.pca.png" \\
        --out-json      "${name}.pca.json" \\
        --min-samples   ${params.pca_min_samples} \\
        --min-variants  ${params.pca_min_variants} \\
        --max-variants  ${params.pca_max_variants}

    cp "${name}.pca.png" "${name}_pca_mqc.png"
    """
}
