// Aggregate per-VCF QC + bcftools stats + per-VCF PCA into a single
// interactive HTML report.
//
// MultiQC's built-in bcftools-stats parser picks up `*.bcftools_stats.txt`
// (Ti/Tv per sample, indel-length distribution, SNP / indel counts,
// substitution-by-type). The `*_qc_mqc.json` sidecars from qc_validate.py
// supply our custom QC checks (build, contig sanity, INFO/FORMAT presence,
// warning / error counts) as a custom-content table. The `*_pca_mqc.png`
// images from plot_pca.py land in the report as embedded scatter plots
// — one per input VCF, captioned with the file basename.
//
// Invoked as `python3 -m multiqc` because the `multiqc` console-script
// shim isn't always on PATH (Debian's python3-multiqc install pattern
// installs the package into site-packages without registering an entry
// point). The module form is equivalent and works regardless.
//
// Skip-path PCA outputs (the zero-byte PNGs plot_pca.py writes when the
// input is too small / deps missing) are filtered out before staging
// because MultiQC would otherwise emit a "could not be parsed" warning
// for each empty file.
process MULTIQC {
    publishDir params.outdir, mode: 'copy'

    input:
    path bcftools_stats
    path qc_mqc_jsons
    path pca_mqc_pngs

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    mkdir -p _mqc_in
    for f in ${bcftools_stats} ${qc_mqc_jsons}; do
        cp "\$f" _mqc_in/
    done
    # Skip zero-byte PCA PNGs (plot_pca.py emits these when an input is
    # too small for a meaningful 2-component projection).
    for f in ${pca_mqc_pngs}; do
        if [ -s "\$f" ]; then cp "\$f" _mqc_in/; fi
    done

    # `--filename` strips the .html and uses the basename for the data dir
    # too (e.g. `multiqc_report` → `multiqc_report.html` + `multiqc_report_data/`).
    # We want the canonical `multiqc_report.html` + `multiqc_data/` pair, so
    # rename the data directory after the run rather than letting MultiQC
    # name it from the filename basename.
    python3 -m multiqc --force --outdir . --filename multiqc_report.html _mqc_in/
    if [ -d multiqc_report_data ]; then mv multiqc_report_data multiqc_data; fi
    """
}
