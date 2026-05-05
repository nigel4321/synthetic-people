# variant-scan — Nextflow pipeline

Scan a directory of VCF files for a target variant, producing four markdown reports plus a carrier-level TSV plus a cohort MultiQC HTML:

1. **`qc_report.md`** — per-file validation run before anything else happens: file integrity, tabix index, header parseability, sample/variant presence, plus soft checks that the file looks like human genomic data (recognised reference build, human chromosome names, `FORMAT=GT`, `AF`/`AC`+`AN` present). Hard failures abort the pipeline in strict mode (the default).
2. **`metadata_report.md`** — cohort-level and per-file overview for researchers orienting themselves to an unfamiliar dataset.
3. **`variant_report.md`** — which files contain the target variant within an acceptable allele-frequency range, with per-population breakdowns.
4. **`carriers_report.md`** — how many individuals carry the alt allele, split by heterozygous / homozygous, with integrity checks against cohort AC.
5. **`carriers.tsv`** — one row per (file, sample) where the individual carries the alt allele. Suitable for downstream joins with a sample/population panel.
6. **`multiqc_report.html`** + `multiqc_data/` — interactive cohort report from `bcftools stats` (Ti/Tv per sample, indel-length distribution, SNP/indel/multi-allelic counts, substitution-by-type) plus a custom-content table of `qc_validate.py`'s checks (build, contig sanity, INFO/FORMAT presence, warning/error counts) — suitable for supplementary material of a methods paper or a single-page cohort overview.

Designed for a standalone VM (local executor). Scales to hundreds of VCFs via per-file parallelism.

## Requirements

- Nextflow ≥ 22.10 — `curl -s https://get.nextflow.io | bash`
- `bcftools` ≥ 1.9 on PATH
- Python 3.8+ — `bin/` scripts use stdlib only; the MultiQC stage additionally needs `multiqc>=1.18` importable as a module (`pip install multiqc` or `apt install python3-multiqc`). Invoked as `python3 -m multiqc` so the console-script shim doesn't need to be on PATH.
- Each input VCF must be bgzipped (`.vcf.gz`) with a matching tabix index (`.vcf.gz.tbi`)

## Usage

```bash
nextflow run main.nf \
    --input '/data/vcfs/*.vcf.gz' \
    --outdir results \
    -params-file variants/rs12913832.yaml
```

All variant fields can also be passed as CLI flags (`--variant_name`, `--variant_chrom`, etc.) instead of a params file.

### Parameters

| Parameter | Required | Description |
|---|---|---|
| `--input` | yes | Glob for input VCFs (quote it) |
| `--variant_name` | yes | Label used in report headings |
| `--variant_chrom` | yes | Chromosome; `chr`-prefix optional, pipeline resolves |
| `--variant_pos` | yes | 1-based position |
| `--variant_ref` | yes | Reference allele (strict match) |
| `--variant_alt` | yes | Alternate allele (strict match) |
| `--variant_min_af` | no | Lower AF bound (default 0.0) |
| `--variant_max_af` | no | Upper AF bound (default 1.0) |
| `--outdir` | no | Output directory (default `results`) |
| `--strict_qc` | no | Abort pipeline if any file hard-fails QC (default `true`). Set `false` to run downstream stages even on broken inputs; QC issues are still reported in `qc_report.md`. |

### Adding new variants

Copy `variants/rs12913832.yaml` and edit. One file per variant keeps the parameter history auditable.

## Why match by position + allele rather than rsID

rsIDs are not guaranteed to be present. The 1000 Genomes Phase 3 release, for example, was frozen against a 2013 dbSNP snapshot — many records carry `.` in the ID column even for variants that today have a well-known rsID. Matching by coordinate + REF + ALT is portable across releases and correctly excludes overlapping structural variants (e.g. `<CN0>`, `<CN2>`) at the same position.

## Pipeline layout

```
nextflow_pipeline/
├── main.nf                          # entry workflow
├── nextflow.config                  # executor, resources, params
├── modules/
│   ├── qc_validate.nf               # per-VCF QC + validation (first stage); also emits MultiQC sidecar
│   ├── bcftools_stats.nf            # per-VCF `bcftools stats -s -` for MultiQC's native parser
│   ├── multiqc.nf                   # cohort-level interactive HTML report
│   ├── inspect_vcf.nf               # per-VCF metadata collection
│   ├── scan_variant.nf              # per-VCF variant lookup
│   └── reports.nf                   # markdown aggregation
├── bin/                             # scripts auto-added to PATH by Nextflow
│   ├── qc_validate.py               # --mqc-out emits a MultiQC custom-content sidecar
│   ├── inspect_vcf.py
│   ├── scan_variant.py              # also emits per-file carriers.tsv
│   ├── build_qc_report.py
│   ├── build_metadata_report.py
│   ├── build_variant_report.py
│   └── build_carrier_report.py
├── tests/                           # unittest-based test suite
└── variants/
    └── rs12913832.yaml              # example variant spec
```

## Execution model

Per input VCF, Nextflow runs QC first, then three independent parallel tasks:

- `QC_VALIDATE` — runs before any downstream processing. Hard checks: file readable, tabix index present, header parseable, at least one sample column, at least one variant record. Soft checks (warnings only): reference names a known human build (hg19/GRCh37/hg38/GRCh38), indexed contigs are recognised human chromosomes, `FORMAT=GT` declared, INFO declares `AF` or both `AC`+`AN`, variant positions fit inside plausible human chromosome lengths. In strict mode (default) a hard failure aborts the whole pipeline; soft warnings never abort. Emits both `<name>.qc.json` (consumed by `build_qc_report.py`) and `<name>_qc_mqc.json` (a MultiQC custom-content sidecar).
- `INSPECT_VCF` — contigs, sample count, variant count, reference build heuristic, pipeline/date tags, sample-list hash (to detect cohort mismatches across files).
- `SCAN_VARIANT` — tabix region query, strict REF/ALT match, AF extraction from INFO (or recomputed via `bcftools +fill-tags` if missing), classification into one of seven statuses, **and per-sample genotype extraction** for any file where the variant is present.
- `BCFTOOLS_STATS` — per-VCF `bcftools stats -s -` output. Consumed by MultiQC's built-in bcftools-stats parser to produce per-sample Ti/Tv, indel-length distribution, substitution-by-type, and singleton stats — none of which the markdown reports compute.

Each `SCAN_VARIANT` task emits two files:

- `<name>.variant.json` — cohort-level scan result (status, AF, per-pop AFs, carrier counts)
- `<name>.carriers.tsv` — one row per carrier sample (header-only when the variant is absent from this file, keeping the output channel stable)

Aggregation runs once all per-file fan-outs finish:

- `QC_REPORT` → `qc_report.md` (pass/fail summary + per-file detail table)
- `METADATA_REPORT` → `metadata_report.md`
- `VARIANT_REPORT` → `variant_report.md` (highlights files where the variant is in range)
- `CARRIER_REPORT` → `carriers.tsv` + `carriers_report.md` (concatenates per-file carrier TSVs, totals het/hom counts, verifies `het + 2·hom == cohort AC`)
- `MULTIQC` → `multiqc_report.html` + `multiqc_data/` (cohort-level interactive report combining `bcftools stats` and the QC sidecar)

### Carrier extraction scope

Carrier genotypes are extracted for **every file where the variant is present** — regardless of whether the cohort AF falls within `[min_af, max_af]`. The AF gate only affects *file* classification in `variant_report.md`; individual carriers are reported unfiltered so you can see raw genotype data even in cohorts where cohort-level AF is unusual. Filter downstream by joining `carriers.tsv` against `variant.json` status if needed.

## Variant statuses

Each file's scan emits one of:

| Status | Meaning |
|---|---|
| `not_applicable` | Target chromosome not in this VCF |
| `position_empty` | Chromosome present, but no variant call at the position |
| `absent` | Variant(s) at position but none with the specified REF/ALT |
| `present_in_range` | Variant found and AF ∈ [min_af, max_af] |
| `present_below_threshold` | Variant found but AF < min_af |
| `present_above_threshold` | Variant found but AF > max_af |
| `present_af_unknown` | Variant found but AF could not be determined |

## Profiles

Tune parallelism via `-profile`:

- `standard` — default, `queueSize=16`
- `small_vm` — `queueSize=4`
- `big_vm` — `queueSize=32`

## Output directory

```
results/
├── qc_report.md             # per-file validation: pass/fail, warnings, detail table
├── metadata_report.md       # cohort overview (one section per input file)
├── variant_report.md        # file-level classification + AF per population
├── carriers_report.md       # individual-level summary: het, hom, totals
├── carriers.tsv             # full per-carrier table (file, sample, GT, dosage)
├── multiqc_report.html      # cohort-level interactive HTML (bcftools stats + QC sidecar)
└── multiqc_data/            # MultiQC's parsed-data drop, including multiqc_data.json (machine-readable)
```

Per-file intermediate JSONs and carrier TSVs stay in Nextflow's `work/` directory and can be inspected for debugging.

## Joining carriers to populations

`carriers.tsv` contains the sample IDs but no population labels. For 1000 Genomes data, join against the sample panel:

```bash
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

awk 'NR==FNR {pop[$1]=$3; next}
     FNR==1 {print $0"\tsuper_pop"; next}
     {print $0"\t"pop[$2]}' \
    integrated_call_samples_v3.20130502.ALL.panel \
    results/carriers.tsv > results/carriers_with_pop.tsv
```

This adds a `super_pop` column (AFR, AMR, EAS, EUR, SAS) for downstream ancestry-stratified analysis.

## Integrity check

The `carriers_report.md` summary includes `het + 2·hom`; this should equal the cohort `AC` reported in `variant_report.md`. Mismatches indicate missing/partial genotypes or a bug — the reconciliation is a cheap end-to-end sanity check.
