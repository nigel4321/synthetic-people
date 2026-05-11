# Synthetic People — Tutorial / User Guide

A practical walkthrough of the `synthetic_people` pipeline aimed at
scientists and academic users. The companion `README.md` is the
reference; this file is the recipe book.

Every example assumes you have already followed the **Install** section
of `README.md`: a project venv at `.venv/` with the requirements
installed and `bcftools` / `tabix` / `bgzip` on `PATH`. All commands
are written to be run from the repository root
(`1000genomes/`).

---

## Table of contents

1. [Who this guide is for](#1-who-this-guide-is-for)
2. [The 60-second hello world](#2-the-60-second-hello-world)
3. [Anatomy of an output directory](#3-anatomy-of-an-output-directory)
4. [Inspecting a generated VCF](#4-inspecting-a-generated-vcf)
5. [Research scenarios with full recipes](#5-research-scenarios-with-full-recipes)
   - [5.1 Benchmark an SNV / indel caller](#51-benchmark-an-snv--indel-caller)
   - [5.2 Test an ancestry / admixture inference tool](#52-test-an-ancestry--admixture-inference-tool)
   - [5.3 Validate a clinical-variant interpretation pipeline](#53-validate-a-clinical-variant-interpretation-pipeline)
   - [5.4 Stress-test a structural variant caller](#54-stress-test-a-structural-variant-caller)
   - [5.5 Build a somatic-variant ground truth](#55-build-a-somatic-variant-ground-truth)
   - [5.6 Population genetics — PCA, LD decay, SFS](#56-population-genetics--pca-ld-decay-sfs)
   - [5.7 Compare populations (CEU vs YRI vs CHB)](#57-compare-populations-ceu-vs-yri-vs-chb)
   - [5.8 Multi-chromosome cohorts](#58-multi-chromosome-cohorts)
   - [5.9 Tune the sequencing-error model](#59-tune-the-sequencing-error-model)
   - [5.10 Reproduce a published benchmark exactly](#510-reproduce-a-published-benchmark-exactly)
6. [Working with the truth-set BEDs](#6-working-with-the-truth-set-beds)
7. [Validating a batch with `validate_batch.py`](#7-validating-a-batch-with-validate_batchpy)
8. [Reproducibility and seeding](#8-reproducibility-and-seeding)
9. [Performance and scaling](#9-performance-and-scaling)
10. [Configuration files (optional)](#10-configuration-files-optional)
11. [Troubleshooting](#11-troubleshooting)
12. [Glossary](#12-glossary)

---

## 1. Who this guide is for

You are a researcher who needs **synthetic genomes with known truth**
to:

- benchmark or stress-test a bioinformatics tool (caller, annotator,
  ancestry inference, PRS pipeline);
- prototype an analysis without sharing real patient data;
- teach the structure of a VCF, an admixture pulse, an LD curve, or a
  variant-calling truth set.

This pipeline emits **per-person, single-sample VCFs** with realistic
linkage disequilibrium (`msprime` + `stdpopsim`), realistic site-
frequency spectra, real ClinVar / dbSNP / COSMIC annotations injected
at real chromosome coordinates, optional structural variants, a
calibrated sequencing-error model, and a paired **truth-set BED** for
every person so you can grade downstream calls without ambiguity.

What it is **not**: a read-level simulator. There are no FASTQs, no
BAMs, no aligner artefacts. The pipeline operates in *call space* — it
emits the VCFs a perfect (or perfectly-noisy) caller would produce.

---

## 2. The 60-second hello world

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 5 --seed 42 \
    --chromosomes 22 --chr-length-mb 0.5 \
    --output-dir out_hello
```

That generates 5 single-sample VCFs over a 0.5 Mb prefix of chr22 in
~30 s on a laptop, including:

- realistic LD from a coalescent simulation under
  `OutOfAfrica_3G09` / CEU (the defaults);
- one ClinVar pathogenic variant per person, flagged `HIGHLIGHT`;
- ~20% of cohort sites overlaid with real dbSNP rsIDs;
- a handful of structural variants per person;
- a tiny sequencing-error layer (~0.1% flips, 0.05% dropouts).

The companion smoke script wraps the same kind of run plus the
validation suite into a single command, useful as a "does my install
actually work?" check:

```bash
bash synthetic_people/scripts/smoke.sh
```

It generates a 5-person cohort under `out_smoke/`, runs
`validate_batch.py`, and asserts every advertised deliverable is on
disk. Override `OUT_DIR`, `N_PEOPLE`, `SEED`, or `PYTHON` by env-var.

---

## 3. Anatomy of an output directory

After a run you have:

```
out_hello/
├── person_0001.vcf.gz             ← bgzipped + tabix-indexed VCF
├── person_0001.vcf.gz.tbi
├── person_0002.vcf.gz
├── ...
├── manifest.json                  ← cohort + per-person metadata
├── summary/
│   └── sfs.tsv                    ← cohort allele-count histogram
└── truth/
    ├── person_0001.golden.bed     ← curated "should be called" set
    ├── person_0001.noise.bed      ← every flip / dropout we injected
    └── ...
```

Add `--admixture` and you also get:

```
out_hello/ancestry/
├── person_0001.bed                ← chrom start end h1_pop h2_pop
└── ...
```

Run `validate_batch.py out_hello/` and you also get:

```
out_hello/validation/
├── report.md                      ← human-readable summary
├── summary.json                   ← machine-readable summary
├── ld_decay.png
├── af_histogram.png
├── indel_lengths.png
└── pca.png
```

The `manifest.json` is the canonical index. Read it whenever you need
to know what happened in a run:

```bash
.venv/bin/python -c "import json; m=json.load(open('out_hello/manifest.json'));
print('mode :', m['mode']);
print('seed :', m['seed']);
print('chrs :', m['chromosomes']);
print('n    :', m['n_people']);
print('FDR  :', m['errors']['realised_fdr'])"
```

Per-person entries surface the highlighted ClinVar variant, ancestry
fractions (admixture mode), realised error counts, and the paths to
the truth BEDs.

---

## 4. Inspecting a generated VCF

Every output VCF passes `bcftools view -h` cleanly. Quick sanity
checks any time you want to look inside:

```bash
# Header — confirms reference build, declared INFO/FORMAT tags
bcftools view -h out_hello/person_0001.vcf.gz | head -40

# Variant count
bcftools view -H out_hello/person_0001.vcf.gz | wc -l

# The highlighted ClinVar variant
bcftools view -H -i 'INFO/HIGHLIGHT=1' out_hello/person_0001.vcf.gz

# All ClinVar-tagged sites
bcftools view -H -i 'INFO/CLNSIG!="."' out_hello/person_0001.vcf.gz

# Just the SVs
bcftools view -H -i 'INFO/SVTYPE!="."' out_hello/person_0001.vcf.gz

# Mean DP
bcftools query -f '%INFO/AC\t[%DP\n]' out_hello/person_0001.vcf.gz | \
    awk '{sum+=$2; n++} END{print sum/n}'

# Het / hom-alt / hom-ref / missing breakdown
bcftools query -f '[%GT\n]' out_hello/person_0001.vcf.gz | sort | uniq -c
```

---

## 5. Research scenarios with full recipes

Each recipe shows the generation command, what it produces, and a
worked example of using the output. Adapt the cohort size (`--n`)
and chromosome span (`--chr-length-mb`) to your compute budget.

### 5.1 Benchmark an SNV / indel caller

**Goal:** measure the precision and recall of a variant caller against
ground truth without depending on Genome-in-a-Bottle or other curated
real datasets.

**Generate:**

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 30 --seed 100 \
    --chromosomes 22 --chr-length-mb 5.0 \
    --error-rate 0.005 \
    --dropout-rate 0.001 \
    --output-dir bench_caller
```

This raises the per-call flip rate to 0.5% and dropout to 0.1% so the
truth-vs-call divergence is large enough to grade meaningfully. The
emitted VCF *is* the noisy caller output; the BEDs tell you exactly
which calls are real and which are model noise.

**Grade your caller's output `your_calls/person_0001.vcf.gz` against
the truth:**

```bash
# Convert your caller's output to a BED of called sites
bcftools query -f '%CHROM\t%POS0\t%END\n' your_calls/person_0001.vcf.gz \
    > called.bed

# True positives = called ∩ golden
bedtools intersect -u -a called.bed \
    -b bench_caller/truth/person_0001.golden.bed | wc -l

# False negatives = golden ∖ called
bedtools intersect -v \
    -a bench_caller/truth/person_0001.golden.bed -b called.bed | wc -l

# False positives = called ∖ (golden ∪ noise) — i.e. things the model
# never emitted, so any hit here is genuinely the caller's fault.
cat bench_caller/truth/person_0001.golden.bed \
    bench_caller/truth/person_0001.noise.bed | sort -k1,1 -k2,2n \
    > all_emitted.bed
bedtools intersect -v -a called.bed -b all_emitted.bed | wc -l
```

The `noise.bed` is the killer feature: a "false positive" that lands
on a noise row is not the caller's fault — it's a known-bad call we
seeded. Subtract those before reporting precision.

### 5.2 Test an ancestry / admixture inference tool

**Goal:** validate a tool like RFMix, Loter, or your own ancestry
deconvolution against per-haplotype, per-segment truth.

**Generate a UK-style admixed cohort:**

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 50 --seed 7 \
    --chromosomes 22 --chr-length-mb 10.0 \
    --admixture \
    --eur-frac 0.60 --sas-frac 0.25 --afr-frac 0.15 \
    --output-dir admix_test
```

This runs a `demes`-defined three-source demography (EUR + SAS + AFR)
joining at a single admixture pulse 20 generations ago into a UK deme.
For each person it writes:

- `person_NNNN.vcf.gz` — the genomic data;
- `ancestry/person_NNNN.bed` — per-segment haplotype ancestry truth.

**Truth BED columns:**

```
chrom  start    end       hap1_pop  hap2_pop
22     0        1234567   EUR       AFR
22     1234567  5000000   EUR       EUR
...
```

Per-segment truth lets you compute haplotype-level (not just
person-level) accuracy.

**Inspect realised cohort proportions vs requested:**

```bash
.venv/bin/python -c "
import json
m = json.load(open('admix_test/manifest.json'))
print('Requested:', m['ancestry_proportions'])
for p in m['people'][:5]:
    print(f\"  {p['sample_id']}: {p['ancestry_fractions']}\")
"
```

The cohort-mean realised ancestry will sit within finite-cohort
sampling noise of the requested mix; finer-grained agreement requires
a larger `--n`.

### 5.3 Validate a clinical-variant interpretation pipeline

**Goal:** confirm your interpretation pipeline correctly flags the
ACMG-style pathogenic variants buried inside a normal-looking cohort
background.

**Generate:**

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 20 --seed 42 \
    --chromosomes 22 --chr-length-mb 5.0 \
    --clinvar-sig "Pathogenic,Likely_pathogenic,Pathogenic/Likely_pathogenic" \
    --clinvar-inject-density 0.02 \
    --rsid-density 0.30 \
    --output-dir clinical_test
```

`--clinvar-inject-density 0.02` overlays ~2% of cohort sites with real
ClinVar pathogenic records at their true chromosome coordinates;
`--rsid-density 0.30` ensures ~30% of sites carry real dbSNP rsIDs so
the cohort looks realistic to an annotator.

**Confirm what's in there:**

```bash
# Pathogenic variants flagged in person_0001
bcftools view -H -i 'INFO/CLNSIG~"athogenic"' \
    clinical_test/person_0001.vcf.gz \
    | head

# All highlighted variants across the cohort
for f in clinical_test/person_*.vcf.gz; do
    bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/CLNSIG\t%INFO/CLNDN\n' \
        -i 'INFO/HIGHLIGHT=1' "$f"
done
```

The truth BED tells you exactly which variants your pipeline must
return as pathogenic:

```bash
grep "flag=HIGHLIGHTED\|flag=CLINVAR" \
    clinical_test/truth/person_0001.golden.bed \
    | head
```

Tune `--clinvar-sig` to broaden or narrow the included clinical
significance terms (e.g. add `risk_factor` or `Affects`).

### 5.4 Stress-test a structural variant caller

**Goal:** measure recall of an SV caller across a length spectrum.

**Generate a cohort heavy on SVs:**

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 10 --seed 11 \
    --chromosomes 22 --chr-length-mb 20.0 \
    --svs-per-person 30 \
    --sv-length-min 50 \
    --sv-length-max 100000 \
    --output-dir sv_bench
```

Each person picks up 30 SVs (mix ~50/30/20 % DEL / DUP / INV) drawn
log-uniformly across 50 bp – 100 kb. Anchor REFs are placeholder bases
(see "Known gaps" in the README), so use the `SVTYPE` / `SVLEN` /
`END` INFO tags rather than the REF column.

**Extract the SV truth subset:**

```bash
grep "flag=SV" sv_bench/truth/person_0001.golden.bed > person_0001.sv_truth.bed

# How many SVs by type?
awk -F'\t' '{
    match($4, /svtype=([^;]+)/, arr); count[arr[1]]++
} END { for (k in count) print k, count[k] }' \
    person_0001.sv_truth.bed
```

**Recall by length bucket:**

```bash
# Bin truth SVs by length, then ask which buckets your caller recovers.
awk -F'\t' '{
    match($4, /svlen=(-?[0-9]+)/, arr); l = arr[1]; if (l < 0) l = -l;
    if (l < 100)        bucket = "1.short_50_100";
    else if (l < 1000)  bucket = "2.med_100_1000";
    else if (l < 10000) bucket = "3.large_1k_10k";
    else                bucket = "4.huge_10k+";
    print bucket "\t" $0;
}' person_0001.sv_truth.bed
```

### 5.5 Build a somatic-variant ground truth

**Goal:** test a somatic caller (e.g. Mutect2, Strelka2) against
COSMIC-grounded truth.

**Prerequisite:** download a COSMIC VCF from the official portal (free
academic registration required). The pipeline never auto-fetches
COSMIC because of the registration gate.

**Generate:**

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 15 --seed 99 \
    --chromosomes 22 --chr-length-mb 5.0 \
    --somatic --cosmic-vcf /path/to/CosmicCodingMuts.vcf.gz \
    --cosmic-inject-density 0.01 \
    --output-dir somatic_test
```

Each cohort site has a 1% chance of being overwritten with a real
COSMIC record; the per-person VCFs end up carrying `COSMIC_ID` and
`COSMIC_GENE` in the INFO field. Note the pipeline emits a single
"sample" per VCF — for a tumour-normal pair workflow you would
generate two cohorts with shared seeding and treat one as the
matched normal.

**Extract COSMIC truth subset:**

```bash
grep "flag=COSMIC" somatic_test/truth/person_0001.golden.bed | head
```

### 5.6 Population genetics — PCA, LD decay, SFS

**Goal:** generate a plausible population for teaching or testing
demographic-inference workflows. Realistic LD and an empirical-shape
SFS are the point.

**A modest 100-sample CEU cohort:**

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 100 --seed 2024 \
    --chromosomes 22 --chr-length-mb 10.0 \
    --demo-model OutOfAfrica_3G09 --population CEU \
    --output-dir popgen_ceu100
```

**Validate the population-genetic shape:**

```bash
.venv/bin/python synthetic_people/validate_batch.py popgen_ceu100/
```

Look at `popgen_ceu100/validation/`:

- `ld_decay.png` — should show monotone r² decay from ~0.5 at
  100–500 bp down toward zero by ~100 kb;
- `af_histogram.png` — heavy on rare alleles (singleton-dominated);
- `pca.png` — a single tight cluster (no real population structure
  in a single-pop draw);
- `report.md` — Ti/Tv (target ≈ 2.1), Het/Hom, dropouts, indel
  breakdown.

The cohort SFS is also persisted to `popgen_ceu100/summary/sfs.tsv`
(columns: `ac`, `n_sites`) for direct plotting:

```bash
# Quick text-mode SFS plot
sort -n popgen_ceu100/summary/sfs.tsv | head -30
```

### 5.7 Compare populations (CEU vs YRI vs CHB)

**Goal:** generate three single-population cohorts under matched
parameters to demonstrate population-specific LD or AF patterns.

```bash
for POP in CEU YRI CHB; do
    .venv/bin/python synthetic_people/generate_people.py \
        --n 50 --seed 42 \
        --chromosomes 22 --chr-length-mb 5.0 \
        --demo-model OutOfAfrica_3G09 --population "$POP" \
        --output-dir "popgen_${POP}"
    .venv/bin/python synthetic_people/validate_batch.py "popgen_${POP}/"
done
```

Each population uses the same demographic model but different
sampling deme. Expect AFR (`YRI`) to show shorter LD blocks and a
higher proportion of common variation than `CEU` / `CHB`, mirroring
the empirical literature.

### 5.8 Multi-chromosome cohorts

Pass a comma-separated chromosome list, a numeric range, or a mix of
both. Each chromosome runs as a separate `msprime` simulation, then the
cohort sites are concatenated in genome order before per-person VCFs are
written:

```bash
# explicit list
.venv/bin/python synthetic_people/generate_people.py \
    --n 30 --seed 17 \
    --chromosomes 19,20,21,22 --chr-length-mb 5.0 \
    --output-dir multichr

# numeric range — all autosomes
.venv/bin/python synthetic_people/generate_people.py \
    --n 30 --seed 17 \
    --chromosomes 1-22 --chr-length-mb 5.0 \
    --output-dir all-autosomes

# mix ranges with singletons (autosomes + X)
.venv/bin/python synthetic_people/generate_people.py \
    --n 30 --seed 17 \
    --chromosomes 1-22,X --chr-length-mb 5.0 \
    --output-dir autosomes-plus-x
```

Ranges are inclusive and must be numeric (`1-22` works, `X-Y` does not).
Duplicates collapse silently and the resolved order follows the order
you wrote.

Memory cost scales with `n_people × n_chromosomes × chr_length_mb`.
For full-chromosome runs use `--chr-length-mb 0` (slower, RAM-hungrier
on chr1–chr5).

### 5.9 Tune the sequencing-error model

The lightweight model is configurable end-to-end. The defaults
(0.001 flip, 0.0005 dropout) target a realistic 0.15% per-call FDR.

**Disable noise entirely (clean truth-only output):**

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 10 --seed 1 \
    --chromosomes 22 --chr-length-mb 1.0 \
    --error-rate 0 --dropout-rate 0 \
    --output-dir clean
```

The `truth/person_NNNN.noise.bed` files will be empty (zero-byte) and
every called GT matches the truth GT exactly.

**Crank noise up to stress-test a downstream filter:**

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 10 --seed 1 \
    --chromosomes 22 --chr-length-mb 1.0 \
    --error-rate 0.05 --dropout-rate 0.02 \
    --output-dir noisy
```

5% flip + 2% dropout; you'll see the realised FDR climb to ~7% in the
manifest's `errors.realised_fdr` field. Flipped calls land low-GQ
because GQ is recomputed from AD (which still reflects the truth) —
a quality filter targeting `GQ >= 20` should reject most of them.

### 5.10 Reproduce a published benchmark exactly

Same seed + same flags ⇒ byte-identical output across runs (and
across machines, modulo `msprime` / `numpy` / `stdpopsim` versions).
The `manifest.json` records the seed plus every parameter that
matters, so a paper can ship a single command and a manifest file.

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 100 --seed 2024 \
    --chromosomes 22 --chr-length-mb 10.0 \
    --output-dir benchmark_v1

# Persist the exact reproduction command alongside the run
cat > benchmark_v1/REPRODUCE.sh <<'EOF'
.venv/bin/python synthetic_people/generate_people.py \
    --n 100 --seed 2024 \
    --chromosomes 22 --chr-length-mb 10.0 \
    --output-dir benchmark_v1
EOF
```

Pin your `requirements.txt` (or a conda lock) to make the
reproduction byte-stable across collaborators.

---

## 6. Working with the truth-set BEDs

Both BEDs are 0-based half-open BED4. The 4th column is a
semicolon-separated `key=value` payload:

```text
22  17104711  17104712  flag=RSID;id=rs340582;ref=C;alt=G;gt=1|1
22  29673445  29673446  flag=CLINVAR;id=190387;ref=C;alt=T;gt=0|1;clnsig=Pathogenic;clndn=Stromme_syndrome
22  167001    173953    flag=SV;id=.;ref=G;alt=<DEL>;gt=1|1;svtype=DEL;svlen=-6951
```

Each golden row carries exactly one priority-ordered tag:

```
HIGHLIGHTED   the per-person highlighted ClinVar variant
CLINVAR       any other ClinVar-annotated cohort row
COSMIC        COSMIC-injected somatic record
SV            structural variant (DEL / DUP / INV)
RSID          dbSNP-known SNP/indel
```

**Filter golden BED by category:**

```bash
grep "flag=CLINVAR" truth/person_0001.golden.bed
grep -E "flag=(HIGHLIGHTED|CLINVAR)" truth/person_0001.golden.bed
```

**Pull truth GTs as a TSV:**

```bash
awk -F'\t' '{
    match($4, /id=([^;]+)/, id); match($4, /gt=([^;]+)/, gt);
    print $1"\t"$2"\t"$3"\t"id[1]"\t"gt[1]
}' truth/person_0001.golden.bed | head
```

**Noise BED records both the truth and the call:**

```text
22  583912  583913  flag=FLIP;id=.;ref=A;alt=G;truth_gt=0|0;called_gt=0|1
22  712451  712452  flag=DROPOUT;id=.;ref=C;alt=T;truth_gt=1|1;called_gt=./.
```

So a "false positive" that intersects a `flag=FLIP` row is a known
seeded error, not a caller bug. Filter accordingly when you compute
caller precision.

---

## 7. Validating a batch with `validate_batch.py`

Drop-in companion CLI:

```bash
.venv/bin/python synthetic_people/validate_batch.py path/to/out/
```

Walks every `person_*.vcf.gz`, computes per-sample + cohort stats,
builds an `(n_samples × n_variants)` dosage matrix, and emits
`<batch>/validation/{report.md, summary.json, ld_decay.png,
af_histogram.png, indel_lengths.png, pca.png}`. The `report.md` is the
primary consumer artefact:

```bash
glow path/to/out/validation/report.md         # if you have glow
.venv/bin/python -m markdown path/to/out/validation/report.md
less path/to/out/validation/report.md
```

Useful flags:

| Flag | Effect |
|---|---|
| `--af-bins 30` | Number of bins in the AF histogram |
| `--pca-components 2` | Number of PCs to retain |
| `--ld-pairs-per-bin 5000` | SNP-pair sub-sample size for LD decay (lower = faster) |
| `--seed 0` | RNG seed for the LD pair sub-sampling |
| `--out-dir custom_dir/` | Override the default `<batch>/validation/` path |

If `matplotlib` is not installed the PNGs are skipped with a warning;
the JSON / Markdown artefacts still land.

---

## 8. Reproducibility and seeding

- **`--seed N`** — same flags + same seed = byte-identical output
  across runs on the same hardware / dependency versions.
- Omit `--seed` and you get a fresh draw each invocation: different
  sample IDs, different highlighted variants, different SVs.
- The seed feeds *every* downstream RNG — coalescent simulation,
  ClinVar / dbSNP / COSMIC injection, SV draws, sequencing-error
  perturbations. There is no second hidden seed.
- The `manifest.json` records the seed verbatim, so any cohort can be
  re-derived from the run log alone.

For papers / reproducibility studies, also pin:

- The Python package versions (`requirements.txt` or conda lock).
- The htslib version (`bcftools --version`).
- The ClinVar VCF cached under `synthetic_people/cache/` — the
  upstream file is updated weekly.

---

## 9. Performance and scaling

Rough ballpark for a single-population coalescent run (laptop,
8-core, 16 GB RAM):

| Cohort | Chrs | Length / chr | Time | Peak RAM |
|---|---|---|---|---|
| `--n 5`   | 22 only | 0.5 Mb | ~30 s   | ~0.5 GB |
| `--n 50`  | 22 only | 5 Mb   | ~3 min  | ~2 GB   |
| `--n 200` | 22 only | 10 Mb  | ~10 min | ~5 GB   |
| `--n 50`  | 19,20,21,22 | 5 Mb each | ~12 min | ~3 GB |

Scaling rules of thumb:

- `msprime` runtime is roughly `O(n_haps² × length × recombination)`,
  but for sub-100-sample cohorts the demographic-model setup
  dominates.
- Memory grows with `n_people × n_sites_per_chrom`. Cut
  `--chr-length-mb` aggressively for prototyping.
- Admixture mode is ~1.5× slower than single-population coalescent
  (extra demography + per-haplotype tree walks).
- Validation (`validate_batch.py`) builds a full dosage matrix in RAM;
  drop `--ld-pairs-per-bin` to 1000 if you're tight.
- **Beyond `--n 30 000` the in-RAM sites-list path stops working**
  (CPython refcount-COW divergence across forked workers; see
  `PERFORMANCE_PLAN.md` §5d for the diagnosis). The `--cohort-mode
  arrow` path detailed in §9.1.1 below replaces the in-RAM hand-off
  with a memory-mapped Arrow IPC scratch file per chromosome, makes
  peak parent RSS bounded as `n` grows, and is what unblocks the
  100 000+ scale target. `--cohort-mode auto` (the default) picks it
  for `--n >= 100 000` automatically; `pip install pyarrow` first.

### 9.1 Worker processes (Phase 1, updated by Phase 5e)

The cohort BCF write and the per-person VCF fan-out run in parallel
process pools on Linux. Control the fan-out with `--workers`:

| `--workers` | Behaviour |
|---|---|
| `0` (default) | Auto = `os.cpu_count()`. |
| `1` | Serial. Useful for profiling and for environments where parallelism is unwanted. |
| `N > 1` | Up to `N` worker processes parallelising (a) the per-chromosome cohort BCF write across sample slices and (b) the per-person VCF fan-out across persons. |

> **Phase 5e change (2026-05):** `--workers` no longer parallelises
> msprime simulation itself. The pre-5e behaviour ("one msprime
> worker per chromosome") multiplied tree-sequence RAM by the worker
> count and OOM-killed workstation-class hosts at `n=3000+`. After
> 5e, simulation runs serially across chromosomes in the parent
> process — only the cohort BCF write step (per chromosome,
> sample-slice parallel) and the per-person fan-out are scaled by
> `--workers`. The per-chrom cohort wall time goes from ~1270 s
> (serial) to ~450 s with `--workers 8` on a 32 GB host at n=3000;
> the simulation phase (~308 s/chrom) remains serial as the
> wall-time floor.

```bash
# Force serial for an apples-to-apples profiling baseline.
.venv/bin/python synthetic_people/generate_people.py \
    --n 30 --seed 17 --chromosomes 22 --chr-length-mb 5.0 \
    --workers 1

# Cap at 4 workers on a constrained host (default would have used all cores).
.venv/bin/python synthetic_people/generate_people.py \
    --n 30 --seed 17 --chromosomes 19,20,21,22 --chr-length-mb 5.0 \
    --workers 4
```

Determinism: output is identical for any `--workers` value at a fixed
`--seed`. Output **does** differ from pre-Phase-1 runs at the same
seed, however — Phase 1 changed how the master rng is consumed (one
`randint` per chromosome and per person up front). If you need to
reproduce a pre-Phase-1 run exactly, check out the pre-Phase-1 commit.

Non-Linux hosts (macOS / Windows) silently fall back to `--workers=1`
because the parallel path uses fork-based multiprocessing. A warning
is printed when you ask for more workers than the host can provide.

The writer also pipes records straight into `bgzip -c` instead of
writing a plain `.vcf` first; this is transparent to the user but
shaves one disk pass per person.

### 9.1.1 Cohort intermediate — `--cohort-mode` (Phase 5d.1)

> **Install first:** the Arrow path needs `pyarrow`, which isn't in
> `requirements.txt`. Run `.venv/bin/pip install pyarrow` once if
> you'll be running with `--n >= 100 000` or passing
> `--cohort-mode arrow` explicitly. Without it, `--cohort-mode auto`
> falls back to `sites_list` and prints a one-line warning.

Between simulation and the BCF write, the cohort lives in some
intermediate form. `--cohort-mode {sites_list, arrow, arrow-streaming, auto}`
selects which:

| Mode | What happens | When to use |
|---|---|---|
| `sites_list` (Phase B) | Parent holds the per-chromosome sites list in RAM as sparse-carriers Python dicts. Workers fork-COW share the list. | small `--n` (≤ ~30 000) — fastest path, no pyarrow dep |
| `arrow` (Phase 5d.1) | Parent materialises the per-chromosome sites list and streams it into `out/cohort/.arrow/cohort.chr<N>.arrow` (Apache Arrow IPC, batch-size 256 sites). Workers `pa.memory_map` and consume their sample slice as a zero-copy numpy view. The Arrow scratch is created + deleted per chromosome. | medium `--n` (≥ 100 000) where workers benefit from mmap, but parent's sites list still fits in RAM |
| `arrow-streaming` (Phase 5d.1 PR 3 — option 3) | Parent **streams directly from the msprime tree sequence into the Arrow file**, never materialising the full sites list. Overlay injection (ClinVar / rsID / COSMIC) and `annotate_clinvar` are applied inline. Parent peak drops from ~17 GB to a few hundred MB at WGS scale. | large workloads where the materialised parent peak overflows host RAM (e.g. WGS at `--n 3000+` on a 32 GB host) |
| `auto` (default) | First checks whether predicted materialised parent peak exceeds 50 % of host RAM (via psutil); if yes picks `arrow-streaming`. Otherwise picks `arrow` for `--n ≥ 100 000` else `sites_list`. If pyarrow is missing, falls back to `sites_list` with a one-line warning. | nearly always |

Why three paths: at small n the sites-list path is fastest. At
medium n the mmap-share win matters for workers but parent fits in
RAM, so the materialised+arrow path is fine. At very large `n` or
`chr_length_mb`, the materialised parent itself overflows host RAM
(memprof28 measured 17 GB at n=3000 × 70 Mb; WGS chr1 at n=3000
extrapolates to ~60 GB). The `arrow-streaming` path was added so
parent never holds the full sites list — required for WGS-scale
runs on commodity hosts, and the foundation for the n=1M target.

**Output is byte-identical across all three modes at a fixed
`--seed`.** The choice is purely about resource use; downstream
BCF + per-person VCF outputs are independent of which path
generated them.

Disk and IO at scale (only relevant when the run picks `arrow`):

- **Per-chromosome scratch:** roughly `2 × n_samples × 5000 × chr_length_mb` bytes for the Arrow file plus the partial BCFs. At n=1M × 70 Mb that's ~80 GB Arrow + ~80 GB partial BCFs in flight per chromosome. The Arrow file is deleted on success after the merge; partials too. The cli runs a pre-flight check at startup and aborts cleanly if free disk under `out/cohort/` can't hold one chromosome's scratch.
- **Final cohort BCFs accumulate:** ~80 GB × 22 chromosomes ≈ 1.7 TB at n=1M. Reserve **2–3 TB total** before kicking off a million-sample run.
- **Disk type matters more than CPU.** Per-chrom data moves through disk twice (write Arrow, read mmap) plus tree-sequence dump and partial-BCF writes — ~400 GB of I/O per chromosome at n=1M. NVMe at ~3 GB/s does that in ~140 s; a spinning disk at ~150 MB/s takes ~45 minutes per chromosome, ~16 hours just for I/O across 22 chromosomes. **Run cohort jobs at scale on NVMe, not HDD.**

To force a path explicitly:

```sh
$ python generate_people.py \
    --n 200000 --seed 42 --chromosomes 22 \
    --cohort-mode arrow \
    --workers 8
```

Output is byte-identical between `sites_list` and `arrow` at the same seed (verified by `tests/test_cohort_arrow_cli.py::CohortModeArrowParityTest`).

### 9.2 Overlay prefetch (Phase 2)

The ClinVar / dbSNP / COSMIC loaders are bcftools-driven and bound on
subprocess I/O. They are submitted to a small thread pool *before* the
coalescent simulation runs, so the loader work overlaps with msprime
instead of running serially after it. There's no flag — prefetch is
always on for the non-legacy path.

You'll see a line like

```
  prefetching overlay loaders in background: clinvar_index, rsid_pool
```

right before the simulation starts, and `Awaiting <name>...` lines
inside the overlay block as each future is resolved. If a loader
finished while the simulation was running, the `Awaiting` line returns
immediately.

### 9.3 Output shape — `--mode` (Phase 5a)

`--mode` controls *what* the run writes to disk. Default is
`per-person`, the layout every existing user is used to:

```bash
generate_people.py --n 100 --seed 42 --chromosomes 22 \
    --chr-length-mb 5 --mode per-person   # default; flag optional
```

Three values:

- `per-person` (default) — emit one bgzipped+tabixed VCF per person,
  identical to today's behaviour. No cohort BCF is written.
- `cohort` — stream the cohort chromosome-by-chromosome straight
  onto disk as `out/cohort/cohort.chr<N>.bcf` + CSI sidecars. Skips
  the per-person fan-out — derive per-person VCFs later via
  `bcftools view -s SAMPLE_ID out/cohort/cohort.chr<N>.bcf` (or
  `bcftools concat` across chromosomes).

  Phase 5b2 made `--mode per-person` and `both` use the same
  streamed pipeline: per-person VCFs are derived from the cohort
  BCFs via `bcftools view -s` rather than from an in-memory
  `cohort_sites` list. The cohort BCFs land alongside the per-person
  VCFs as intermediates and are listed in the manifest's
  `cohort_bcfs[]` field; users who only want per-person can `rm -rf
  out/cohort` after the run.
- `both` — emit both deliverables in the same run.

Cohort mode is the scaling-friendly format for large `--n`. Two
things compose to keep peak RAM down:

1. The streaming cohort flow (Phase 5b1) holds at most one
   chromosome's worth of cohort sites in memory at a time — the
   chromosome is simulated, overlaid, written to its own
   `cohort.chr<N>.bcf`, and freed before the next chromosome is
   simulated.
2. The per-person fan-out is skipped entirely — no fork-pool of N
   worker processes each touching the cohort_sites payload and
   dirtying COW pages.

The same host that OOM-kills on `--n 30 --chromosomes 1-22
--chr-length-mb 70` under `--mode per-person` finishes cleanly under
`--mode cohort`. Quick comparison at `n=500 × 3 chroms × 5 Mb`: peak
RSS drops from ~1.9 GB (per-person) to ~0.85 GB (cohort streamed),
wall time from ~2:27 to ~0:19.

Example: derive a single sample's VCF from cohort BCFs after the
fact:

```bash
# Single chromosome
bcftools view -s HG12345 -Oz \
    out/cohort/cohort.chr22.bcf > out/person_HG12345.chr22.vcf.gz

# All chromosomes — concat per-chrom slices
bcftools concat -Oz -o out/person_HG12345.vcf.gz \
    <(bcftools view -s HG12345 -Ou out/cohort/cohort.chr1.bcf) \
    ...
tabix -p vcf out/person_HG12345.vcf.gz
```

Phase 5b2 added a resume contract for these long runs: a
`cohort.meta.json` file alongside the per-chrom BCFs records the
run's params + sample IDs + per-person seeds + per-chromosome
overlay seeds + the list of chromosomes whose BCF has finished
writing. If a run is killed by OOM / SIGINT / node failure, just
re-run with the same flags and the same `--output-dir`; completed
chromosomes get reused, only the missing ones are re-simulated.
Mismatched params surface a clear error rather than silently
re-using incompatible state. Pass `--no-resume` to wipe everything
and start from scratch.

The manifest's `shape` field records which `--mode` produced the
run. `samples[]` is always present at the top level (any mode) so
listing the cohort doesn't need a per-mode code path; `people[]`
appears in `per-person` and `both`; `cohort_bcfs[]` appears in
`cohort` and `both` (a list — singleton in 5a, populated with
per-chromosome paths once Phase 5b lands).

### 9.4 Progress logging on long runs

The cohort BCF write loop and the per-person fan-out both emit
throttled (~20 s cadence) heartbeat lines so a multi-hour run has
visible progress without flooding stderr. You'll see things like:

```
  cohort BCF: 1,234,567/2,285,875 sites (24,500/s)
  person VCFs: 1,234/100,000 written (5/s, elapsed 240s, eta 19000s)
```

Small cohorts that finish in under twenty seconds skip the
intermediate logs and just print the final summary.

### 9.5 Chunked simulation — `--chr-chunk-mb` (Phase 5f)

At cohort sizes around `n=3000` and full-chromosome lengths,
msprime's working memory during a single chromosome's simulation
can exceed available RAM on workstation-class hosts (8-16+ GB at
`n=3000 × chr1 × 70 Mb × OutOfAfrica_3G09`). Phase 5f splits each
chromosome into independent sub-chunks so per-chunk working memory
fits the host:

```bash
generate_people.py --n 3000 --seed 45 --chromosomes 1-22 \
    --chr-length-mb 70 --output-dir ~/out
# auto-picks --chr-chunk-mb based on free RAM and --workers
```

The default (`--chr-chunk-mb 0`) auto-picks a chunk size at run
start using `psutil.virtual_memory().available` and the
configured worker count, aiming for the per-chunk working set
(estimated from `n × chunk_size × demo_model_factor`) to stay
under ~50% of free RAM. The chosen value is logged so the user
sees what was picked:

```
  --chr-chunk-mb auto-picked 8.74 Mb (available RAM 16.0 GB,
  --workers 4, n=3000)
```

To pin a specific chunk size — useful for reproducibility across
heterogeneous hosts, or to force smaller chunks for safety —
pass an explicit value:

```bash
generate_people.py … --chr-chunk-mb 5
```

**Cross-chunk LD caveat.** Chunks simulate independently, so
linkage disequilibrium decays sharply at chunk boundaries.
Each chunk simulates ~5-10% past its declared end (boundary
smoothing — variants past the chunk end are dropped at write
time, but their presence makes the central region's coalescent
context less truncated), but this isn't true cross-chunk LD
recovery. Effects on common analyses:

| Analysis | Impact under chunking |
|---|---|
| Ti/Tv | unaffected |
| Allele frequency spectrum | unaffected |
| Per-person genotype lists | unaffected (each chunk's haplotypes stay consistent across the cohort) |
| ClinVar / dbSNP / COSMIC overlay placement | unaffected (overlays operate within chunks) |
| Short-range LD (≤ chunk size) | preserved within chunks |
| Long-range LD (> chunk size) | NOT realistic — analyses requiring chr-scale haplotype block structure should not use chunked mode |

If your analysis depends on chr-scale LD, use
`--chr-chunk-mb N` with `N ≥ chr_length_mb` (or pass `0` on a
host with enough RAM that auto-pick keeps the full chromosome
in one chunk).

### 9.6 Diagnosing memory pressure — `--profile-memory`

When a run OOMs and the auto-picked chunk size still doesn't bring
peak RAM under the host's limit, the next step is to look at the
RSS curve over time and figure out which phase is the offender.
Pass `--profile-memory PATH` to spawn a background thread that
samples this process's RSS once per second to a TSV at `PATH`,
plus drops labelled checkpoint rows at every key transition
(ClinVar fetch, overlay loads, per-chromosome simulation, BCF
write, per-person fan-out start/end). Every write is flushed +
fsynced so an OOM kill preserves the trace up to the kernel reap.

```bash
generate_people.py --n 3000 --seed 45 --chromosomes 1-22 \
    --chr-length-mb 70 --output-dir ~/out \
    --profile-memory ~/out/memprof.tsv
```

The TSV has six columns:

| Column | What it is |
|---|---|
| `elapsed_s` | seconds since profiler started |
| `rss_mb` | parent process RSS in MB |
| `vms_mb` | parent process VMS in MB |
| `children_rss_mb` | sum of all live descendants' RSS (workers, bcftools subprocesses) |
| `total_rss_mb` | parent + children — usually what you want for OOM diagnosis since the kernel budgets against the whole tree |
| `label` | empty for periodic samples, caller-supplied string for checkpoint marks |

A quick plot:

```python
import pandas as pd, matplotlib.pyplot as plt
df = pd.read_csv("memprof.tsv", sep="\t")
plt.plot(df.elapsed_s, df.total_rss_mb, label="total")
plt.plot(df.elapsed_s, df.rss_mb, label="parent only")
for _, row in df[df.label.notna() & (df.label != "")].iterrows():
    plt.axvline(row.elapsed_s, color="red", alpha=0.3)
    plt.text(row.elapsed_s, row.total_rss_mb, row.label, rotation=90)
plt.legend(); plt.xlabel("seconds"); plt.ylabel("RSS (MB)")
plt.show()
```

The flag is opt-in and zero-cost when not used. Needs
`pip install psutil` (already in `requirements.txt` since Phase
5f). Marks fired from forked workers are silently ignored — only
the parent process writes to the TSV.

---

## 10. Configuration files (optional)

Once your invocation grows past a handful of flags it's easier to
keep them in a file than to retype them every run. The tool reads
**`generate_people_config.yaml`** automatically when it exists in the
directory you run from. CLI flags still win over config values, so
you can keep a stable baseline file and override one or two settings
on the command line for a single experiment.

### Bootstrap a starter config with `--print-config`

The fastest way to create your first config file is to ask the
tool to print one. `--print-config` emits a fully-valid YAML
config with every field set to its built-in default, each preceded
by a `# description` comment drawn from the schema:

```bash
.venv/bin/python generate_people.py --print-config > generate_people_config.yaml
```

The emitted file is a **no-op as-is** — running the tool with no
other flags after that command behaves identically to running with
no config at all. From there, edit only the values you want to
change; the comments tell you what each field controls and the
`# yaml-language-server: $schema=...` line wires up IDE
auto-complete (see [Editor integration](#editor-integration-vs-code-intellij-helix-) below).

A snippet of what gets emitted:

```yaml
# generate_people_config.yaml
#
# Starter config emitted by `generate_people --print-config`.
# Every field is set to its built-in default …

# yaml-language-server: $schema=./generate_people_config.schema.json

schema_version: 1

cohort:
  # Cohort size (number of person VCFs).
  n: 10
  # Reference build assembly.
  build: GRCh38
  # Master RNG seed; omit for fresh randomness each run.
  seed: null
  # Chromosomes spec: list / range / mix (e.g. '22', '19-22,X').
  chromosomes: '22'
  # Simulated prefix per chromosome in Mb; 0 = full length.
  chr_length_mb: 5.0

# … one section per area: simulation, overlays, structural_variants,
# sequencing_errors, performance, output, admixture, legacy_background
```

`--print-config` writes only to stdout and exits 0 without touching
any other state, so it composes cleanly with shell redirection,
diffing against an existing config (`diff <(generate_people --print-config) generate_people_config.yaml`),
or piping into another tool.

### Quick start (hand-written)

If you'd rather hand-write a minimal config, create
`generate_people_config.yaml` next to where you run from:

```yaml
# yaml-language-server: $schema=./generate_people_config.schema.json
schema_version: 1

cohort:
  n: 3000
  seed: 42
  chromosomes: "1-22"
  chr_length_mb: 70

performance:
  workers: 8
  cohort_mode: arrow
```

Run with no flags — the tool picks it up:

```bash
.venv/bin/python generate_people.py
```

You should see on stderr:

```
  Loading values from config file: generate_people_config.yaml
  Effective non-default values:
    n                            = 3000                 [config]
    seed                         = 42                   [config]
    chromosomes                  = '1-22'               [config]
    chr_length_mb                = 70.0                 [config]
    workers                      = 8                    [config]
    cohort_mode                  = 'arrow'              [config]
```

Override one value for a single run:

```bash
.venv/bin/python generate_people.py --workers 4
```

stderr now reports the source mix:

```
    workers                      = 4                    [cli, overrides config value 8]
    n                            = 3000                 [config]
    ...
```

### What every field looks like

The complete set of keys, types, defaults, and bounds is in
[`generate_people_config.schema.json`](generate_people_config.schema.json).
Sections in the YAML mirror the CLI groupings:

| Section | Covers |
|---|---|
| `cohort` | `n`, `build`, `seed`, `chromosomes`, `chr_length_mb` |
| `simulation` | `demo_model`, `population`, `rec_rate`, `mu` |
| `overlays.clinvar` / `overlays.rsid` / `overlays.cosmic` | density + source overrides for each annotation source |
| `structural_variants` | per-person count + length bounds |
| `sequencing_errors` | GT flip rate + dropout rate |
| `performance` | `workers`, `cohort_mode`, `cohort_arrow_batch_size`, `fanout_batch_size`, `chr_chunk_mb`, `no_resume`, `profile_memory` |
| `output` | `dir`, `cache_dir`, `mode` |
| `admixture` | `enabled` + EUR/SAS/AFR fractions (must sum to 1.0 when enabled) |
| `legacy_background` | M4 legacy 1000G-pool sampler config |

Every section is optional. A minimal valid config is just
`schema_version: 1` — every default applies and the run behaves
identically to no config at all.

### Editor integration (VS Code, IntelliJ, Helix, …)

The `# yaml-language-server: $schema=...` comment at the top of the
example wires up live validation. With the
[YAML extension](https://marketplace.visualstudio.com/items?itemName=redhat.vscode-yaml)
installed in VS Code, you get:

- Red squiggles on typos (`cohort_mode: arow` flagged immediately).
- Hover tooltips with each field's description.
- Auto-complete for the section names and enum values.

The schema file lives next to the docs as
`synthetic_people/generate_people_config.schema.json`; it's
auto-generated from the Pydantic models in
`syntheticgen/config.py` and a CI test fails if the two ever drift.

### Validation errors

The config is fully validated before any simulation work runs. If
something's wrong, every problem is reported at once so you don't
have to fix one, re-run, see the next, fix it, re-run, etc.:

```
Config validation failed in generate_people_config.yaml:
  - cohort.n: input should be greater than or equal to 1
  - performance.cohort_mode: input should be 'auto', 'sites_list', or 'arrow'
  - overlays.clinvar.inject_density: input should be less than or equal to 1
```

### Precedence rules

Three sources, in strictly decreasing priority:

1. **CLI flag** you explicitly typed
2. **Config-file value** (if present)
3. **Built-in default**

`--no-config` skips the auto-discovery if you want a config-free
run from a directory that happens to have a config file in it.
Passing `--config /path/to/other.yaml` overrides discovery and
loads that file instead. Missing `--config` path is a hard error;
missing auto-discovery file just falls back to CLI + defaults
silently.

### Versioning and forward compatibility

`schema_version: 1` is required. The loader rejects configs with
unknown schema versions with a clear message, so a future incompat-
ible change can never silently mis-interpret an old config — you'll
be told to update the schema and check the changelog. Today the
only supported value is `1`.

---

## 11. Troubleshooting

**"`bcftools not found`"** — install htslib:
`apt install bcftools tabix` (Linux) or
`brew install htslib bcftools` (macOS).

**"`ModuleNotFoundError: No module named 'msprime'`"** — your venv
is missing the M5+ deps. Re-run
`.venv/bin/pip install -r synthetic_people/requirements.txt`.

**"Optional Python deps (not required for M1)"** in `--check-deps`
output — informational. The defaults need `numpy` / `msprime` /
`stdpopsim` / `matplotlib` / `scikit-allel`; the legacy path
(`--legacy-background`) needs only `numpy`.

**Run hangs or eats RAM** — drop `--chr-length-mb` and `--n`.
Coalescent simulation cost scales with both. For full-chromosome
sims of large autosomes (`--chr-length-mb 0` on chr1–chr5) you'll
want a real machine, not a laptop.

**ClinVar download fails** — the pipeline tries
`https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`
once. If you're behind a proxy, download manually and place it as
`synthetic_people/cache/clinvar_GRCh38.vcf.gz` (and `.tbi`).

**`--somatic` rejected** — pass `--cosmic-vcf PATH`. We never
auto-fetch COSMIC; you need a registered download.

**`--art` rejected** — the heavy ART-based read simulation path is
gated until the GRCh38 reference FASTA is wired in. Use the default
`--error-rate` / `--dropout-rate` lightweight model for now.

**Single tight PCA cluster, no structure** — that's the right
behaviour for a single-population run. Use `--admixture` if you want
visible PCA structure; or generate two cohorts with different
`--population` and concatenate.

**Truth BED rows look out of order** — they are sorted by
`(contig_index, chrom, start)`, which matches genome order under
`BUILDS[args.build]["contigs"]` (chr1, chr2, ..., chrX, chrY, chrM).
A simple alphabetic `sort` will *not* match — use
`sort -k1,1V -k2,2n` if you need to merge with another tool.

---

## 12. Glossary

| Term | Meaning |
|---|---|
| **AD** | Allelic depth (per-allele read count). FORMAT tag `AD`. |
| **Admixture pulse** | A demographic event where one or more source populations contribute ancestry to a new deme at a single point in time. The pipeline uses a 20-generation-old single pulse for the UK cohort. |
| **Cohort site** | A genomic position emitted into all VCFs in a run; per-person genotypes vary, but `pos`/`ref`/`alt` are shared. |
| **DP** | Read depth. FORMAT tag `DP`. |
| **GQ** | Genotype quality (Phred-scaled). Recomputed from AD, so flipped GTs have low GQ. |
| **Highlighted variant** | The single ClinVar-pathogenic variant flagged `INFO/HIGHLIGHT=1` per person — the headline event each generated person carries. |
| **LD** | Linkage disequilibrium — non-random association of alleles at nearby loci. Measured here as r² across distance bins. |
| **Local ancestry** | The source population of a *segment* of one haplotype (vs. a person's overall ancestry fractions). The `ancestry/` BEDs encode this. |
| **r²** | Squared correlation between alt-dosage vectors at two sites. Standard LD measure. |
| **SFS** | Site-frequency spectrum — distribution of derived allele counts across the cohort. The pipeline targets a power-law `1/k^α` shape (default α = 2.0). |
| **stdpopsim** | The community catalogue of validated demographic models for human population genetics; provides `OutOfAfrica_3G09`, `Africa_1T12`, etc. |
| **SVTYPE / SVLEN / END** | Standard VCF 4.2 INFO tags for structural variants. |
| **Ti/Tv** | Transition / transversion ratio; ~2.1 for whole-genome human data. |

For the architecture of the package, the per-milestone build history,
and the developer-facing test layout, see `README.md` and
`IMPLEMENTATION_PLAN.md`.
