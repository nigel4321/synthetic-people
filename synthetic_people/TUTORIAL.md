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
10. [Troubleshooting](#10-troubleshooting)
11. [Glossary](#11-glossary)

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

---

## 10. Troubleshooting

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

## 11. Glossary

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
