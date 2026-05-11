# synthetic_people

Generate cohorts of synthetic whole-genome VCFs with realistic variant
content: an **LD-correct coalescent-simulated background** under a
configurable demographic model, **ClinVar-pathogenic variants** highlighted
per person, and full VCF 4.2 **per-sample quality metrics** (GT / DP /
GQ / AD).

The spec is `SYHTHETIC_PROJECT.md`; incremental build plan and
per-milestone status is in `IMPLEMENTATION_PLAN.md`.

As of **M11** every batch ships a self-contained delivery: per-person
`person_NNNN.vcf.gz` + `.tbi`, a pair of BED4 truth tracks
(`out/truth/person_NNNN.golden.bed` for the curated golden set,
`out/truth/person_NNNN.noise.bed` for every M9 flip / dropout), local
ancestry BED (admixture mode), `manifest.json` cataloguing it all, and
the `validation/` artefacts described next.
A smoke test (`scripts/smoke.sh`) exercises the full pipeline in <2
min on a laptop.
**M10** validates every batch end-to-end with
`validate_batch.py`, which walks the per-person VCFs and writes
`summary.json`, a Markdown report, and four PNGs (LD decay, AF
histogram, indel length distribution, cohort PCA) under
`<batch>/validation/`. Acceptance criteria from `SYHTHETIC_PROJECT.md`
§6 — Ti/Tv ≈ 2.1, monotone LD decay, PCA cluster structure on
admixture-mode batches — are visible in those artefacts. **M9** added
a configurable lightweight sequencing-error model: per-call genotype
flips and coverage dropouts applied at a target FDR (~0.1% by
default). Flipped calls land low-GQ because GQ is recomputed from
AD, which still reflects the truth; dropouts emit `./.` with DP /
GQ / AD all zero. **M8** added a
handful of structural variants per person — `<DEL>` / `<DUP>` /
`<INV>` symbolic ALTs with `SVTYPE` / `SVLEN` / `END` / `CIPOS` INFO
tags — alongside the SNV/indel cohort background. **M7** grounds
cohort sites against
public variant databases (ClinVar pathogenic records and dbSNP rsIDs
at real chromosome coordinates). **M6** added the `--admixture` mode
(EUR + SAS + AFR → UK pulse with per-person local-ancestry BED
truth); **M5** is the default single-population coalescent path; the
M4 1000G-pool + power-law SFS sampler is retained behind
`--legacy-background`.

---

## Install

System binaries (not pip-installable):

```
bcftools tabix bgzip      # htslib — header parse, bgzip, tabix-index
```

Python deps go in a project venv (repo root):

```bash
sudo apt install python3.12-venv          # only if ensurepip missing
python3 -m venv .venv
.venv/bin/pip install -r synthetic_people/requirements.txt
```

`requirements.txt` pins:

| Package | Used by | Purpose |
|---|---|---|
| `numpy>=1.24` | M4+ | sampling, histogram arithmetic |
| `msprime>=1.3`, `tskit>=0.5` | M5, M6 | coalescent simulation + tree sequences (M6 also uses `record_migrations` for local ancestry) |
| `stdpopsim>=0.2` | M5 | human demographic catalogue (`OutOfAfrica_3G09`, etc.) |
| `demes` (transitive via msprime) | M6 | UK-cohort admixture demography graph |
| `matplotlib>=3.7`, `scikit-allel>=1.3` | M10 | LD decay r² + cohort PCA + plot artefacts |
| `scikit-learn>=1.3` (transitive) | M10 | PCA decomposition for the cohort matrix |
| `pydantic>=2.0`, `PyYAML>=6.0` | config file | Validate and load `generate_people_config.yaml` (see TUTORIAL.md §10) |

**Optional dep — `pyarrow` (cohort scale ≥ n=100 000):**

```bash
.venv/bin/pip install pyarrow
```

Required only when `--cohort-mode` resolves to `arrow` (the auto-pick
above `--n 100000`, or any explicit `--cohort-mode arrow`). The default
`--cohort-mode auto` falls back to the in-RAM `sites_list` path with
a one-line warning when pyarrow is missing, so runs at smaller `n`
keep working without it. See TUTORIAL.md §9.1.1 for the mode picker
details and at-scale disk requirements.

Probe the environment:

```bash
.venv/bin/python synthetic_people/generate_people.py --check-deps
```

First run only: ClinVar is downloaded (~50 MB) into
`synthetic_people/cache/` and re-used thereafter.

---

## Usage

### Default: coalescent-simulated cohort (M5)

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 50 --seed 42 \
    --chromosomes 22 --chr-length-mb 5.0 \
    --demo-model OutOfAfrica_3G09 --population CEU
```

- `--n` — cohort size (number of person VCFs).
- `--chromosomes` — list, range, or both (e.g. `22`, `19,20,21,22`,
  `1-22`, `1-3,5,19-22,X`). Ranges must be numeric and inclusive.
- `--chr-length-mb` — simulated prefix per chromosome, 0 = full length.
- `--demo-model` — stdpopsim model id; `none` falls back to constant-size
  Ne = 10 000 msprime draw.
- `--population` — sampling population for the demo model (CEU, YRI, CHB
  for `OutOfAfrica_3G09`).
- `--rec-rate`, `--mu` — only consulted when `--demo-model=none`.

### Admixture: EUR + SAS + AFR → UK pulse (M6)

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 50 --seed 42 \
    --chromosomes 22 --chr-length-mb 5.0 \
    --admixture --eur-frac 0.60 --sas-frac 0.25 --afr-frac 0.15
```

- `--admixture` — overrides `--demo-model` / `--population`; runs the
  three-source UK-cohort demography (M6).
- `--eur-frac` / `--sas-frac` / `--afr-frac` — per-source ancestry
  proportions. Must sum to 1.0 and be non-negative. Defaults
  60/25/15.

In addition to the per-person VCFs, this mode emits one
`out/ancestry/person_NNNN.bed` truth track per individual (columns:
`chrom  start  end  hap1_pop  hap2_pop`) and an `out/manifest.json`
including realised per-person ancestry fractions.

### Database grounding: ClinVar + dbSNP + COSMIC (M7)

Both the coalescent and admixture paths apply the same overlay
sequence to the cohort sites before per-person emission. Defaults give
~20% rsID density and a small ClinVar-injection fraction at realistic
chromosome coordinates; tune with:

| Flag | Effect |
|---|---|
| `--rsid-density 0.20` | Fraction of cohort sites rewritten to a known dbSNP variant (real `pos` / `ref` / `alt` / `rsID`). `0` disables. |
| `--dbsnp-vcf PATH` | Optional override for the rsID source. Default is the cached ClinVar VCF whose `INFO/RS` tag carries dbSNP rs numbers — no extra download. Pass a real dbSNP VCF (rsIDs in the ID column) to use a richer pool. |
| `--clinvar-inject-density 0.01` | Fraction of cohort sites overwritten with a random ClinVar pathogenic record. Lands `CLNSIG` / `CLNDN` on a handful of records per person. `0` disables (the per-person highlighted variant still lands). |
| `--somatic --cosmic-vcf PATH` | Opt-in COSMIC overlay (`COSMIC_ID` / `COSMIC_GENE` INFO tags). COSMIC is registration-gated, so we never auto-fetch — supply a local file. |
| `--cosmic-inject-density 0.005` | Fraction of cohort sites rewritten with COSMIC records when `--somatic` is on. |

Overlays operate on disjoint sites — ClinVar-injected rows are
reserved against rsID injection so each row carries at most one
overlay. Cohort GT blocks (the LD signal) are preserved across
injection; only `pos` / `ref` / `alt` / `id` and the new INFO tags
change. Run summary and `out/manifest.json` record realised counts of
each overlay.

### Structural variants (M8)

Every per-person VCF picks up a handful of SVs by default:

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 5 --seed 42 \
    --chromosomes 22 --chr-length-mb 5.0 \
    --svs-per-person 5 --sv-length-min 50 --sv-length-max 10000
```

| Flag | Effect |
|---|---|
| `--svs-per-person 3` | Number of SVs emitted per person. `0` disables SV emission. |
| `--sv-length-min 50` | Minimum SV length in bp (log-uniform). |
| `--sv-length-max 10000` | Maximum SV length in bp; positions are drawn so `END = POS + length` stays inside the simulated span. |

Type mix: ~50% `DEL`, ~30% `DUP`, ~20% `INV`. `INFO/SVLEN` is negative
for `DEL` and positive for `DUP` / `INV`; `INFO/END = POS + |SVLEN|`;
`INFO/CIPOS = -50,50` (every SV is currently flagged "imprecise").
Anchor `REF` is a random standard base — the real GRCh38 reference
isn't on disk in M8; M11 will wire the reference FASTA in for
exact-anchor reporting.

`out/manifest.json` gains an `svs` block recording per-person count,
length range, and cohort-total SVs emitted; per-person entries record
their `n_svs`.

### Sequencing errors (M9)

Every batch passes through a lightweight per-call noise model after
the truth-state DP/AD have been drawn:

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 5 --seed 42 \
    --chromosomes 22 --chr-length-mb 5.0 \
    --error-rate 0.001 --dropout-rate 0.0005
```

| Flag | Effect |
|---|---|
| `--error-rate 0.001` | Probability of a per-call genotype flip (false positive / negative). 0 disables. |
| `--dropout-rate 0.0005` | Probability of a coverage dropout. Emits `./.:0:0:0,0`. 0 disables. |
| `--art` | Heavy path: ART read simulation + `bcftools call`. Currently rejected with a clear message; needs the M11 GRCh38 reference FASTA. |

Mechanics: GT is perturbed *after* the AD draw, so the recomputed GQ
reflects the disagreement between reads and call. Hom→het is the
dominant flip direction (0.7 weight), matching empirical caller
behaviour; het splits 50/50 between hom-ref and hom-alt; multi-allelic
hets collapse one allele to REF. `out/manifest.json` gains an `errors`
block with requested rates, realised counts, and realised FDR; each
per-person entry carries its own `errors` sub-dict so M11's truth-set
BED can grade calls against the truth.

### Validation suite (M10)

After generating a batch, run `validate_batch.py` against the output
directory:

```bash
.venv/bin/python synthetic_people/validate_batch.py path/to/out/
```

Walks every `person_*.vcf.gz` under the batch dir and produces, under
`<batch>/validation/`:

| Artefact | Contents |
|---|---|
| `summary.json` | Per-sample stats + cohort aggregates (Ti/Tv, Het/Hom, AF histogram, indel/SV breakdown, LD-decay bins, PCA projection) |
| `report.md` | Markdown report linking the plots and surfacing the per-sample table |
| `ld_decay.png` | Mean r² vs. distance, log-x |
| `af_histogram.png` | Allele frequency distribution |
| `indel_lengths.png` | Indel length distribution (insertions positive, deletions negative) |
| `pca.png` | Cohort PCA scatter; admixture-mode batches colour each person by dominant ancestry component |

The validator also reads the batch's `manifest.json` if present, so
admixture-mode reports surface requested vs realised ancestry.
Without matplotlib installed, the JSON / Markdown artefacts still
land and a one-line warning skips the plots.

### Truth-set BED tracks (M11)

Every run drops two BED4 files per person under `out/truth/`:

| File | Contents |
|---|---|
| `person_NNNN.golden.bed` | Every record matching the spec's "golden truth" set, tagged with priority `HIGHLIGHTED` > `CLINVAR` > `COSMIC` > `SV` > `RSID`. The 4th column carries a semicolon-separated `flag=…;id=…;ref=…;alt=…;gt=…;…` payload so a downstream caller can split on `flag=` to grade against the model. |
| `person_NNNN.noise.bed` | One row per M9 noise event: `flag=FLIP` or `flag=DROPOUT` plus `truth_gt=…;called_gt=…`. Lets a caller's per-call accuracy be graded against the known noise model. |

Rows are sorted by `(contig_order, chrom, start)` so the BEDs are
`sort -k1,1 -k2,2n`-friendly with no follow-up shell sort. BED
coordinates are 0-based half-open. Manifest entries gain
`golden_bed`, `noise_bed`, `n_golden`, `n_noise` per person; the
top-level `errors` block continues to record the cohort-wide realised
FDR.

### Smoke test (M11)

Quick end-to-end exercise of generation + validation, useful for CI:

```bash
bash synthetic_people/scripts/smoke.sh
```

Generates a 5-person × 0.5 Mb chr22 cohort with default error /
dropout rates, runs the validation suite, and asserts every
deliverable lands on disk. `OUT_DIR`, `N_PEOPLE`, `SEED`, and
`PYTHON` can be overridden by env-var.

### Legacy: 1000G-pool + power-law SFS (M4)

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 50 --seed 42 --legacy-background --build GRCh37 \
    --background-glob "ALL.chr22.phase3_*.vcf.gz" \
    --n-background 500 --sfs-alpha 2.0
```

Retained for comparison, offline-only use, and any workflow that wants
to draw coordinates directly from the local 1000G Phase 3 data.
Legacy-only flags (`--background-glob`, `--n-background`, `--af-min`,
`--sfs-alpha`) are marked `[legacy]` in `--help`.

### Reproducibility

- `--seed N` — same inputs + same seed → byte-identical output, and
  the same regardless of `--workers`. Note: Phase 1 changed how the
  master rng is consumed, so output at a given seed differs from
  pre-Phase-1 runs.
- Omit `--seed` — each invocation produces different people (different
  sample IDs, different highlighted variants, different genotypes).

### Performance / parallelism

- `--workers N` — parallelise the cohort BCF write (sample-slice
  parallelism, post-Phase-5e) and the per-person VCF fan-out across
  `N` worker processes. `0` (default) means auto
  (`os.cpu_count()`); `1` means serial. Linux only — non-Linux hosts
  silently fall back to serial because the parallel path uses
  `fork`-based multiprocessing. **msprime simulation itself runs
  serially across chromosomes** regardless of `--workers` — the
  pre-Phase-5e parallel-chromosome path multiplied tree-sequence
  RAM by the worker count and is no longer used. See TUTORIAL.md
  §9.1 for details.
- `--cohort-mode {sites_list, arrow, auto}` — selects the cohort
  intermediate between simulation and the BCF write. `sites_list`
  (Phase B) keeps the cohort in parent RAM and fork-shares to
  workers; works up to `n ~= 30000`. `arrow` (Phase 5d.1) writes a
  per-chrom Arrow IPC scratch file and workers consume it via
  `pyarrow.memory_map` (zero-copy). `auto` (default) picks `arrow`
  for `--n >= 100000` and `sites_list` below. The Arrow path
  requires `pip install pyarrow`; without it, `auto` falls back to
  `sites_list` with a warning. **At `n=100k+` plan on 2-3 TB free
  scratch on a fast SSD/NVMe filesystem** — the Arrow scratch is
  ~80 GB per chromosome at n=1M and gets created + deleted per
  chromosome, but the merged cohort BCFs accumulate at ~80 GB ×
  22 = ~1.7 TB. The cli runs a pre-flight disk-space check at
  startup when Arrow mode is selected and aborts cleanly if free
  disk can't fit one chromosome's scratch.
- `--cohort-arrow-batch-size N` — sites per Arrow record batch when
  `--cohort-mode` resolves to `arrow`. Default `256`, identified
  empirically by Spike 2b; smaller batches lower parent peak RSS
  linearly at a small (~7%) write-throughput cost. Don't change
  unless you're characterising a specific host.
- The writer streams records straight into `bgzip -c` (no plain `.vcf`
  intermediate), so per-person disk I/O is one pass instead of two.
- The overlay loaders (ClinVar, dbSNP, COSMIC) are bcftools-driven and
  bound on subprocess I/O. They are submitted to a small thread pool
  *before* the coalescent simulation runs, so the loader work overlaps
  with msprime instead of running serially after it. No flag —
  prefetch is always on for the non-legacy path.
- `--mode {per-person, cohort, both}` — output shape selector,
  default `per-person`. As of Phase 5b2 all three modes flow through
  the streamed pipeline on the standard coalescent path: cohort
  chromosomes are simulated → overlaid → written to
  `out/cohort/cohort.chr<N>.bcf` → freed before the next chromosome
  is simulated. Peak RAM is bounded by one chromosome's working set
  rather than the whole cohort.
  - `cohort`: writes the per-chrom BCFs and stops. Per-person VCFs
    can be derived later via `bcftools view -s SAMPLE_ID`.
  - `per-person` (default): per-person VCFs are derived from the
    streamed cohort BCFs via the same `bcftools view -s` pipeline.
    The cohort BCFs land as intermediates alongside the per-person
    VCFs; users who only want per-person can `rm -rf out/cohort`
    after the run, or pass `--mode per-person` plus a wrapper that
    cleans up.
  - `both`: keeps both deliverables explicitly.
- `--no-resume` — on the streamed coalescent path, ignore any
  existing `out/cohort/cohort.meta.json` + cohort BCFs and start a
  fresh simulation. Default behaviour is to resume a prior run when
  its params match — useful for multi-hour cohort runs that get
  interrupted by OOM, SIGINT, or node failure. Param mismatches
  surface a clear error rather than silently re-using incompatible
  state.
- `--chr-chunk-mb N` — split each chromosome's msprime simulation
  into independent sub-chunks of `N` Mb. Bounds per-chunk peak
  RAM for large cohorts. `0` (default) auto-picks a chunk size at
  run start from `psutil.virtual_memory().available` and
  `--workers`, aiming for the per-chunk working set to stay under
  ~50% of free RAM. Positive `N` overrides. Chunks simulate
  independently, so cross-chunk LD is lost — short-range LD inside
  each chunk is preserved, long-range (>chunk size) is not. See
  PERFORMANCE_PLAN §"Phase 5f" for the full memory model.
- Long-running runs print throttled progress lines (~20 s cadence)
  during the cohort BCF write and the per-person fan-out — `cohort
  BCF: 12,345/100,000 sites (4,500/s)` and `person VCFs: 1,234/100,000
  written (5/s, elapsed 240s, eta 19000s)` — so a multi-hour 100k run
  has visible heartbeat without flooding stderr.

---

## Data sources

### ClinVar (highlighted variants)

Downloaded from NCBI once, cached in `./cache/`:

- GRCh37: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz`
- GRCh38: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`

Filtered on `CLNSIG`. Default `Pathogenic`, `Likely_pathogenic`,
`Pathogenic/Likely_pathogenic`; override via `--clinvar-sig`. Delete
`cache/clinvar_<build>.vcf.gz*` to force a refresh.

### 1000 Genomes (legacy background only)

Default glob: `../ALL.chr*.phase3_*.genotypes.vcf.gz`. For each source,
5 000 variants with `MAX(INFO/AF) >= --af-min` (default 0.05) are
reservoir-sampled. Only coordinates + allele strings carry through —
source AFs are ignored downstream because M4's SFS sampler redraws
them per site.

### stdpopsim HomSap catalogue (coalescent default)

`OutOfAfrica_3G09` (CEU / YRI / CHB) drives the demographic history by
default. `Africa_1T12`, `OutOfAfrica_2T12`, etc. are also available.
First coalescent run may download model metadata; subsequent runs are
offline.

---

## Output layout

```
out/
├── person_0001.vcf.gz   # (--mode per-person | both)
├── person_0001.vcf.gz.tbi
├── person_0002.vcf.gz
├── ...
├── cohort/              # streamed per-chrom BCFs (any --mode on coalescent path)
│   ├── cohort.chr1.bcf  # per-chromosome cohort BCF (Phase 5b1+)
│   ├── cohort.chr1.bcf.csi
│   ├── cohort.chr2.bcf
│   ├── ...
│   └── cohort.meta.json # resume-state record (Phase 5b2)
├── manifest.json        # shape, samples[], per-person summary, cohort_bcfs path list
├── ancestry/            # admixture mode only: per-person local-ancestry BEDs
│   ├── person_0001.bed
│   └── ...
├── truth/               # M11: golden + noise BED tracks per person
│   ├── person_0001.golden.bed
│   ├── person_0001.noise.bed
│   └── ...
├── summary/
│   └── sfs.tsv          # cohort AC histogram (columns: ac, n_sites)
└── validation/          # M10 (created by validate_batch.py)
    ├── report.md
    ├── summary.json
    ├── ld_decay.png
    ├── af_histogram.png
    ├── indel_lengths.png
    └── pca.png
```

Which artefacts land depends on `--mode`:

| `--mode` | Per-person VCFs | Cohort BCFs (under `cohort/`) | Manifest fields |
|---|---|---|---|
| `per-person` (default) | yes | per-chrom intermediates: `cohort.chr<N>.bcf` | `shape=per-person`, `samples[]`, `cohort_bcfs[]`, `people[]` |
| `cohort` | — | per-chrom: `cohort.chr<N>.bcf` | `shape=cohort`, `samples[]`, `cohort_bcfs[]`, no `people[]` |
| `both` | yes | per-chrom: `cohort.chr<N>.bcf` | `shape=both`, `samples[]`, `cohort_bcfs[]`, `people[]` |

The top-level `samples[]` list always lands so a downstream tool
wanting "every sample ID" reads one field regardless of `--mode`.
`cohort_bcfs` is a list with one entry per chromosome (parity with
the 1000G layout). All three modes flow through the streamed
pipeline as of Phase 5b2 — `cohort_bcfs` is always populated. The
legacy and admixture paths are exceptions: they keep the in-memory
cohort accumulator and write a single combined `cohort.bcf` only
when `--mode` includes the cohort deliverable.

To derive a per-person VCF from cohort BCFs later (e.g. when `--mode
cohort` was used to keep memory bounded on a 100k-person run):

```bash
# Single chromosome
bcftools view -s HG12345 -Oz \
    out/cohort/cohort.chr22.bcf > out/person_HG12345.chr22.vcf.gz

# All chromosomes — concatenate per-chrom slices
bcftools concat -Oz -o out/person_HG12345.vcf.gz \
    <(bcftools view -s HG12345 -Ou out/cohort/cohort.chr1.bcf) \
    <(bcftools view -s HG12345 -Ou out/cohort/cohort.chr2.bcf) \
    ...
tabix -p vcf out/person_HG12345.vcf.gz
```

Phase 5b2 will provide this derivation as a built-in CLI step.

Per-person BED columns (admixture mode only):

```
chrom    start    end    hap1_pop    hap2_pop
22       0        1234567   EUR         AFR
22       1234567  5000000   EUR         EUR
...
```

`hap1_pop` / `hap2_pop` ∈ {EUR, SAS, AFR, OOA, ANC}. With the default
~600-year-old admixture pulse, the vast majority of segments resolve to
the three source demes; OOA / ANC appear only when an unusually deep
lineage's tree-walk has not yet found an EUR/SAS/AFR ancestor by t=20
gens — reported faithfully when it happens.

Per-person VCF:

- **Header**: `VCFv4.2`, `##reference=<accession URL>`, `##contig` lines
  for every standard chromosome with `assembly=GRCh37` or `GRCh38`, full
  `INFO`/`FORMAT`/`ALT` declarations (AC, AN, AF, SVTYPE, SVLEN, END,
  CIPOS, HIGHLIGHT, CLNSIG, CLNDN; GT/DP/GQ/AD; `<DEL>`/`<DUP>`/
  `<INV>`/`<INS>`).
- **Records**: one ClinVar-highlighted variant flagged `HIGHLIGHT`
  (plus `CLNSIG` / `CLNDN` when present) + the cohort background
  projected to that person (hom-ref calls dropped).
- **Per-record FORMAT** (`GT:DP:GQ:AD`): DP ~ Poisson(λ ≈ 30, per-sample
  jitter σ = 3), AD ~ Binomial(DP, genotype_p) with ~5% ref-bias on
  hets, support-weighted Phred GQ capped [0, 99].
- **Multi-allelic** (legacy path): comma-separated ALT / AC / AF / AD;
  samples can carry `1|2`-style hets.

Compatibility: output passes `nextflow_pipeline/bin/qc_validate.py
--strict` cleanly (GRCh37/GRCh38 recognised, human contigs, GT/DP/GQ/AD
declared, AF + AC/AN present).

---

## Architecture

```
synthetic_people/
├── SYHTHETIC_PROJECT.md
├── IMPLEMENTATION_PLAN.md
├── README.md
├── requirements.txt
├── generate_people.py        # 15-line shim → syntheticgen.cli:main
├── syntheticgen/
│   ├── builds.py             # GRCh37 + GRCh38 contig tables, ClinVar URLs
│   ├── clinvar.py            # M1 + M7 ClinVar fetch + candidate load + cohort overlay/injection
│   ├── dbsnp.py              # M7 rsID injection (default source: ClinVar INFO/RS)
│   ├── cosmic.py             # M7 COSMIC overlay (--somatic, registration-gated)
│   ├── background.py         # 1000G coordinate pool loader (reservoir)
│   ├── cohort.py             # M4 shared-site cohort + haplotype slotting
│   ├── sfs.py                # M4 P(k) ∝ 1/k^α sampler + histogram
│   ├── coalescent.py         # M5 msprime + stdpopsim driver
│   ├── admixture.py          # M6 UK-cohort demes pulse + local ancestry
│   ├── sv.py                 # M8 structural variant generator (DEL/DUP/INV)
│   ├── errors.py             # M9 lightweight per-call noise (GT flips + dropouts)
│   ├── truth.py              # M11 golden + noise BED4 truth-set writer
│   ├── validate.py           # M10 stats / LD decay / PCA primitives
│   ├── plots.py              # M10 matplotlib plot helpers
│   ├── titv.py               # M3+ Ti/Tv calibrator for de-novo SNVs
│   ├── quality.py            # M2 DP / GQ / AD simulation
│   ├── header.py             # VCF header assembly
│   ├── writer.py             # bgzip + tabix single-sample write
│   └── cli.py                # argparse + orchestration
├── validate_batch.py         # M10 top-level CLI → out/validation/
├── tests/
│   ├── test_quality.py       # M2 + M3 (generalised to N alleles)
│   ├── test_multiallelic.py  # M3
│   ├── test_titv.py          # M3+
│   ├── test_sfs.py           # M4
│   ├── test_cohort.py        # M4
│   ├── test_coalescent.py    # M5 (skips cleanly without msprime/stdpopsim)
│   ├── test_admixture.py     # M6 (skips cleanly without msprime/demes/tskit)
│   ├── test_overlays.py      # M7 (pure-Python; no bcftools / network)
│   ├── test_sv.py            # M8 (pure-Python; no bcftools / network)
│   ├── test_errors.py        # M9 (pure-Python; no bcftools / network)
│   ├── test_validate.py      # M10 (gates plot/PCA tests on matplotlib/sklearn)
│   └── test_truth.py         # M11 (pure-Python; no bcftools / network)
├── scripts/
│   └── smoke.sh              # M11 end-to-end CI smoke test
└── out/                      # generated VCFs + truth/ + summary/ + ancestry/ + validation/
```

---

## Functionality by milestone

### M1 — package scaffolding, GRCh38 default, rich header

Split the original single-file `generate_people.py` into the
`syntheticgen/` package (`generate_people.py` is now a thin shim).
Default `--build` flipped to GRCh38; GRCh37 retained. `##contig` lines
carry `assembly=GRCh38` (or `GRCh37`). FORMAT declares the full
GT/DP/GQ/AD tag set from day one; INFO declares SVTYPE / SVLEN / END /
CIPOS plus symbolic `<DEL>`/`<DUP>`/`<INV>`/`<INS>` ALTs — populated in
later milestones, declared early so the header stops changing shape.
`--check-deps` audits htslib binaries and Python deps.

### M2 — per-variant quality metrics (DP / GQ / AD)

`syntheticgen/quality.py` simulates:

- **DP** ~ Poisson(λ = 30, per-sample jitter σ = 3).
- **AD** ~ Binomial(DP, p) with p = 0.0 / 0.475 / 1.0 by genotype (0.475
  on hets models empirical ref-bias in short-read WGS).
- **GQ** — support-weighted Phred, capped at 99, with depth-dependent
  ceiling `10 · log10(DP) · 6`.

Writer emits `GT:DP:GQ:AD`. 5-person batch verified: DP mean 29.4,
AD `sum == DP` on 100% of rows, het alt-fraction 0.475.

### M3 — indels, multi-allelic, Ti/Tv in target range

Indels flow through unchanged from the 1000G source (already
left-aligned, parsimonious — verified no indel has a common prefix
> 1 char across a 50-person batch).

Multi-allelic opened end-to-end: loader keeps per-allele AFs, draw
samples two haplotypes from a categorical over
`{REF, alt_1, …, alt_k}` so `1|2` hets can occur, per-allele AC / AF /
AD, `Number=R` for AD. 50-person chr22 batch Ti/Tv = 2.11 naturally
from the 1000G source.

### M3+ — Ti/Tv calibrator

`syntheticgen/titv.py` — `choose_alt(ref, rng, target=2.1)` draws a
non-REF base weighted so long-run Ti/Tv converges on `target`
(transition partner weight `target`; each transversion weight `0.5`).
Landed defensively ahead of M5, where de-novo SNV generation drops the
unbiased ratio to ~0.5. `is_transition` / `titv_ratio` helpers for
downstream validation.

### M4 — cohort-level generation + power-law SFS

Pivot from per-person independent HWE draws to one-pass cohort
simulation:

- `syntheticgen/sfs.py` — `draw_minor_count(n_hap, α)` samples k ∈
  {1, …, 2N-1} with `P(k) ∝ 1/k^α` (default α = 2.0). Steeper than
  Watterson (α = 1.0) to match gnomAD-like singleton-dominated spectra.
  `draw_allele_counts` handles multi-allelic via rejection so total
  AC ≤ 2N-1.
- `syntheticgen/cohort.py` — `assign_haplotypes` places alt alleles
  into specific 2N haplotype slots without replacement. Diploid GTs per
  person come from pairing consecutive slots, so realised AC matches
  drawn AC exactly (no HWE-resampling smoothing at the site level).

SFS histogram persisted to `out/summary/sfs.tsv`; singleton count +
fraction printed to the run log. 50-person chr22 legacy batch: 317
singletons / 511 alt observations = **62% singleton fraction** (clears
the >50% exit threshold).

### M5 — coalescent backbone (msprime + stdpopsim)

`syntheticgen/coalescent.py` drives msprime through stdpopsim's engine;
the chosen demographic model supplies population-size history and
per-chromosome metadata. REF/ALT bases for binary tree-sequence
mutations are synthesised via the M3+ Ti/Tv calibrator. Output matches
the M4 cohort-site dict shape, so the writer is unchanged.

New flags: `--chromosomes`, `--chr-length-mb`, `--demo-model`
(default `OutOfAfrica_3G09`), `--population` (default CEU), `--rec-rate`,
`--mu`, `--legacy-background`.

200-sample × 10 Mb chr22 exit check: 28 054 variable sites in 16 s;
8 805 common (MAF ≥ 5%); monotonic LD decay:

| distance bin | mean r² |
|---|---|
| 100–500 bp | 0.55 |
| 0.5–1 kb | 0.46 |
| 1–5 kb | 0.35 |
| 5–20 kb | 0.20 |
| 20–100 kb | 0.05 |
| 100–500 kb | 0.006 |
| ≥ 500 kb | <0.003 |

r² < 0.1 reached by ~20 kb — well inside the "<0.1 by 1 Mb" plan
threshold. Short-range anchor sits at ~0.5 rather than the plan's 0.9
because recombination is uniform; wiring in `HapMapII_GRCh38` hit a
stdpopsim "missing data" error on sub-chromosome regions and is deferred
to M6/M10.

### M6 — UK-cohort admixture + local ancestry truth

`syntheticgen/admixture.py` builds a `demes`-defined demography with
three source demes (EUR, SAS, AFR) joining at a single admixture pulse
into a UK deme **20 generations (~600 years) ago**. Source population
sizes mirror the Gutenkunst OOA_3G09 parameterisation
(ANC = 12.3 k, AFR = 12.3 k, OOA bottleneck = 2.1 k, EUR / SAS = 10 k);
present-day UK Ne = 50 k. Mutations come from `BinaryMutationModel`
with REF/ALT bases drawn through the M3+ Ti/Tv calibrator, so the
output `sites` list has the exact shape M4 / M5 produce (writer is
unchanged).

Local ancestry: `msprime.sim_ancestry(..., record_migrations=True)`
records every t = 20 migration. For each haplotype-sample we walk the
tree at every breakpoint to find the lineage node spanning the pulse
time, then look up which source deme it migrated into. Adjacent
same-ancestry segments are merged; haplotype pairs are then
intersected into per-person joint
`(start, end, h1_pop, h2_pop)` rows, written one BED per person to
`out/ancestry/person_NNNN.bed`.

`out/manifest.json` lists each person with VCF path, BED path,
highlighted ClinVar variant, background record count, and realised
ancestry fractions. Top-level fields capture the requested
`ancestry_proportions` and the run mode (`coalescent` / `admixture-uk`
/ `legacy-background`).

20-person × 5 Mb chr22 exit check (seed 42, default 60/25/15):
13 549 variable sites; 43 ancestry segments across the cohort
(mean 2.1 segments/person — biologically expected because 20
generations × 5 Mb yields ≈ 1 recombination breakpoint per haplotype);
cohort-mean realised ancestry **EUR = 0.456, SAS = 0.352, AFR = 0.192**
— within finite-cohort sampling noise of the requested 0.60 / 0.25 /
0.15 mix. Per-person VCFs pass `qc_validate.py --strict` with 0
errors / 0 warnings.

A 30-person × 1 Mb stand-alone proportions check
(`tests/test_admixture.py::test_ancestry_fractions_track_requested_proportions`)
lands EUR ≈ 0.6, SAS ≈ 0.25, AFR ≈ 0.15 within ±15%. The literal PCA
acceptance test in spec §6 lands in M10.

### M7 — ClinVar / dbSNP / COSMIC grounding

`syntheticgen/clinvar.py` gains `load_clinvar_index`,
`annotate_clinvar` (collision-only) and `inject_clinvar`
(coordinate-replacing). Coalescent positions live in `[1, sim_length]`
while ClinVar sits at real chromosome coordinates (chr22 ClinVar
spans 15.5 M – 50.8 M), so collision-only annotation almost never
fires; `inject_clinvar` is the practical mechanism for landing
CLNSIG / CLNDN at realistic positions. Cohort GT blocks survive
injection — only `pos` / `ref` / `alt` / `id` and the INFO tags are
overwritten.

`syntheticgen/dbsnp.py` exposes `load_rsid_pool` and `inject_rsids`.
The default rsID source is the cached ClinVar VCF whose `INFO/RS` tag
already carries dbSNP rs numbers (thousands of records per chromosome,
no extra download); `--dbsnp-vcf PATH` accepts any dbSNP-style file
where rsIDs sit in the ID column. `_normalise_rsid` handles both
shapes and bare-digit / prefixed / semicolon- or comma-listed values.

`syntheticgen/cosmic.py` overlays a user-supplied COSMIC VCF behind
`--somatic --cosmic-vcf PATH`; never auto-fetches because COSMIC is
registration-gated. `inject_cosmic` lands COSMIC_ID / COSMIC_GENE INFO
tags onto a configurable fraction of sites.

The three overlays operate on disjoint cohort rows: each pass reserves
already-claimed indices so no row carries conflicting annotations.
Header gains COSMIC_ID / COSMIC_GENE INFO declarations alongside
CLNSIG / CLNDN; writer carries every annotation field through the
per-person record onto the emitted INFO field.

5-person × 1 Mb chr22 exit check (`--demo-model none`, seed 42,
default densities): 1,299 cohort sites; **13 ClinVar pathogenic
injections** at real chr22 coordinates (e.g.
chr22:29673446 `Pathogenic / Neurofibromatosis,_type_2`);
**259 rsID injections** (~20% of records, e.g.
chr22:15528207 `rs3924507 C>T`). Per-person VCFs carry **117–135
rsIDs** and **5–7 CLNSIG-bearing records** each. All 5 VCFs pass
`qc_validate.py --strict` with 0 errors / 0 warnings.

### M8 — Structural variants

`syntheticgen/sv.py` emits a handful of SVs per person — DEL / DUP /
INV — with proper VCF 4.2 symbolic ALTs and the full SV INFO tag set.
`generate_person_svs(rng, chromosomes, chrom_length_bp, n_svs,
length_min_bp, length_max_bp)` draws lengths log-uniformly on
[min, max] (default 50 bp – 10 kb) with type weights 50/30/20%
(DEL/DUP/INV). Anchor `REF` is a single standard base placeholder —
the real GRCh38 FASTA isn't loaded in M8. SVs are emitted as part of
each person's `background` record list and flow through the standard
writer, which detects `variant["svtype"]` and emits
`SVTYPE / SVLEN / END / CIPOS` after the AC/AN/AF/HIGHLIGHT/CLNSIG
fields.

3-person × 1 Mb chr22 exit check (seed 42,
`--svs-per-person 5`): each VCF carries exactly 5 SVs;
`bcftools view -i 'INFO/SVTYPE!="."'` returns the SV records;
`bcftools stats` counts them under "number of others". Spot check:
`chr22:535206 A→<DEL>` with `SVTYPE=DEL;SVLEN=-586;END=535792;
CIPOS=-50,50;GT=0|1`. All three VCFs pass `qc_validate.py --strict`
with 0 errors / 0 warnings.

### M9 — Sequencing error modelling

`syntheticgen/errors.py` injects per-call genotype flips and coverage
dropouts at configurable rates. The perturbation is applied **after**
the truth-state AD has been drawn, so a flipped GT lives with a low
recomputed GQ — the realistic mis-call signal where the reads
disagree with the call. `maybe_flip_gt` weights biallelic flips
toward hom→het (0.7) over hom→opposite-hom (0.3); het splits 50/50
between hom-ref and hom-alt; multi-allelic `1|2`-style hets collapse
one allele to REF. `maybe_dropout` zeros DP/AD/GQ and emits `./.`.

Writer threads `error_rate`, `dropout_rate`, and an optional `stats`
dict; CLI accumulates per-person counters into a manifest `errors`
block. Heavy-path `--art` is gated and rejected for now: it needs
the M11 GRCh38 reference FASTA to feed read simulation.

3-person × 1 Mb chr22 exit check (seed 42, `--error-rate 0.01
--dropout-rate 0.005`): **realised FDR 1.55%** vs requested 1.50%;
16 flips + 13 dropouts over 1,871 calls. Mis-call signal visible in
the output — e.g. `chr22:758982 0|0:40:0:0,40` (called hom-ref but
all 40 reads on the alt, GQ correctly drops to 0). Dropouts emit
`./.:0:0:0,0`. All three VCFs pass `qc_validate.py --strict` with
0 errors / 0 warnings.

### M10 — Validation suite

`syntheticgen/validate.py` provides the analytics primitives;
`syntheticgen/plots.py` keeps all matplotlib code in one gated module;
`validate_batch.py` is the top-level CLI. The validator walks every
`person_*.vcf.gz` under the batch directory, computes per-sample +
cohort stats (Ti/Tv, Het/Hom-alt, AF, indel lengths, SVs, singletons,
dropouts), then builds a genotype-dosage matrix to compute LD decay
and PCA. Output lands as `summary.json`, `report.md`, and four PNGs
under `<batch>/validation/`.

LD decay uses log-spaced distance bins (100 bp – 500 kb) and samples
SNP pairs reproducibly under a seedable RNG. PCA mean-imputes missing
genotypes, prunes zero-variance columns, and runs `sklearn`'s PCA on
the cohort matrix. Admixture-mode batches automatically label each
sample by dominant ancestry component so the spec-§6 "cluster or
bridge clusters" criterion is testable visually.

Single-population exit check (30 people × 5 Mb chr22, `--demo-model
none`, seed 42): **Ti/Tv = 1.822** ✓, **Het/Hom-alt = 2.004** ✓,
LD decay monotone (0.327 → 0.123 across 100 bp – 500 kb), PCA PC1
captures 5.3% of variance — appropriately low for a single-pop
constant-Ne draw.

Admixture exit check (30 people × 5 Mb chr22, default 60/25/15
EUR/SAS/AFR, seed 42): **Ti/Tv = 1.882** ✓, LD decay monotone
(0.338 → 0.104), **PCA PC1 captures 19.6% of variance** — the clear
ancestry signal the spec calls for, with EUR / SAS / AFR-dominant
labels visible as separable clusters in `pca.png`.

### M11 — Delivery packaging

`syntheticgen/truth.py` adds `TruthBedWriter`, which emits two BED4
truth tracks per person under `out/truth/`:

- `person_NNNN.golden.bed` — every record matching the spec's
  "golden truth" set, tagged with priority `HIGHLIGHTED` >
  `CLINVAR` > `COSMIC` > `SV` > `RSID`. Each row carries a
  semicolon-separated `flag=…;id=…;ref=…;alt=…;gt=…;…` payload in
  the BED4 name column so a downstream caller can split on `flag=`
  to grade against the model. ClinVar rows additionally surface
  `clnsig` / `clndn`; COSMIC rows surface `cosmic_id` / `cosmic_gene`;
  SV rows surface `svtype` / `svlen`.
- `person_NNNN.noise.bed` — one row per M9 noise event (`flag=FLIP`
  or `flag=DROPOUT`) with `truth_gt=…;called_gt=…` so a caller's
  per-call accuracy can be graded against the known noise model.

Rows are buffered in memory and sorted by `(contig_order, chrom,
start)` on close, so the BED is `sort -k1,1 -k2,2n`-friendly with
no follow-up shell sort. BED coordinates are 0-based half-open
(SNV at 1-based pos 1000 → `[999, 1000)`; a 4-base deletion at
pos 1000 → `[999, 1003)`).

The writer is created per-person in `cli.py` and threaded through
`write_person_vcf`, so it sees both golden and noise events as the
per-record loop runs. Manifest entries gain `golden_bed`,
`noise_bed`, `n_golden`, and `n_noise` fields.

`scripts/smoke.sh` runs an end-to-end 5-person × 0.5 Mb chr22 cohort
plus the validation suite and asserts every advertised deliverable
lands on disk (VCF + tbi, both BEDs per person, manifest, summary,
all four validation PNGs, report.md, summary.json). Defaults to a
deterministic seed and finishes in <2 min on a laptop.

5-person × 0.5 Mb chr22 exit check (`scripts/smoke.sh`, seed 42,
default error / dropout rates): every deliverable lands without
manual intervention; realised FDR **0.293%** over 1,365 calls
(2 flips + 2 dropouts). Per-person golden BEDs carry 56–69 rows
each (highlighted + injected ClinVar + injected rsIDs + 3 SVs);
noise BEDs carry 0–2 rows. Manifest exposes the new
`golden_bed` / `noise_bed` / `n_golden` / `n_noise` fields per
person.

---

## Test suite

206 tests across twelve files; all passing with deps installed.

```bash
cd synthetic_people && ../.venv/bin/python -m unittest discover -s tests -v
```

Without msprime / stdpopsim / demes / tskit / matplotlib / sklearn
installed, the corresponding tests skip cleanly and **165/165**
remaining still pass (`test_overlays.py`, `test_sv.py`,
`test_errors.py`, and `test_truth.py` are pure-Python; numpy-only
validate tests still run if numpy is on PATH; matplotlib/sklearn-gated
subset of `test_validate.py` skips when those deps are absent).

| File | Count | Coverage |
|---|---|---|
| `test_quality.py` | 12 | Poisson DP distribution, bi/multi-allelic AD (`sum==DP`, ref-bias on `0\|1` / `0\|2`, 50/50 split on `1\|2`), GQ range, depth-dependent cap, genotype consistency |
| `test_multiallelic.py` | 5 | Bi-allelic HWE reduction, multi-allelic categorical, `1\|2`-style het rate, per-alt dosage vectors |
| `test_titv.py` | 14 | Transition-partner table, transversion enumeration, case-insensitivity, `titv_ratio` corner cases (empty / no-Tv / indel skip), `choose_alt` uniformity, convergence at targets 0.5 / 1.0 / 2.1 / 3.0 (±5%), parameter validation |
| `test_sfs.py` | 16 | `draw_minor_count` range + near-uniform at α→0, singleton fraction >55% at default α = 2.0, `draw_allele_counts` total bound, histogram aggregation, TSV round-trip, parameter validation |
| `test_cohort.py` | 14 | `assign_haplotypes` exact-count preservation, random-slot placement, overflow rejection, cohort reproducibility under seed, every-site-variable invariant, coord-sharing across people, hom-ref drop-out |
| `test_coalescent.py` | 10 | Output shape, monotone positions, realised AC = declared AC, no fixed sites, seed reproducibility, Ti/Tv ∈ [1.7, 2.6], multi-chromosome, error on unknown chromosome, stdpopsim end-to-end (`skipUnless` on msprime/stdpopsim import) |
| `test_admixture.py` | 13 | Demography proportion validation, UK 3-ancestor topology, sites + per-person segments shape, full-chromosome coverage, realised AC = declared AC, BED round-trip, ancestry-fraction normalisation + empty input, multi-chromosome, seed reproducibility, aggregate ancestry tracks requested 60/25/15 within ±15% (`skipUnless` on msprime/demes/tskit import) |
| `test_overlays.py` | 23 | ClinVar `annotate` (collision match, alt mismatch, no-match returns 0); ClinVar `inject` (density count, GT-block preservation, post-sort invariant, off-chromosome skip, zero-density no-op); rsID `_normalise_rsid` (ID-prefixed, bare-digit, INFO/RS fallback, semicolon and comma lists, missing-returns-empty); rsID `inject_rsids` (density, GT preservation, sort invariant, reserve_indices exclusion, zero-density no-op); COSMIC inject (ID + gene + REF/ALT swap, zero-density no-op); ClinVar + rsID overlay disjointness via reserve_indices |
| `test_sv.py` | 22 | `_draw_length` log-uniform skew + bounds + collapsed-range + invalid-range; `_build_sv_record` SVLEN sign by type, anchor base validity, CIPOS default, unknown-SVTYPE rejection; `generate_person_svs` count, zero-returns-empty, well-formed records, length bounds honoured, END inside chrom span, multi-chromosome distribution, type distribution within ±0.07 of (0.50, 0.30, 0.20), seed reproducibility, different-seed divergence, too-small-chrom and empty-chromosome rejection |
| `test_errors.py` | 18 | `maybe_flip_gt` zero/negative-rate no-op, full-rate always-flips, realised flip rate ~1% on 10k draws, biallelic-only flip targets, hom→het 0.7-weight bias, het 0|1 50/50 split between hom-ref and hom-alt, `1\|2` multi-allelic collapse to REF, unparseable GT pass-through, seed reproducibility; `maybe_dropout` zero/full-rate, realised rate, seed reproducibility; `new_error_stats` shape, `merge_stats` in-place add + missing-key seeding; default-constants lock-in |
| `test_validate.py` | 39 | `_parse_info` empty / single / multiple / flag forms; SNV / indel / SV classification; GT dosage hom-ref/het/hom-alt/multi-allelic/missing; dropout detection; Ti/Tv aggregation + zero-Tv + empty corner cases; Het/Hom ratio + zero-hom + zero-both; indel/SV aggregation; AF histogram bin placement + empty; `_r2_pair` perfect-corr / perfect-anticorr / uncorrelated / constant-vector / few-samples / missing-mask; `ld_decay` shape + short-vs-long ordering; cohort PCA on a clear 2-cluster signal (PC1 > 95% variance) + insufficient-columns guard; PNG smoke tests for every plot helper (LD / AF / indel / PCA-handles-None) |
| `test_truth.py` | 20 | `classify_golden` priority (HIGHLIGHTED > CLINVAR > COSMIC > SV > RSID), `.` clnsig treated as missing, unannotated row returns None, GOLDEN_CATEGORIES priority lock-in; `golden_bed_line` half-open SNV interval, deletion ref-extends-end, SV uses explicit end, payload escapes tab / semicolon / newline, ClinVar payload carries clnsig / clndn; `noise_bed_line` flip records both GTs, dropout records `./.`; `TruthBedWriter` sorts by `(contig_order, chrom, start)` on close, count tracking, context-manager protocol, empty writer creates empty files, parent-dir creation |

Per-milestone exit check: `nextflow_pipeline/bin/qc_validate.py --vcf
<person.vcf.gz> --name <id> --out <out.json> --strict` (exit 1 on any
hard failure).

---

## Known gaps

Tracked in `IMPLEMENTATION_PLAN.md`:

- **Heavy `--art` path** — ART read simulation + `bcftools call`
  is gated and exits with a clear message. Wiring it on requires the
  GRCh38 reference FASTA on disk (multi-GB download); deferred until
  there's a concrete need beyond the lightweight noise model.
- **Exact SV anchor REF** — currently a random standard base
  placeholder; same FASTA dependency as `--art`.
- **HapMap recombination map** — coalescent path uses uniform
  recombination, so short-range r² caps at ~0.55 vs the spec's
  aspirational 0.9. `HapMapII_GRCh38` hit a stdpopsim "missing data"
  error on sub-chromosome regions; revisit when running full-chrom
  sims.

## CLI reference

| Flag | Purpose | Default |
|---|---|---|
| `--n` | Cohort size | `10` |
| `--seed` | RNG seed; omit for fresh randomness each run | `None` |
| `--build` | `GRCh37` or `GRCh38` | `GRCh38` |
| `--output-dir` | Per-person VCF output | `./out` |
| `--cache-dir` | ClinVar download cache | `./cache` |
| `--check-deps` | Print dependency status and exit | `False` |
| `--clinvar-sig` | Comma-separated CLNSIG values | `Pathogenic,Likely_pathogenic,Pathogenic/Likely_pathogenic` |
| `--clinvar-inject-density` | [M7] Fraction of cohort sites overwritten with random ClinVar pathogenic records | `0.01` |
| `--rsid-density` | [M7] Fraction of cohort sites overwritten with a known dbSNP variant + rsID | `0.20` |
| `--dbsnp-vcf` | [M7] Override rsID source. Default = cached ClinVar VCF (INFO/RS) | `None` |
| `--somatic` | [M7] Enable COSMIC overlay (requires `--cosmic-vcf`) | `False` |
| `--cosmic-vcf` | [M7] Path to COSMIC-format VCF (registration required) | `None` |
| `--cosmic-inject-density` | [M7] Fraction of cohort sites overwritten with COSMIC records when `--somatic` | `0.005` |
| `--svs-per-person` | [M8] Number of SVs (DEL/DUP/INV) per person | `3` |
| `--sv-length-min` | [M8] Minimum SV length in bp (log-uniform draw) | `50` |
| `--sv-length-max` | [M8] Maximum SV length in bp | `10000` |
| `--error-rate` | [M9] Per-call probability of a GT flip (lightweight noise model) | `0.001` |
| `--dropout-rate` | [M9] Per-call probability of a coverage dropout (`./.:0:0:0,0`) | `0.0005` |
| `--art` | [M9, heavy] ART read simulation + `bcftools call`. Currently rejected; needs M11 reference FASTA | `False` |
| `--chromosomes` | [coalescent] List, range, or mix (e.g. `22`, `19,20,21,22`, `1-22`, `1-3,5,19-22,X`) | `22` |
| `--chr-length-mb` | [coalescent] Simulated prefix per chrom | `5.0` |
| `--demo-model` | [coalescent] stdpopsim model id; `none` for uniform | `OutOfAfrica_3G09` |
| `--population` | [coalescent] Sampling population | `CEU` |
| `--rec-rate` | [coalescent, `--demo-model=none`] Uniform recombination rate | `1e-8` |
| `--mu` | [coalescent, `--demo-model=none`] Mutation rate | `1.29e-8` |
| `--admixture` | Run M6 EUR + SAS + AFR → UK pulse, write per-person ancestry BED | `False` |
| `--eur-frac` | [admixture] EUR proportion | `0.60` |
| `--sas-frac` | [admixture] SAS proportion | `0.25` |
| `--afr-frac` | [admixture] AFR proportion (sum must be 1.0) | `0.15` |
| `--legacy-background` | Use M4 1000G-pool + power-law SFS sampler | `False` |
| `--background-glob` | [legacy] Source glob(s) for common variants | 1000G files in parent dir |
| `--n-background` | [legacy] Shared background site count | `500` |
| `--af-min` | [legacy] Minimum AF when loading the pool | `0.05` |
| `--sfs-alpha` | [legacy] Power-law exponent for the SFS | `2.0` |
| `--workers` | [perf] Worker processes for the parallel cohort BCF write (sample-slice, post-5e) and the per-person fan-out. `0` = auto (`os.cpu_count()`), `1` = serial. msprime simulation itself is single-threaded and runs serially across chromosomes regardless of `--workers`. Linux only. | `0` |
| `--cohort-mode` | [perf, Phase 5d.1] Cohort intermediate between simulation and BCF write: `sites_list` (in-RAM, n ≤ 30 000), `arrow` (materialised list + mmap'd Arrow IPC scratch; n ≥ 100 000 with parent peak ≲ 50% host RAM), `arrow-streaming` (parent streams directly from msprime tree to Arrow — required when predicted parent peak > 50% host RAM, e.g. WGS at n=3000+), or `auto`. The Arrow paths need `pip install pyarrow`. See TUTORIAL.md §9.1.1. | `auto` |
| `--cohort-arrow-batch-size` | [perf, Phase 5d.1] Sites per Arrow record batch when `--cohort-mode` resolves to `arrow`. Default identified empirically (parent peak RSS ≈ 9.5 bytes × batch_size × n_haplotypes per batch). Don't change unless characterising a specific host. | `256` |
| `--fanout-batch-size` | [perf] Persons grouped into one `bcftools query` invocation during the per-person VCF fan-out. Larger batches reduce subprocess overhead; smaller batches reduce the parent's batch RSS. Binding RAM ceiling is `(parent_baseline + B × per_person) × workers`. | `4` |
| `--mode` | [perf] Output shape: `per-person` (default), `cohort` (skip per-person fan-out), or `both`. All three flow through the streamed cohort pipeline — derive per-person VCFs later via `bcftools view -s` against the per-chrom cohort BCFs. | `per-person` |
| `--no-resume` | [perf] Ignore any existing `cohort.meta.json` + cohort BCFs and start a fresh simulation. Default behaviour resumes a prior run when its params match. | `False` |
| `--chr-chunk-mb` | [perf] Split each chromosome's msprime simulation into independent sub-chunks of this size (Mb). `0` = auto-pick from `psutil.virtual_memory().available` and `--workers`. Cross-chunk LD is lost (chunks simulate independently); short-range LD within chunks is preserved. | `0` (auto) |
| `--profile-memory` | [perf, diagnostic] Sample RSS to a TSV at this path, ~1 s cadence plus labelled marks at phase boundaries (sim start, chrom yielded, BCF written, etc.). Useful for OOM postmortems and peak-RSS plots. | `None` (off) |
| `--config` | Load values from this YAML config file (overrides built-in defaults; CLI flags still win over config values). If omitted, the tool auto-discovers `./generate_people_config.yaml` in the current directory. See TUTORIAL.md §10. | `None` |
| `--no-config` | Skip the auto-discovery of `generate_people_config.yaml` in the current directory. No effect when `--config` is set. | `False` |
| `--print-config` | Print a starter `generate_people_config.yaml` (every field at its built-in default, each documented with a leading `# description` comment) to stdout and exit. Redirect to a file to bootstrap a config: `generate_people --print-config > generate_people_config.yaml`. See TUTORIAL.md §10. | `False` |
