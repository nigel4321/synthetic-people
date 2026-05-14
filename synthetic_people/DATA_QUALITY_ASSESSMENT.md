# Data-quality assessment & roadmap

**Date:** 2026-05-11
**Base commit:** `ea89da5` (post-Phase-5d.1 streaming-cohort merges)
**Scope:** scientific fidelity of generated VCFs against the spec in
[`SYHTHETIC_PROJECT.md`](SYHTHETIC_PROJECT.md) and against real human WGS data.
**Out of scope:** code quality, performance, CI/CD, packaging.

This document records a deliberate review of where the generator's
output diverges from real genomes, organised by what would matter to
a scientist comparing synthetic batches against 1000 Genomes / gnomAD
or using them to benchmark a variant caller. The Phase 5d.1
streaming-cohort work has removed the memory ceiling that motivated
several earlier simplifications; this document maps which of those
can now be revisited.

The assessment is intentionally critical. The tool already meets the
structural spec (valid VCF 4.2, ClinVar/dbSNP/COSMIC overlays at real
coordinates, validated SFS and Ti/Tv, admixture truth tracks); the
gaps below are about the next layer of realism, not the foundation.

---

## 1. Scorecard: spec vs current state

| Spec requirement (§) | Status | Comment |
|---|---|---|
| VCF 4.2, contig headers for 1–22, X, Y, MT (§2.1) | ⚠ partial | Contig lengths exist in `syntheticgen/builds.py:6-24`; default is `--chromosomes 22`. X/Y/MT can be **declared** but not **simulated correctly** (see §3 below). |
| Phased GT, DP, GQ, AD (§2.2) | ✓ | `syntheticgen/quality.py` draws DP ~ Poisson(λ≈30) + per-sample jitter; AD carries the 0.475 het-alt ref-bias matching Illumina + BWA-MEM. |
| SNPs, indels, SVs, multi-allelic (§2.2) | ⚠ partial | SNPs ✓; SVs ✓ (DEL/DUP/INV); **indels only via overlays, none from coalescent**; **no true multi-allelics** (`BinaryMutationModel`). |
| LD via coalescent or HMM (§3.1) | ⚠ partial | msprime ✓, but with a **flat per-genome recombination rate** (no deCODE / HapMap recombination map). Chunked simulation explicitly drops cross-chunk LD. |
| Ti/Tv ≈ 2.1, SFS power-law (§3.2) | ✓ | Validated end-to-end in `syntheticgen/validate.py`. |
| Database grounding: ClinVar / dbSNP / COSMIC (§3.2) | ✓ structurally / ⚠ semantically | Injected at real coordinates with real `INFO/RS` and `CLNSIG`, **but injected genotypes do not respect the source DB's known population AF** (see §4.2 below). |
| UK admixture EUR + SAS + AFR (§4) | ⚠ partial | Single-pulse only (`syntheticgen/admixture.py:37`, `PULSE_TIME = 20.0` generations); real UK demography is multi-pulse + continuous post-WWII migration. |
| Local-ancestry BED truth (§4.2) | ✓ | Per-person ancestry BEDs from msprime `record_migrations`. |
| Validation suite: PCA vs 1000G, LD decay, Ti/Tv (§6) | ⚠ partial | PCA is computed on the synthetic cohort alone — **no projection against real 1000G reference samples**. LD decay ✓. Ti/Tv ✓. |
| Reference build aligned to GRCh38 (§2.1) | ✓ (since M12, 2026-05-14) | Pass `--reference-fasta` / `cohort.reference_fasta` and REF bases come from the FASTA via `pysam.FastaFile`. Legacy fabricated-REF path retained without the flag, preserving seed-pinned test behaviour. |
| Error modeling via ART / SimNGS (§5) | ⚠ substitute | Lightweight GT-flip + dropout model in `syntheticgen/errors.py` — context-free. ART would produce reads + recall variants — a different (much heavier) workflow. |

---

## 2. The headline

Today's `synthetic_people` is structurally what the spec asks for —
valid VCFs at real coordinates, real `CLNSIG`/`RS` annotations,
validated SFS and Ti/Tv, admixture truth tracks. **It is not yet
what a scientist reviewing it for benchmarking or population
genetics research would call high-fidelity**, primarily because of:

1. ~~fabricated `REF` bases (no FASTA loaded)~~ — **closed by M12
   (2026-05-14)**;
2. sex-chromosome ploidy unmodeled (X/Y/MT treated as autosomes);
3. overlay genotype/AF inconsistency (ClinVar variants land at
   arbitrary cohort AF);
4. uniform genome-wide recombination and mutation rates (no
   hotspots, no CpG enrichment).

M12 closes #1; the remaining three are tractable now that Phase
5d.1 has removed the memory ceiling that motivated them. **M13
(sex chromosomes)** is the next-largest perception gap with real
WGS, followed by M14 (context-aware μ) which is now unblocked by
M12's real REF.

---

## 3. The X / Y / MT problem — most consequential gap

The codebase has every contig length in `syntheticgen/builds.py`,
and the CLI accepts `--chromosomes X,Y,MT`, but **every simulation
site assumes `n_haplotypes = 2 × n_people` unconditionally**
(`syntheticgen/coalescent.py:398`, `syntheticgen/admixture.py:158`).
Consequences:

- **Males are simulated with two X haplotypes** (impossible) and a
  **diploid Y** (impossible).
- **Y is biallelic** rather than haploid — the non-recombining
  haploid Y is the cleanest patrilineal marker in real data and is
  the single most distinctive sex-chromosome property.
- **PAR1 (~2.78 Mb) and PAR2 (~330 kb)** — regions where X and Y do
  recombine in males and must carry identical variants between
  X-PAR and Y-PAR — are not modeled. Whole regions of the
  pseudo-autosomal genome are silently autosomal in the output.
- **MT** (clonal, maternally inherited, with heteroplasmy) is
  treated as standard diploid — wrong inheritance, wrong allele
  dosage, wrong copy number.
- **No `--sex` flag exists**, so the cohort is implicitly "all
  pseudo-2N individuals". A scientist running case/control on an
  X-linked trait would get nonsense.

The validation suite is silent on every one of these: no Y-haploid
check, no PAR-consistency check, no MT-clonality check, no X
male/female ploidy check. This is the spec-required deliverable
most cleanly absent.

---

## 4. Variant content: where the cohort diverges from a real cohort

### 4.1 ~~Reference base is fabricated~~ — closed by M12 (2026-05-14)

**Resolved** by M12. Pass `--reference-fasta <path>` (or
`cohort.reference_fasta` in YAML) and the cli looks up REF
from the FASTA via `pysam.FastaFile`. Without the flag the
legacy `rng.choice("ACGT")` path still applies — preserves
seed-pinned tests + lets quick smoke runs skip the FASTA.

The wiring touches four REF-picking sites (materialised path,
streaming pass-1 meta, streaming pass-2, admixture inline);
all share a `_pick_ref` helper that consumes one rng draw
unconditionally so the seed stream is invariant between
fasta-on and fasta-off paths (only the REF/ALT content of
the variant differs, not the overall rng trajectory).

Tier 1's REF-check gate (which previously failed on every
record by design) now passes end-to-end on a real run with
`--reference-fasta`. Empirical proof point captured in
`tests/test_reference.py::ReferenceEndToEndTest`.

### 4.2 Overlay AF inconsistency

`syntheticgen/clinvar.py:138` injects a ClinVar pathogenic variant
by overwriting `pos`/`ref`/`alt`/`id`/`INFO` on a randomly-chosen
cohort site **while keeping that site's coalescent-drawn genotype
block**. Effect: a ClinVar variant with a real population frequency
of, say, 1×10⁻⁴ can land at an arbitrary cohort AF (often 5–30%).
Pipelines that compare per-variant AF against gnomAD or ClinVar's
own population data will see the synthetic batches as internally
inconsistent.

The fix isn't a lookup-replace on the genotype block (that would
destroy LD), but to **rejection-sample the cohort site** whose AF
already matches the source DB's known frequency band, and to flag
the matched AF range in the per-batch manifest.

### 4.3 No indels from the coalescent

`BinaryMutationModel()` at `syntheticgen/coalescent.py:8` produces
SNPs only. Real WGS sits at ~10–15% indels by variant count. The
current output's indel content is 100% overlay-driven, capped at
whatever density the overlay flags allow. msprime's
`MutationModel` accepts a Jukes-Cantor or infinite-sites
finite-alphabet kernel that produces multi-allelics and indels at
similar runtime — binary was a complexity reduction, not a
performance one.

### 4.4 No CpG / context-dependent mutation rate

`DEFAULT_MU = 1.29e-8` uniformly across the genome
(`syntheticgen/coalescent.py:28`). Real CpG sites mutate ~10–15×
faster and are transition-biased. Combined with the fabricated
`REF`, there is no way to even retrofit this without first loading
the reference. Once the reference is loaded, a per-trinucleotide-
context Poisson rate from the published mutation-rate tables
(Karczewski 2020 / Roulette 2023) is one Python file.

### 4.5 Structural variants are IID, not segregating

`syntheticgen/sv.py:89` (`generate_person_svs`) draws SVs uniformly
per-person with no parent/offspring relationship — fine in a
coalescent IID cohort but wrong for any family/trio study. There
is no `--trio` or pedigree feature. SV breakpoints are also drawn
from a uniform position prior, not enriched at low-copy repeats
where >70% of real human DEL/DUP breakpoints land.

### 4.6 Sequencing-error model is context-free

`syntheticgen/errors.py:34-50` has uniform GT-flip rate and uniform
dropout. Real WGS error rates vary 10–100× by sequence context
(homopolymers, GC extremes, segmental duplications). For
benchmarking variant callers — one of the stated use cases in
spec §1 — flat error rates will *over-estimate* caller performance
because callers' filter sets are tuned for context-dependent errors.

### 4.7 Coverage is regionally uniform

`syntheticgen/quality.py:30-37` draws DP ~ Poisson(λ ≈ 30) globally.
Real WGS shows ~3× DP variation by region (centromeres, segmental
duplications, GC-extreme exons). A caller validated on this data is
being validated on the easy half of the genome.

---

## 5. Population genetics: the demographic model

- **Single-pulse admixture** (`syntheticgen/admixture.py:37`,
  `PULSE_TIME = 20`) is the deepest demographic simplification.
  Real UK admixed populations have at minimum two pulses
  (post-WWII Caribbean + 1960s–2000s South Asian) plus continuous
  gene flow. A scientist studying local-ancestry tract lengths
  would get the wrong tract-length distribution — the exponential
  decay from a single pulse is the cleanest detectable signal, and
  it's the wrong one.
- **No selection.** `OutOfAfrica_3G09` is neutral. No purifying
  selection at conserved sites, no balancing at HLA, no positive
  at LCT / EDAR. This shapes the SFS at functional sites.
- **No population substructure within continental groups.** "CEU"
  is one node. Real cohorts show fine-grained substructure (PoBI:
  ~17 UK clusters within "British").

stdpopsim — already a dependency — has both selection (DFE
catalogues) and per-population recombination maps available; the
code just doesn't invoke them. The cost is plumbing flags through,
not new science.

---

## 6. Performance compromises that can now be revisited

Phase 5d.1 (streaming-cohort + Arrow intermediate + chunked
simulation) collectively unblocked WGS-scale on 32 GB hosts. The
memory ceiling that drove several earlier simplifications no
longer binds:

| Compromise | Why it was made | Why it can now be revisited |
|---|---|---|
| Fabricated `REF` bases | Loading ~3 GB FASTA was significant when parent peak RSS was ~17 GB at n=3000 × 70 Mb | Streaming + Arrow cap parent RSS well below the prior ceiling; FASTA mmap (`pysam.FastaFile`) adds ~50 MB resident regardless of cohort size |
| `BinaryMutationModel` | Simpler downstream code (no multi-allelic GT logic) | Cohort writer already handles `INFO/AC,AN,AF` arrays; Jukes-Cantor / K80 is a 1-day extension |
| Cross-chunk LD loss | Full chr1 msprime exceeded 32 GB at WGS-n3000 | The streaming generator yields per-chunk site dicts; whole-chrom simulation is feasible if we accept slower wall time. Worth offering an explicit `--no-chunking` opt-in for accuracy-first runs. |
| No real reference / no CpG-aware μ | Required the FASTA | Same unblock — FASTA loading enables both |
| Uniform recombination rate | Avoided loading deCODE / HapMap maps | stdpopsim ships these per model; one keyword (`use_recombination_map=True`) per `simulate` call |
| No selection | Neutral-only was simpler | stdpopsim has DFE catalogues per population |
| Lightweight error model | ART/SimNGS need real reads | Keep the lightweight model as default but make it **context-aware** by sampling kmer / homopolymer / GC features from the now-loaded reference |

### 6.1 Empirical scaling ceiling for `arrow-streaming` (2026-05-12)

A user-driven run at **n=1 000 000**, full WGS (chromosomes 1-22 + X
@ 70 Mb per chrom), `--cohort-mode arrow-streaming`, `--mode
per-person`, `--workers 12` was killed by the kernel at **63 GB
parent RSS** roughly 115 minutes into chrom 1. memprof phase marks
showed the kill happened *during chrom 1's streaming-write pass*,
after `chrom 1 ts ready (1 081 579 sites)` at 2.6 GB and before any
worker fan-out (`children_rss_mb` stayed at 0 throughout).

Root cause (see `coalescent.py:611-619` for the architectural
budget the design bet against): the streaming heap in
`_stream_cohort_pass2` holds **full site dicts** whose `carriers`
field is a Python `list[tuple[int, int]]` — one tuple per non-zero
haplotype. At n=1M:

- A singleton (AF ≈ 1/2M) costs ~80 bytes per site.
- A common-AF site (AF ≈ 0.3) costs ~600 K tuples × ~80 bytes ≈
  **~50 MB per site**.

The buffer depth is bounded by overlay-injection-position pressure
(roughly `O(sqrt(N_inject))` if positions are uniform, larger when
they cluster). At the canonical 0.2 rsid density × real dbSNP
positions on chr1, the buffer plausibly holds hundreds of
common-AF sites concurrently — pushing RAM well above the host
ceiling.

**The streaming guarantee that Phase 5d.1 validated at n=3000
extrapolated cleanly to n=10 000; it does *not* extrapolate
linearly to n=1M**, because the per-site `carriers` payload
itself scales linearly with n. The fix is not architectural —
it's representational (pack `carriers` into a numpy array, ~10–20×
smaller) — and is captured as deferred work in
`PERFORMANCE_BUDGETS.md` § "Known scaling ceiling."

Today's supported envelope, updated:

- **n ≤ ~100 000, full WGS, `--mode cohort`** — comfortable.
- **n ≤ ~500 000, full WGS, `--mode cohort`** — feasible with
  reduced overlay densities (rsid ≤ 0.05) on a 64 GB host.
- **n ≥ ~1 000 000** — needs the carriers-packing fix (or a
  user-side workaround: `chr_length_mb` ≤ 10, overlay density 0,
  `--mode cohort`).

---

## 7. Proposal — prioritised roadmap

Ordered by **scientific impact per implementation cost**, not
chronologically.

### ~~M12 — Reference-aware foundation~~ — **shipped 2026-05-14**

- ~~Load GRCh38 primary FASTA via `pysam.FastaFile`~~ — done.
  New `syntheticgen/reference.py` module wraps loading +
  validation; `--reference-fasta` cli flag added (also exposed
  as `cohort.reference_fasta` in the YAML config).
- ~~Replace `rng.choice("ACGT")` with FASTA lookup~~ — done.
  All four REF-picking sites (`_tree_sequence_to_sites`,
  `_tree_sequence_to_sites_meta`, `_stream_cohort_pass2`, and
  the admixture inline producer) now go through a shared
  `_pick_ref` helper that prefers the FASTA when provided,
  falls back to `rng.choice` otherwise. The fallback path
  preserves the legacy seed-stream behaviour for tests and
  development runs that don't have a real FASTA.
- ~~Validation gate~~ — Tier 1's REF-check passes end-to-end
  when M12 wiring is correct. New
  `tests/test_reference.py::ReferenceEndToEndTest` runs the
  cli with a tiny synthetic FASTA and verifies every emitted
  REF matches the FASTA at that POS.
- **Rng-stream invariance**: `_pick_ref` always consumes one
  `rng.choice("ACGT")` draw regardless of whether the FASTA
  is used — so downstream rng consumers (overlay sampling,
  error model, AC/AN-driven calls) see the same rng state in
  both paths. Only the REF/ALT *content* of the variant changes
  between the two paths; the simulation's overall rng
  trajectory is identical.
- **What this unblocks**: M14 (CpG-aware μ — needs real
  trinucleotide context) and Tier 2 #5 (mutation spectrum —
  needs real REF to bin into the 96 context channels). Both
  are now mechanically actionable.

### M13 — Sex chromosomes & MT

- New `--sex` flag (`m`, `f`, or per-person draws).
- Per-person ploidy table; X non-PAR is haploid in males, Y is
  haploid in males and absent in females, MT is haploid and
  clonally inherited from a maternal-line sample.
- PAR1 / PAR2 simulated as a single template and copied to both X
  and Y in males.
- New validation gates: Y heterozygosity ≈ 0 in males, female Y
  absence, MT GT homogeneity.

### M14 — Realistic mutation & recombination

- Switch `BinaryMutationModel` → `JC69MutationModel` (or `K80`),
  keep Ti/Tv target.
- Per-trinucleotide-context μ from Karczewski / Roulette tables
  (now possible because `REF` is real).
- Pass `genetic_map="HapMapII_GRCh38"` to stdpopsim where the
  model supports it.

### M15 — Overlay AF consistency

- Index ClinVar / dbSNP records by population-AF band (rare, low,
  common).
- Rejection-sample cohort sites that match the source DB's known
  band before overwriting `pos`/`ref`/`alt`.
- Record realised vs target AF correlation in `manifest.json`.

### M16 — Demographic richness *(optional, use-case-driven)*

- Multi-pulse UK admixture: two-pulse default + an
  `--admixture-pulses` config block.
- Per-population sub-structure (PoBI clusters) as
  `--population british:cornwall` etc.
- `--selection-dfe` keyword to enable stdpopsim DFEs.

### M17 — Validation: against real 1000G

- The spec asks for PCA against 1000G samples — currently we do
  PCA on the synthetic cohort alone.
- Project synthetic samples onto pre-computed 1000G principal-
  component axes (chr19–22 phase3 VCFs are already available
  locally per `CLAUDE.md`).
- Acceptance gate: synthetic AFR/SAS/EUR samples land within 1σ
  of the corresponding 1000G clusters; admixed samples bridge them.

### M18 — Trios & pedigrees *(if family-aware tools are a use case)*

- Significant architectural lift: today every person is an IID
  coalescent draw on the shared backbone. Trio simulation needs a
  different generator that does parent → offspring meiosis.
- Defer unless requested — most variant-caller benchmarking
  doesn't need it.

---

## 8. What I would NOT change

- **Streaming-cohort architecture.** It's the right scaling
  primitive and survives every M12–M17 change unchanged.
- **Lightweight error model as default.** The right choice for
  fast iteration; ART would 10× the runtime and make benchmarking
  unwieldy. Keep it, just make it context-aware.
- **Per-person VCFs derived from a shared cohort backbone.** This
  is what makes M5+ tractable at WGS-n3000+. Trio support (M18) is
  a parallel generator, not a replacement.

---

## 9. References

- Karczewski et al. 2020. *The mutational constraint spectrum
  quantified from variation in 141,456 humans.* Nature 581.
- Roulette: <https://www.nature.com/articles/s41588-023-01573-x>
- 1000 Genomes Project Consortium 2015. *A global reference for
  human genetic variation.* Nature 526.
- stdpopsim catalogue: <https://popsim-consortium.github.io/stdpopsim-docs/>
- HapMap II / deCODE recombination maps:
  <https://github.com/popsim-consortium/stdpopsim/tree/main/stdpopsim/catalog/HomSap>

---

## Appendix A: Validation-coverage audit (2026-05-12)

A walk through `validate_batch.py` + `syntheticgen/validate.py`,
itemising what the existing acceptance suite actually proves and
where the silent gaps are. Companion to §1 — that scorecard is
about the *output content*; this appendix is about what the
*validator can detect*. The two are connected: if a fidelity gap
isn't visible to the validator, a regression in that area can
ship silently.

### A.1 What `validate_batch.py` checks today

Per-person walk (`summarise_vcf`, `syntheticgen/validate.py:183-231`)
produces a `SampleStats` per VCF, then aggregates across the cohort:

| Check | Source | Notes |
|---|---|---|
| Per-person record count | `n_records` | trivial sanity gate |
| SNV / indel / SV classification | `_classify_record` | by ALT shape (`<>` → SV, single-base → SNV, else indel) |
| Ti/Tv ratio | `titv_from_stats` | cohort-aggregate; Markdown report flags "outside [1.7, 2.6]" — wide band |
| Het / Hom-alt ratio | `het_hom_ratio` | per-sample + cohort; reported but uninterpreted |
| Per-record AF histogram | `af_histogram` | uses single-sample `INFO/AF`, 20 linear bins on [0, 1] |
| Indel length distribution | `aggregate_indel_lengths` | bp = `len(ALT) − len(REF)` |
| SV-by-type counts | `aggregate_sv_summary` | DEL/DUP/INV tallies |
| Singleton count | `s.singletons` | `INFO/AC == 1` in single-sample VCF — not a true cohort singleton |
| Dropout count | `n_dropout` | GT contains `.` |
| LD decay r² vs distance | `ld_decay` | log-spaced bins, ~5K pair sample per bin, single curve |
| Cohort PCA (synthetic-only) | `cohort_pca` | sklearn PCA on `(n_samples, n_variants)` dosage matrix |
| Admixture-mode PCA labels | `_default_pca_labels` | dominant-ancestry colouring from manifest |

Artefacts: `summary.json` + Markdown `report.md` + four PNGs
(`ld_decay`, `af_histogram`, `indel_lengths`, `pca`).

### A.2 What the validator does NOT check

Two categories: things the existing data *could* be checked
against today, and things blocked on infrastructure that doesn't
exist yet.

**Category 1 — checkable from current output, just not done:**

| Missing check | Fidelity gap it would catch |
|---|---|
| **REF allele matches GRCh38 at POS** | every fabricated REF (§4.1) — `bcftools norm --check-ref` would fail today on every record |
| **GT phasing preserved across overlay injections** | silent phase loss in `inject_clinvar` / `inject_rsids` / `inject_cosmic` |
| **Mutation spectrum (96-channel SNV context)** | the "no CpG μ" gap from §4.4 — real WGS is dominated by C>T at CpG |
| **Population-stratified AFs** | demographic-model misuse — the validator currently collapses to one univariate histogram |
| **Hardy-Weinberg equilibrium per site** | overlay-genotype mishandling, downstream of injection bugs |
| **F-statistic (inbreeding coefficient)** | hidden inbreeding in the simulation; admixture-mode F patterns |
| **Per-region variant density (per-Mb)** | flat density from uniform μ (§4.4) |
| **DP/GQ/AD distribution sanity** | silent regression in `quality.py`'s Poisson(λ≈30) + 0.475 ref-bias model |
| **Realised vs requested overlay density** | density-target drift in `--rsid-density` / `--clinvar-inject-density` / `--cosmic-inject-density` |
| **Per-chromosome statistics** | chrom-specific regressions (e.g. X-only bug after M13) invisible to today's cohort-wide aggregates |
| **Ti/Tv tolerance tightening** | the current `[1.7, 2.6]` band passes any biologically-plausible noise; real WGS is 2.0–2.1 ± 0.05 |
| **Realised admixture tract-length distribution** | ancestry BEDs are written but the validator never reads them |

**Category 2 — needs new infrastructure:**

| Missing check | Blocked on | Tied to roadmap |
|---|---|---|
| Sex-chromosome ploidy (Y haploid in males, X non-PAR haploid in males, MT clonal) | The simulation treats X/Y/MT as autosomes (§3) | **M13** |
| PAR1/PAR2 X↔Y consistency | Same | **M13** |
| PCA projected onto 1000G reference axes | Need cached 1000G phase3 VCFs + projection matrix | **M17** (spec §6.1 ask) |
| LD-block boundary fidelity at hotspots | Need real recombination map for comparison | **M14** |
| Mendelian consistency | No trio architecture today | **M18** |

### A.3 How the gaps connect to the roadmap

Mapping missing checks back to §7's M-milestones:

- **M12 (reference-aware FASTA)** unblocks: REF-matches-GRCh38 check.
- **M13 (sex chromosomes)** unblocks: Y-haploid, PAR consistency, MT clonality.
- **M14 (mutation + recombination)** unblocks: mutation spectrum, per-region density, LD-hotspot fidelity.
- **M15 (overlay AF consistency)** unblocks: HWE-per-site, realised-vs-target overlay density (the latter is actually checkable already without M15).
- **M17 (validation vs 1000G)** unblocks: PCA projection, population-stratified AF.
- **M18 (trios)** unblocks: Mendelian consistency.

The unblocked Category 1 checks (per-region density, per-chrom
stats, Ti/Tv tightening, rsID coverage, DP/GQ/AD sanity, mutation
spectrum modulo REF, phasing consistency, F-statistic, realised
admixture tract lengths) are all implementable today against the
current output. Adding them would 2–3× the validator's
discrimination power before any new code feature ships.

### A.4 Recommended additions, prioritised

Choose-your-own-adventure, ordered by catch-rate per cost.

#### Tier 1 — cheap, high-discrimination, no blockers

1. **REF-matches-GRCh38 gate**. Wrap a
   `bcftools norm --check-ref e -f reference.fa` in the validator
   (skip cleanly if no FASTA). Even without M12 in place, this
   proves the validator catches the bug if/when M12 lands
   incomplete. ~1 hour.
2. **Realised overlay density counters**. Walk the VCF `ID`
   column for rsIDs (for example, `ID` values starting with `rs`),
   plus `INFO/CLNSIG` and `INFO/COSMIC_ID`; tally non-empty
   fractions; compare against manifest's requested densities.
   ~1 hour.
3. **Per-chromosome breakouts**. Re-emit all the aggregate stats
   (Ti/Tv, het/hom, AF histogram, indel lengths) per chromosome,
   not just cohort-wide. ~2 hours.
4. **Ti/Tv tolerance tightening**. Drop the report band from
   `[1.7, 2.6]` to `[2.0, 2.2]` for WGS. The wider band hides
   drift. ~5 min + recalibration if it fails.

#### Tier 2 — moderately cheap, illuminates the model

5. **Mutation spectrum (96-channel)**. For each SNV with a real
   REF (post-M12), bin into the 96 trinucleotide contexts.
   Compare against COSMIC SBS1. Today's spectrum is degenerate;
   post-M14 it should match reality. **Blocked on M12** for real
   REF. ~1 day.
6. ~~**Per-region variant density**~~ — **shipped 2026-05-13**.
   1 Mb bins per chrom, with a coefficient-of-variation
   diagnostic. Today flat (CV ≈ 0); post-M14 expect 0.5–1.0 on
   most chroms.
7. ~~**DP/GQ/AD distribution sanity**~~ — **shipped 2026-05-13**.
   Sampled at ~50K records per VCF; summary stats (mean / median
   / stdev / p10 / p90) compared against the targets baked into
   `quality.py` (DP=30, AD ref-fraction=0.475 at hets).
8. ~~**F-statistic / inbreeding coefficient**~~ — **shipped
   2026-05-13**. Computed from the existing
   `build_genotype_matrix` cohort dosage matrix; expected F ≈ 0
   for outbred cohorts. |F| > 0.05 flags drift.
9. ~~**Realised admixture tract-length distribution**~~ —
   **shipped 2026-05-13**. Parses per-person ancestry BEDs;
   reports mean / median tract length per population. Activates
   only when ancestry BEDs are present (i.e., admixture mode).

#### Tier 3 — wait for the corresponding feature

10. Sex-chromosome ploidy checks — wait for M13.
11. PCA-vs-1000G projection — wait for M17.
12. Mendelian consistency — wait for M18.

### A.5 Recommendation

**Ship Tier 1 first, as one PR, before any M12+ code work.**
Reasoning:

- Total cost ~5 hours.
- Becomes the regression net for M12+. When M12 lands and the
  REF-check gate passes, that's empirical evidence the FASTA
  wiring actually works.
- Each Tier 1 check is independent — no architectural risk.
- Today's discrimination power is weak enough that some M12+
  features could silently regress and the suite would pass.

After Tier 1: revisit M12–M18 ordering with the sharper
validator in hand. The priority of M12 vs M15 vs M17 may shift
once we can actually *measure* overlay AF realism, mutation
spectrum, etc. — instead of intuiting it.
