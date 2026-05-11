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
| Reference build aligned to GRCh38 (§2.1) | ⚠ broken | Coordinates are GRCh38, but **`REF` bases are fabricated** uniformly from `{A,C,G,T}` (`syntheticgen/coalescent.py:440`). A real tool ingesting these VCFs cannot re-align them. |
| Error modeling via ART / SimNGS (§5) | ⚠ substitute | Lightweight GT-flip + dropout model in `syntheticgen/errors.py` — context-free. ART would produce reads + recall variants — a different (much heavier) workflow. |

---

## 2. The headline

Today's `synthetic_people` is structurally what the spec asks for —
valid VCFs at real coordinates, real `CLNSIG`/`RS` annotations,
validated SFS and Ti/Tv, admixture truth tracks. **It is not yet
what a scientist reviewing it for benchmarking or population
genetics research would call high-fidelity**, primarily because of:

1. fabricated `REF` bases (no FASTA loaded);
2. sex-chromosome ploidy unmodeled (X/Y/MT treated as autosomes);
3. overlay genotype/AF inconsistency (ClinVar variants land at
   arbitrary cohort AF);
4. uniform genome-wide recombination and mutation rates (no
   hotspots, no CpG enrichment).

All four are tractable now that Phase 5d.1 has removed the memory
ceiling that motivated them. M12 (reference-aware foundation) and
M13 (sex chromosomes) would close the largest perception gap with
real WGS on their own.

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

### 4.1 Reference base is fabricated

`syntheticgen/coalescent.py:440-442` does `ref = rng.choice("ACGT")`.
The VCF is structurally valid but **its `REF` field does not match
GRCh38 at that position**. Anything that re-validates against the
reference — `bcftools norm --check-ref`, `vcfvalidator`, any aligner
— will fail. This is acknowledged in the `README.md` as "M11 will
wire in the FASTA" but is the current state of `main`.

The original memory rationale for skipping FASTA loading (~3 GB for
GRCh38 primary) no longer applies post-Phase-5d.1: `pysam.FastaFile`
mmaps the reference, adding ~50 MB resident regardless of cohort
size — orders of magnitude below the cohort-intermediate we already
build.

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

---

## 7. Proposal — prioritised roadmap

Ordered by **scientific impact per implementation cost**, not
chronologically.

### M12 — Reference-aware foundation *(unblocks several others)*

- Load GRCh38 primary FASTA via `pysam.FastaFile` (cached alongside
  ClinVar in `cache/`).
- Replace `rng.choice("ACGT")` with `fasta.fetch(chrom, pos-1, pos)`.
  ALT picking unchanged.
- Validation gate: round-trip every emitted VCF through
  `bcftools norm --check-ref e -f reference.fa` in CI on a small
  batch. Currently this would fail at 100%.

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
