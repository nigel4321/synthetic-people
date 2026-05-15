"""Reference build metadata: contig lengths, ClinVar URLs, reference FASTA."""

from __future__ import annotations


GRCH37_CONTIG_LENGTHS = {
    "1":  249250621, "2":  243199373, "3":  198022430, "4":  191154276,
    "5":  180915260, "6":  171115067, "7":  159138663, "8":  146364022,
    "9":  141213431, "10": 135534747, "11": 135006516, "12": 133851895,
    "13": 115169878, "14": 107349540, "15": 102531392, "16":  90354753,
    "17":  81195210, "18":  78077248, "19":  59128983, "20":  63025520,
    "21":  48129895, "22":  51304566, "X": 155270560, "Y":  59373566,
    "MT": 16569,
}

GRCH38_CONTIG_LENGTHS = {
    "1":  248956422, "2":  242193529, "3":  198295559, "4":  190214555,
    "5":  181538259, "6":  170805979, "7":  159345973, "8":  145138636,
    "9":  138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16":  90338345,
    "17":  83257441, "18":  80373285, "19":  58617616, "20":  64444167,
    "21":  46709983, "22":  50818468, "X": 156040895, "Y":  57227415,
    "MT": 16569,
}

# M13: pseudoautosomal region (PAR) coordinates per build. The two
# PARs (PAR1 at the chr telomere, PAR2 at the opposite end) are
# regions where chrX and chrY share IDENTICAL sequence and recombine
# in male meiosis just like autosomes — hence "pseudo-autosomal".
# Outside the PARs, chrX non-PAR is X-only haploid in males and
# diploid in females; chrY non-PAR is Y-only haploid in males and
# absent in females.
#
# Coordinates are 1-based, inclusive of both endpoints. Sources:
#   - GRCh37 / GRCh38: UCSC Table Browser ``par`` track + GRC
#     authoritative coordinates published with each assembly.
#   - On GRCh38, only PAR1 has matching bp ranges between chrX and
#     chrY (both 10,001-2,781,479). PAR2 sits near the telomere of
#     each chromosome and therefore lands at different bp on chrX
#     (~155.7-156.0 Mb, near the end of the ~156 Mb chrX) vs chrY
#     (~56.9-57.2 Mb, near the end of the ~57 Mb chrY).

GRCH37_PAR_REGIONS = {
    "X": [(60_001, 2_699_520), (154_931_044, 155_260_560)],
    "Y": [(10_001, 2_649_520), (59_034_050, 59_363_566)],
}

GRCH38_PAR_REGIONS = {
    "X": [(10_001, 2_781_479), (155_701_383, 156_030_895)],
    "Y": [(10_001, 2_781_479), (56_887_903, 57_217_415)],
}

BUILDS = {
    "GRCh37": {
        "assembly": "GRCh37",
        "clinvar_url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/"
                       "vcf_GRCh37/clinvar.vcf.gz",
        # Annotation URL for the VCF ``##reference=`` header line.
        # Identifies the assembly to downstream tooling without
        # implying the file was actually fetched from this URL.
        "reference": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
                     "GCA_000001405.14_GRCh37.p13/"
                     "GCA_000001405.14_GRCh37.p13_genomic.fna.gz",
        # M12 reference FASTA. Ensembl primary assembly is the
        # natural fit because its FASTA records are named ``1``,
        # ``2``, …, ``22``, ``X``, ``Y``, ``MT`` — matching the
        # ``contigs`` table above exactly. NCBI's full-assembly
        # FASTA uses RefSeq accession IDs (``CM000663.2`` etc.)
        # which would need a separate mapping table.
        "reference_fasta_url": (
            "https://ftp.ensembl.org/pub/grch37/release-113/"
            "fasta/homo_sapiens/dna/"
            "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
        ),
        "contigs": GRCH37_CONTIG_LENGTHS,
        "par_regions": GRCH37_PAR_REGIONS,
    },
    "GRCh38": {
        "assembly": "GRCh38",
        "clinvar_url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/"
                       "vcf_GRCh38/clinvar.vcf.gz",
        "reference": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
                     "GCA_000001405.15_GRCh38/"
                     "GCA_000001405.15_GRCh38_genomic.fna.gz",
        "reference_fasta_url": (
            "https://ftp.ensembl.org/pub/release-113/"
            "fasta/homo_sapiens/dna/"
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        ),
        "contigs": GRCH38_CONTIG_LENGTHS,
        "par_regions": GRCH38_PAR_REGIONS,
    },
}


# M13: per-chromosome ploidy lookup. Reads PAR regions from BUILDS
# so a future build addition (e.g. CHM13) only needs the data; the
# logic here doesn't change.

_AUTOSOMES = frozenset({str(i) for i in range(1, 23)})


def is_in_par(chrom: str, pos: int, build: str) -> bool:
    """True iff ``pos`` (1-based) on ``chrom`` falls inside a PAR.

    Only chrX / chrY have PARs; returns False for any other chrom.
    The PAR ranges in BUILDS are inclusive of both endpoints.
    """
    if chrom not in ("X", "Y"):
        return False
    par_regions = BUILDS.get(build, {}).get("par_regions", {})
    for lo, hi in par_regions.get(chrom, []):
        if lo <= pos <= hi:
            return True
    return False


def ploidy_for(
    chrom: str, sex: str, build: str = "GRCh38",
    pos: int | None = None,
) -> int:
    """Return the ploidy of ``chrom`` for a person of the given sex.

    Sex must be ``"m"`` or ``"f"``. ``pos`` is only consulted for
    chrX (in males, to distinguish PAR from non-PAR) and chrY (in
    males, same reason — chrY is absent in females regardless of
    position).

    Return values:

      * **2** — diploid (autosomes always; chrX in females always;
        chrX or chrY PAR positions in males when ``pos`` lands in
        a PAR range).
      * **1** — haploid (chrX non-PAR in males; chrY non-PAR in
        males; MT in everyone; chrX in males when ``pos`` is None
        because the non-PAR is the dominant case).
      * **0** — chromosome absent (chrY in females, at any position).

    Unknown chroms default to **2** (defensive — a user-supplied
    custom contig shouldn't crash the simulator).

    When ``pos`` is None for chrX/chrY, returns the **non-PAR**
    ploidy. Callers that care about PAR semantics MUST pass ``pos``.
    """
    if sex not in ("m", "f"):
        raise ValueError(
            f"ploidy_for: sex must be 'm' or 'f', got {sex!r}",
        )
    if chrom in _AUTOSOMES:
        return 2
    if chrom == "MT":
        return 1
    if chrom == "X":
        if sex == "f":
            return 2
        # Male X: PAR diploid, non-PAR haploid. Without a position
        # we conservatively return the non-PAR answer (1) — callers
        # that need PAR semantics MUST pass pos.
        if pos is not None and is_in_par("X", pos, build):
            return 2
        return 1
    if chrom == "Y":
        if sex == "f":
            return 0
        # Male Y: same PAR / non-PAR split as X.
        if pos is not None and is_in_par("Y", pos, build):
            return 2
        return 1
    # Unknown chrom — default to autosomal diploid to avoid breaking
    # the simulator on custom contigs.
    return 2
