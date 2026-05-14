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
    },
}
