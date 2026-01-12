# dORFannotator

**dORFannotator** is a high-performance bioinformatics tool for annotating 3' untranslated region (UTR) variants that affect downstream Open Reading Frames (dORFs). It identifies and classifies how genetic variants in the 3'UTR can create, destroy, or modify downstream ORFs, which may have functional consequences for gene regulation and protein expression.

Currently, dORFannotator will annotate whether a small variation (1-5bp) including SNVs and indels in 3'UTR would have any of the following molecular consequences:

- **dStart_gained**: Variant creates a new start codon in the 3'UTR (e.g. A→G mutation creates ATG from ATA)

- **dStop_gained**: Variant creates a new stop codon (e.g. G→T mutation creates TAA)

-  **dStart_lost**: Variant distrup an existing start codon (e.g. ATG→ATC)

-  **dStop_lost**: Variant distrup an existing stop codon (e.g. TAA→CAA)

-  **dStart_change**: Start codon changes to another valid start codon (ATG→GTG)

-  **dStop_change**: Stop codon changes to another valid stop codon (TAA→TAG)

-  **dMissense**: Single nucleotide variant within a dORF result same amino acid change

-  **dIndel_frameshift**: Insertion/deletion causing frameshift in dORF

-  **dIndel_inframe**: Insertion/deletion maintaining reading frame

## Installation

### Prerequisites

dORFannotator requires Python 3.11+ and the following dependencies:

- **pysam** >= 0.19.0
- **cyvcf2** >= 0.30.0
- **biopython** >= 1.79
- **pandas** >= 1.3.0
- **numpy** >= 1.20.0

### Install

```bash
# Install Dependencies
pip install pysam cyvcf2 biopython pandas numpy

# Clone the repository
git clone https://github.com/liusihan/dORFannotator.git
cd dORFannotator
python setup 

# install via pip
pip install dorfannotator

# The tool is a standalone Python script, no additional installation needed
# Just ensure dependencies are installed
```

## Usage
To run dORFannotator, you could the following command line:
```bash
python dORFannotator.py \
    --genome <reference.fa> \
    --vcf <variants.vcf.gz> \
    --dorf <known_dORFs.bed> \
    --annotation <annotation.gtf> \
    --out <output_prefix>
```
To get a full list of options use
```bash
python dORFannotator.py -h
```


### Required Arguments

| Argument | Short | Description |
|----------|-------|-------------|
| `--genome` | `-g` | Reference genome FASTA file (must have `.fai` index) |
| `--vcf` | `-v` | Input VCF file (must be bgzipped with `.tbi` index) |
| `--dorf` | `-d` | Known dORF BED file |
| `--annotation` | `-a` | Gene annotation GTF/GFF file |
| `--out` | `-o` | Output file prefix |

Format of known dORF file (BED):

- Standard BED format
- Minimum columns: `chrom`, `start`, `end`, `strand`, `transcript`
- Optional 6th column: dORF identifier

```
chr1    990561  990719  +       ENST00000379370.2
chr1    991234  991320  +       ENST00000379370.2
```

### Optional Arguments


| Argument | Description | Example |
|----------|-------------|---------|
| `--chr` | Filter by chromosome(s) | `--chr chr1`, `--chr 1` or `--chr chr1,chr2,chrX` |
| `--region` | Filter by genomic region | `--region chr1:1000000-2000000` |
| `--threads` | Number of threads | 1 |
| `--start-codons` | Comma-separated start codons | `ATG,CTG,GTG,TTG,ACG` |
| `--stop-codons` | Comma-separated stop codons | `TAG,TAA,TGA` |
| `--min-length` | Minimum ORF length (bp) | 30 |
| `--max-length` | Maximum ORF length (bp) | 303 |
| `--verbose` | Enable verbose logging | False |
| `--version` | Show version and exit | - |




## Output Format

The output is a tab-separated value (TSV) file named `<output_prefix>_dORFannotator.tsv` with the following columns:

- `chromosome`: Chromosome name
- `genomic_pos`: Genomic position (1-based)
- `ref`: Reference allele
- `alt`: Alternate allele
- `Symbol`: Gene symbol
- `transcript`: Transcript ID
- `strand`: Strand orientation (+ or -)
- `existing_dORF_count`: Number of known dORFs in the transcript's 3'UTR
- `affected_dORF_id`: ID of the affected dORF (if applicable)
- `variant_type`: One of the 8 variant classification types (see below)
- `dORF_type`: `predicted` (newly created) or `existing` (affects known dORF)
- `dORF_genomic_start`: Genomic start position (1-based)
- `dORF_genomic_end`: Genomic end position (1-based)
- `dORF_cdna_start`: cDNA start position (relative to transcript)
- `dORF_cdna_end`: cDNA end position (relative to transcript)
- `dORF_length`: dORF length in base pairs
- `distance_to_CDS_stop`: Distance from CDS stop codon to dORF start (bp)
- `overlap_CDS`: Whether dORF overlaps with CDS (`TRUE`/`FALSE`)
- `dORF_start_codon`: Start codon sequence (e.g., ATG)
- `dORF_stop_codon`: Stop codon sequence (e.g., TAA)
- `dORF_sequence`: Complete dORF nucleotide sequence
- `dORF_AA`: Translated amino acid sequence
- `dORF_GC_percent`: GC content percentage
- `kozak_sequence`: Kozak context sequence
- `kozak_strength`: Kozak strength classification (`Strong`, `Moderate`, `Weak`)
- `aa_change`: Amino acid change (for missense variants)
- `alt_stop_exists`: Whether alternative stop codon exists after variant
- `alt_stop_distance`: Distance to alternative stop (if exists)
- `indel_type`: Type of indel effect (`frameshift` or `inframe`)


## Examples
The directory **examples/** contains some small example files that are useful when getting started. A test run on a set of binary traits can be achieved by the following 2 commands.

```bash
cd examples
## samtools faidx test_reference.fa
## bgzip -c test_variants.vcf > test_variants.vcf.gz
## tabix -p vcf test_variants.vcf.gz

python ../dORFannotator.py \
    --genome test_reference.fa \
    --vcf test_variants.vcf.gz \
    --dorf test_known_dORF.bed \
    --annotation test_annotation.gtf \
    --out test_output
```
The output result from this command is included in example/test_output_dORFannotator.tsv.



## Citation

If you use dORFannotator in your research, please cite:

```
XXXXXX
```

## Contact

- **Author**: Sihan Liu
- **Institution**: Institute of Rare Diseases, West China Hospital, Sichuan University
- **Email**: liusihan@wchscu.cn
- **Issues**: Please report bugs and feature requests on GitHub

## License

This project is licensed under the GNU General Public License v3.0.
