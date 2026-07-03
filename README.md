# dORFannotator: functional annotation of variants affecting downstream open reading frames in 3′UTRs

dORFannotator is a freely available tool for annotating SNVs, MNVs and small indels that create or alter downstream open reading frames (dORFs) and downstream overlapping open reading frames (doORFs) in 3′ untranslated regions (3′UTRs). 

Reference data and precomputed annotation-ready resources associated with the dORFannotator study are available from Zenodo: <https://doi.org/10.5281/zenodo.21148673>.

---

## Table of Contents

- [dORFannotator: functional annotation of variants affecting downstream open reading frames in 3′UTRs](#dorfannotator-functional-annotation-of-variants-affecting-downstream-open-reading-frames-in-3utrs)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
    - [Prerequisites](#prerequisites)
    - [Install from source](#install-from-source)
  - [Usage](#usage)
    - [Build a reference database](#build-a-reference-database)
    - [Annotate TSV variants](#annotate-tsv-variants)
    - [Annotate VCF variants](#annotate-vcf-variants)
    - [Common options](#common-options)
  - [Examples](#examples)
    - [ClinVar example test](#clinvar-example-test)
    - [Use the precomputed database](#use-the-precomputed-database)
  - [Input and output formats](#input-and-output-formats)
    - [Variant input](#variant-input)
    - [Translated dORF evidence input](#translated-dorf-evidence-input)
    - [Annotation output](#annotation-output)
    - [Consequence terms](#consequence-terms)
  - [Contributing](#contributing)
  - [License](#license)
  - [Contact](#contact)

---

## Installation

### Prerequisites

- Python 3.10 or newer
- `pip`

### Install from source

```bash
# Clone the repository
git clone https://github.com/liusihan/dORFannotator.git

# Navigate into the dORFannotator directory
cd dORFannotator

# Install locally via pip
pip install -e .

# Test installation
dORFannotator -h
```

## Usage
dORFannotator has two main subcommands:

```bash
dORFannotator build
dORFannotator annotate
```

### Build a reference database

```bash
dORFannotator build \
  --gtf input.gtf.gz \
  --fasta hg38.fa \
  --translated-dorf translated_dORFs.tsv \
  --out-db output.db
```

### Annotate TSV variants

```bash
dORFannotator annotate \
  --tsv input.tsv.gz \
  --db output.db \
  --out annotated.dORFannotator.tsv
```

### Annotate VCF variants

```bash
dORFannotator annotate \
  --vcf input.vcf.gz \
  --db output.db \
  --out annotated.dORFannotator.tsv
```

### Common options

```bash
dORFannotator annotate \
  --vcf input.vcf.gz \
  --db output.db \
  --out annotated.dORFannotator.tsv \
  --chr chr1,chr2,chrX \
  --threads 8 \
  --batch-size 10000 \
  --mane-only
```

Selected options:

| Option | Description |
|---|---|
| `--vcf` | Input VCF or VCF.GZ file. |
| `--tsv` | Input tab-delimited variant file. |
| `--db` | dORFannotator SQLite reference database. |
| `--out` | Output TSV file. |
| `--chr` | Restrict annotation to one or more chromosomes. |
| `--region` | Restrict annotation to a genomic interval, for example `chr1:1000000-2000000`. |
| `--threads` | Number of CPU threads. |
| `--batch-size` | Number of variants processed per batch. |
| `--mane-only` | Report annotations only for MANE transcripts. |
| `--include-predicted` | Include all sequence-predicted ORF consequences. |
| `--evidence-only` | Report only consequences affecting translated-evidence ORFs. |

By default, `annotate` uses an evidence-first output mode: all consequences affecting translated-evidence ORFs are reported, while sequence-predicted ORFs are limited to strong-Kozak `dStart_gained`, strong-Kozak `dStop_gained`, and `dKozak_changed` events whose alternate Kozak strength is strong. Use `--include-predicted` to report all sequence-predicted ORF consequences.

---

## Examples

### ClinVar example test

This repository includes a ClinVar-derived input file and the expected dORFannotator output:

```text
tests/clinvar_mane_3utr.filtered.tsv.gz
tests/clinvar.dORFannotator.tsv
```

### Use the precomputed database

Download `gencode.v45.db.zip` from the Zenodo record (https://doi.org/10.5281/zenodo.21148673) and extract it:

```bash
mkdir -p db
unzip gencode.v45.db.zip -d db
```

Run the example annotation:

```bash
dORFannotator annotate \
  --tsv tests/clinvar_mane_3utr.filtered.tsv.gz \
  --db db/gencode.v45.db \
  --out tests/test.tsv \
  --include-predicted
```

## Input and output formats

### Variant input

TSV input must contain the following columns:

```text
chrom  pos  ref  alt
```

Header names are case-insensitive. `#CHROM`, `POS`, `REF`, and `ALT` are also accepted. If no recognized header is present, the first four columns are interpreted as `chrom`, `pos`, `ref`, and `alt`. Comma-separated ALT alleles are expanded.

VCF and VCF.GZ files are supported with `--vcf`.

### Translated dORF evidence input

The translated dORF evidence file used by `build` must contain six tab-delimited columns:

```text
chrom  start  end  strand  transcript_id  orf_class
```

Coordinates are 1-based genomic positions. `orf_class` must be `dORF` or `doORF`. Extra columns after `orf_class` are ignored.

### Annotation output

The standalone annotator writes a compact TSV with fixed columns:

```text
chrom
pos
ref
alt
gene
transcript
mane
strand
orf_class
csq
evidence
dorf_count
dorf_start
dorf_end
dist_cds
detail
```

Column descriptions:

| Column | Description |
|---|---|
| `chrom`, `pos`, `ref`, `alt` | Input variant allele. |
| `gene` | Gene symbol associated with the annotated transcript. |
| `transcript` | Transcript identifier. |
| `mane` | Whether the transcript is a MANE transcript. |
| `strand` | Transcript strand. |
| `orf_class` | `dORF` or `doORF`. |
| `csq` | dORFannotator consequence term. |
| `evidence` | `true` if the affected reference ORF matches translated dORF evidence; otherwise `false`. |
| `dorf_count` | Number of interpreted dORFs for the annotated transcript. |
| `dorf_start`, `dorf_end` | 1-based genomic coordinates of the interpreted ORF. |
| `dist_cds` | Transcript-oriented distance from the annotated CDS stop to the interpreted ORF start. |
| `detail` | Consequence-specific semicolon-separated `key=value` fields. |

The `evidence` field indicates support for the reference ORF only. It does not assert experimental support for a variant-created ORF.

### Consequence terms

dORFannotator reports the following primary consequence terms:

| Consequence | Description |
|---|---|
| `dStart_lost` | Variant disrupts a dORF/doORF initiation codon. |
| `dStart_changed` | Variant changes one initiation codon to another accepted initiation codon. |
| `dStart_gained` | Variant creates a new accepted initiation codon and a candidate dORF/doORF. |
| `dStop_lost` | Variant disrupts a dORF/doORF termination codon. |
| `dStop_changed` | Variant changes one stop codon to another stop codon. |
| `dStop_gained` | Variant creates a premature stop codon in a dORF/doORF. |
| `dFrameshift` | Indel changes the reading frame of a dORF/doORF. |
| `dInframe` | Indel changes a dORF/doORF sequence without changing frame. |
| `dMissense` | SNV changes the encoded amino acid in a dORF/doORF. |
| `dSynonymous` | SNV does not change the encoded amino acid in a dORF/doORF. |
| `dKozak_changed` | Variant changes Kozak-context strength around a dORF/doORF initiation codon. |

Kozak strength uses the canonical −3 A/G and +4 G rule:

| Strength | Definition |
|---|---|
| `strong` | Both positions match. |
| `moderate` | One position matches. |
| `weak` | Neither position matches. |

---

## Contributing

Contributions are welcome. Please submit pull requests or open issues on the GitHub repository.

---

## License

This project is licensed under the GNU General Public License v3.0. See `LICENSE`.

---

## Contact

- **Author**: Sihan Liu
- **Institution**: Institute of Rare Diseases, West China Hospital, Sichuan University
- **Email**: liusihan@wchscu.cn
- **Issues**: Please report bugs and feature requests on GitHub