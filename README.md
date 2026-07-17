# Bamurai

A Python toolkit for manipulating BAM and FASTQ files, designed to split reads into smaller fragments, extract statistics, validate files, and manage multi-sample data.

**For BAM/SAM/CRAM files, Bamurai only processes primary alignments. Secondary and supplementary alignments are ignored in all commands.** This approach ensures that each original read from sequencing is counted only once. Secondary and supplementary alignments represent alternative mappings or split alignments of the same read, not additional unique reads. Including them would artificially inflate read counts and statistics, leading to misleading results.

## Description

Bamurai is a command-line tool for splitting reads in BAM/FASTQ files into smaller fragments. It is designed to be fast and efficient, and can be used to split reads into a target length or a target number of pieces per read.

These are the current features of Bamurai:

1. Splitting reads in a file to a target length (`split`)
2. Splitting reads in a file to a target number of pieces per read (`divide`)
3. Getting statistics from a BAM or FASTQ(.gz) file (`stats`)
4. Basic validation of BAM and FASTQ(.gz) files (`validate`)
5. Splitting a file into size-based chunks (`chunk`)
6. Splitting or extracting reads by donor ID using a barcode mapping
   (`split_samples`, `extract_sample`, `assign_samples`)
7. Extracting hashtag oligo (HTO) information from 10x FASTQ pairs (`get_hto`)

The `split` command splits reads into a target length, each read will be split into fragments as close to the target length as possible. Reads shorter than the target length will not be split.

The `divide` command splits reads into a target number of pieces (default 2), each read will be split into the number of pieces specified. A further minimum length can be specified to ensure that reads are not split if the resultant fragments are less than the minimum length (default 100).

The `stats` command will output the following information by default:
```
Statistics for input.bam:
  Total reads: 8160
  Average read length: 30638
  Throughput (Gb): 0.25
  N50: 82547
```

It can be used with the `--tsv` argument to output the statistics in a tab-separated format for computational analysis. Note that `throughput` is reported in raw bases here, rather than the gigabases of the default output.
```bash
file_name       total_reads     avg_read_len    throughput      n50
input.bam      8160    30638   250006998       82547
```

The `validate` command will check the integrity of a BAM or FASTQ(.gz) file and output the following information if the file is valid.:
```bash
input.bam is a valid BAM file with 8160 records.
```

Unlike the read-processing commands, `validate` inspects every record in the file, so the count it reports includes any secondary and supplementary alignments.

## Installation

To install the released version of Bamurai from PyPI

```bash
pip install bamurai
```

To install the latest version of Bamurai from GitHub

```bash
pip install git+https://github.com/Shians/Bamurai.git
```

## Usage

To get help on the command-line interface and list available commands
```bash
bamurai --help
```

To get help on a specific command
```bash
bamurai <command> --help
```

### Splitting reads to target size

To split a file into 10,000 bp reads
```bash
bamurai split input.bam --len-target 10000 --output output.fastq
```

To create a gzipped output file
```bash
bamurai split input.bam --len-target 10000 | gzip > output.fastq.gz
```

### Dividing reads into a target number of pieces

To divide reads into 2 pieces
```bash
bamurai divide input.bam --num-fragments 2 --output output.fastq
```

To divide reads into 2 pieces unless resultant fragments are less than 1000 bp
```bash
bamurai divide input.bam --num-fragments 2 --min-length 1000 --output output.fastq
```

### Getting statistics from a BAM or FASTQ file

To get stats from a BAM file
```bash
bamurai stats input.bam
```

To get stats from a FASTQ file or Gzipped FASTQ file
```bash
bamurai stats input.fastq
bamurai stats input.fastq.gz
```

### Validating BAM or FASTQ files

To validate a BAM file
```bash
bamurai validate input.bam
```

### Splitting a file into chunks

To split a file into chunks of at least 1 GB each
```bash
bamurai chunk input.bam --size 1G
```

Sizes are given in human-readable notation (`1G`, `100M`, `1000K`). Each chunk is at least the requested size, and may be slightly larger because reads are never split across chunks.

Output files are written as `chunk_1.fastq`, `chunk_2.fastq`, and so on. Use `--prefix` to change the name
```bash
bamurai chunk input.bam --size 100M --prefix sample_a
```

### Working with multi-sample BAM files

Bamurai provides commands for processing BAM files with multiple samples based on barcode information.

#### Splitting BAM or FASTQ files by donor ID

To split a BAM or FASTQ file into multiple files, one for each donor ID:

```bash
bamurai split_samples --input input.bam --tsv barcode_to_donor.tsv --output-dir donor_bams
bamurai split_samples --input input.fastq.gz --tsv barcode_to_donor.tsv --output-dir donor_fastqs
```

The TSV file should contain a barcode column and a donor ID column, one row per barcode; see [Barcode and donor columns](#barcode-and-donor-columns) below for how those columns are found.

You can process multiple BAM or FASTQ files at once:

```bash
bamurai split_samples --input input1.bam input2.bam --tsv barcode_to_donor.tsv --output-dir donor_bams
bamurai split_samples --input input1.fastq.gz input2.fastq.gz --tsv barcode_to_donor.tsv --output-dir donor_fastqs
```

#### Extracting reads for a specific donor

To extract all reads belonging to a specific donor from a BAM file:

```bash
bamurai extract_sample --bam input.bam --tsv barcode_to_donor.tsv --donor-id donor1 --output donor1.bam
```

You can also process multiple BAM files at once, combining all donor-specific reads into a single output file:

```bash
bamurai extract_sample --bam input1.bam input2.bam input3.bam --tsv barcode_to_donor.tsv --donor-id donor1 --output donor1.bam
```

This command will extract all reads with barcodes belonging to the specified donor ID and write them to a new BAM file.

#### Tagging reads with their donor ID

The `assign_samples` command labels each read in a BAM with the donor it came from, using the `RG` (read group) tag. Unlike `split_samples`, which writes one file per donor, this keeps every read in a single BAM so that downstream tools can group by read group.

```bash
bamurai assign_samples --bam input.bam --tsv barcode_to_donor.tsv --output assigned.bam
```

This writes `assigned.bam` containing every read from the input, where:

- reads whose cell barcode maps to a donor carry an `RG` tag naming that donor
- reads with an unknown or absent barcode are written through untagged
- the header gains an `@RG` record per donor, with `ID` and `SM` both set to the donor ID

A summary of how many reads were assigned to each donor is printed on completion. Note that `--bam` takes a single BAM, unlike `split_samples` and `extract_sample`, which accept several.

#### Barcode and donor columns

`split_samples`, `extract_sample` and `assign_samples` all read the same barcode-to-donor TSV, and all detect its columns the same way. The barcode column is taken from a column named `barcode` or `cell`, and the donor column from one named `donor_id`. If both `barcode` and `cell` are present the file is ambiguous and the command fails, asking you to choose.

Use `--barcode-column` and `--donor-id-column` when your columns are named something else, or to resolve that ambiguity:

```bash
bamurai extract_sample --bam input.bam --tsv mapping.tsv --donor-id donor1 \
  --barcode-column cell_barcode --donor-id-column sample --output donor1.bam
```

#### Which barcode tag is read

For BAM input, a read's cell barcode is read from the first of these tags that is present:

| Tag | Meaning |
| --- | --- |
| `CB` | Cell barcode, error-corrected (preferred) |
| `CR` | Cell barcode, raw/uncorrected |
| `XC` | Cell barcode, Drop-seq and early 10x convention (not in the SAM spec) |

`BC` and `RX` are deliberately not consulted: per the SAM specification `BC` identifies a *sample* and `RX` holds the UMI, so neither names the cell that a donor is assigned to. A read carrying only those tags is treated as having no barcode.

### Extracting HTO information

The `get_hto` command pulls hashtag oligo (HTO) information out of a pair of 10x FASTQ files. It assumes R1 holds the cell barcode followed by the UMI, and R2 holds the hashtag sequence.

```bash
bamurai get_hto --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
  --bc-len 16 --umi-len 12 --output htos.tsv
```

`--bc-len` and `--umi-len` give the barcode and UMI lengths in R1, and are required. The output is a tab-separated file with the columns `read_name`, `cell_barcode`, `umi`, `hto`, `bc_qual`, `umi_qual` and `hto_qual`, where the three `_qual` columns carry the FASTQ quality strings for the preceding sequences.

Where in R2 the hashtag is read from can be adjusted with `--hashtag-len` (how many bases to take, default 15) and `--hashtag-left-buffer` (how many bases to skip first, default 10).

