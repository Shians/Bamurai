# Bamurai

A Python toolkit for manipulating BAM and FASTQ files, designed to split reads into smaller fragments, extract statistics, validate files, and manage multi-sample data.

## Description

Bamurai is a command-line tool for splitting reads in BAM/FASTQ files into smaller fragments. It is designed to be fast and efficient, and can be used to split reads into a target length or a target number of pieces per read.

These are the current features of Bamurai:

1. Splitting reads in a file to a target length
2. Splitting reads in a file to a target number of pieces per read
3. Getting statistics from a BAM or FASTQ(.gz) file
4. Basic validation of BAM and FASTQ(.gz) files

The `split` command splits reads into a target length, each read will be split into fragments as close to the target length as possible. Reads shorter than the target length will not be split.

The `divide` command splits reads into a target number of pieces, each read will be split into the number of pieces specified. A further minimum length can be specified to ensure that reads are not split if the resultant fragments are less than the minimum length.

The `stats` command will output the following information by default:
```
Statistics for input.bam:
  Total reads: 8160
  Average read length: 30638
  Throughput: 250006998
  N50: 82547
```

It can be used with the `--tsv` argument to output the statistics in a tab-separated format for computational analysis.
```bash
file_name       total_reads     avg_read_len    throughput      n50
input.bam      8160    30638   250006998       82547
```

The `validate` command will check the integrity of a BAM or FASTQ(.gz) file and output the following information if the file is valid.:
```bash
input.bam is a valid BAM file with 8160 records.
```

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
bamurai split input.bam --target-length 10000 --output output.fastq
```

To create a gzipped output file
```bash
bamurai split input.bam --target-length 10000 | gzip > output.fastq.gz
```

### Dividing reads into a target number of pieces

To divide reads into 2 pieces
```bash
bamurai divide input.bam --num_fragments 2 --output output.fastq
```

To divide reads into 2 pieces unless resultant fragments are less than 1000 bp
```bash
bamurai divide input.bam --num_fragments 2 --min_length 1000 --output output.fastq
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

### Working with multi-sample BAM files

Bamurai provides commands for processing BAM files with multiple samples based on barcode information.

#### Splitting BAM or FASTQ files by donor ID

To split a BAM or FASTQ file into multiple files, one for each donor ID:

```bash
bamurai split_samples --input input.bam --tsv barcode_to_donor.tsv --output-dir donor_bams
bamurai split_samples --input input.fastq.gz --tsv barcode_to_donor.tsv --output-dir donor_fastqs
```

The TSV file should contain at least two columns with headers 'barcode' and 'donor_id'. Each row maps a barcode to a donor ID.

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

### Assigning samples to barcodes

The `assign_samples` command assigns donor IDs to barcodes based on a provided TSV mapping file. This is useful for annotating barcodes in single-cell data.

To assign donor IDs to barcodes in a TSV file:

```bash
bamurai assign_samples --barcodes barcodes.tsv --tsv barcode_to_donor.tsv --output assigned_barcodes.tsv
```

- `barcodes.tsv` should contain a list of barcodes (one per line or as a column in a table).
- `barcode_to_donor.tsv` should have at least two columns: 'barcode' and 'donor_id'.
- The output file `assigned_barcodes.tsv` will contain the barcodes with their assigned donor IDs.

You can also specify a custom column name for barcodes in the input file:

```bash
bamurai assign_samples --barcodes barcodes.tsv --tsv barcode_to_donor.tsv --barcode-column cell_barcode --output assigned_barcodes.tsv
```

This will use the column 'cell_barcode' in `barcodes.tsv` as the barcode column.

