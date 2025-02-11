# Bamurai

A Python package for splitting reads BAM files into smaller fragments.

## Description

Bamurai is a command-line tool that splits reads in BAM files to a rough target length and outputs a FASTQ file. This is useful for testing the performance of bioinformatics tools on long-reads data at different read lengths.

## Installation

```bash
pip install git+https://github.com/Shians/Bamurai.git
```

## Usage

To split a file into 10,000 bp reads
```bash
bamurai split input.bam --target-length 10000 --output output.fastq
```

To create a gzipped output file
```bash
bamurai split input.bam --target-length 10000 | gzip > output.fastq.gz
```

To get stats from a BAM file
```bash
bamurai stats input.bam
```

To get stats from a FASTQ file or Gzipped FASTQ file
```bash
bamurai stats input.fastq
bamurai stats input.fastq.gz
```