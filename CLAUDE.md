# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Bamurai is a Python CLI tool for manipulating BAM and FASTQ files. It splits reads into smaller fragments, extracts statistics, validates files, and manages multi-sample data.

**Critical Design Principle**: For BAM/SAM/CRAM files, Bamurai only processes primary alignments. Secondary and supplementary alignments are ignored in all commands to avoid artificially inflating read counts.

## Development Commands

### Running Tests
```bash
# Run all tests
python -m pytest

# Run specific test file
python -m pytest tests/test_core.py

# Run a specific class or test (tests are organised into classes)
python -m pytest tests/test_core.py::TestSplitRead
python -m pytest tests/test_core.py::TestSplitRead::test_single_split

# Run with coverage (config lives in pyproject.toml)
python -m pytest --cov=bamurai --cov-report=term-missing

# Regenerate the static test inputs in tests/data/
# REQUIRED after changing tests/generate_test_data.py or any builder in
# tests/conftest.py -- otherwise the drift guard test fails. See "Testing".
python tests/generate_test_data.py
```

### Building and Installing
```bash
# Install in development mode
pip install -e .

# Install with dev dependencies
pip install -e ".[dev]"

# Build package
python -m build
```

### Running the CLI
```bash
# After installation, use the CLI directly
bamurai --help
bamurai split input.bam --len_target 10000 --output output.fastq
bamurai stats input.bam
bamurai divide input.bam --num_fragments 2 --output output.fastq
```

### Code Quality
```bash
# Format code
black bamurai/

# Type checking
mypy bamurai/
```

## Architecture

### Core Components

**bamurai/core.py**: Central parsing and read utilities
- `Read` dataclass: Represents a FASTQ read with validation
- `parse_reads()`: Unified generator for parsing BAM/FASTQ/FASTQ.gz files
  - Only processes primary alignments for BAM files (filters out secondary/supplementary)
  - Returns `Read` objects with standardized structure
- `split_read()`: Splits a read at specified positions, appends fragment indices to read IDs
  - Returns new `Read` objects; does not mutate the input read
- `qual_to_fastq_numpy()`: Converts pysam quality scores to FASTQ quality strings

**bamurai/cli.py**: Command-line interface with argparse
- Defines all subcommands: split, divide, stats, validate, chunk, split_samples, extract_sample, assign_samples, get_hto
- Uses custom formatter for help text
- Each subcommand calls a corresponding function from its module

### Command Modules

**bamurai/split.py**: Split reads to target length
- `calculate_split_len()`: Determines split locations to achieve fragments close to target length
- Reads shorter than target length are not split

**bamurai/divide.py**: Divide reads into fixed number of pieces
- `calculate_split_pieces()`: Splits reads into N equal pieces
- Respects minimum fragment length constraint

**bamurai/stats.py**: Calculate file statistics
- `file_read_stats()`: Returns total reads, average length, throughput, N50
  - For an empty file returns zeroed stats with `n50: 0` (for display), which
    differs from `calc_n50([])` returning `None`
- `calc_n50()`: Standard N50 calculation; sorts a copy, does not mutate the input list
- Supports TSV output format for computational analysis

**bamurai/validate.py**: File integrity validation
- Checks BAM headers, record structure, sequence/quality consistency
- For FASTQ: validates headers, separators, sequence characters

**bamurai/chunk.py**: Split files into size-based chunks
- Splits BAM/FASTQ files into chunks of at least specified size
- Uses human-readable size notation (1G, 100M, 1000K)

### Multi-Sample Processing

**bamurai/split_samples.py**: Split by donor ID
- `split_bam_by_donor()`: Creates temporary BAM files per donor
- `split_fastq_by_donor()`: Creates temporary FASTQ files per donor
- Uses temporary directory pattern to handle multiple input files
- Concatenates temp files into final donor-specific outputs
- Supports `--barcode-column` and `--donor-id-column` for flexible TSV formats

**bamurai/extract_sample.py**: Extract reads for specific donor
- Can process multiple BAM files at once
- Combines donor-specific reads into single output
- Supports `--barcode-column` and `--donor-id-column` for flexible TSV formats

**bamurai/assign_samples.py**: Assign donor IDs to barcodes
- Maps barcodes to donor IDs using TSV mapping
- Auto-detects barcode columns (same as split_samples and extract_sample)

**bamurai/get_hto.py**: Extract HTO (Hashtag Oligo) information
- Designed for 10x FASTQ files
- R1 contains cell barcode + UMI, R2 contains HTO sequence

### Utility Modules

**bamurai/utils.py**: General utilities
- `smart_open()`: Handles both regular and gzipped files transparently
- `create_progress_bar_for_file()`: Creates appropriate tqdm progress bars
- `count_reads_async_generic()`: Background thread for counting reads to update progress bar total
- `print_elapsed_time_pretty()`: Formats elapsed time based on duration

**bamurai/utils_samples.py**: Multi-sample utilities
- `parse_barcode_donor_mapping()`: Parses TSV with flexible column detection
  - Auto-detects 'barcode' or 'cell' columns for barcodes
  - Auto-detects 'donor_id' column for donor IDs
  - Raises error if both 'barcode' and 'cell' present without explicit specification
  - Accepts optional `barcode_column` and `donor_id_column` parameters for custom column names
- `get_read_barcode()`: Extracts barcode from read tags (CB, XC, or BC)
- `concatenate_bam_files()`: Merges multiple BAM files using pysam

**bamurai/logging_config.py**: Centralized logging configuration
- All commands use consistent logging format
- Call `configure_logging()` at start of command functions

## Key Patterns

### Read Processing Pattern
All read-processing commands follow this pattern:
1. Call `configure_logging()` and create logger
2. Set up output file (stdout or file, handle gzip)
3. Create progress bar with `create_progress_bar_for_file()`
4. Start async read counting with `count_reads_async_generic()`
5. Iterate through `parse_reads()` generator
6. Process each read and update progress
7. Close progress bar and output file
8. Log summary statistics

### File Type Detection
- Check file extensions: `.bam`, `.sam`, `.cram` for alignment files
- Use `is_fastq()` helper from utils for FASTQ detection
- `smart_open()` automatically handles `.gz` extensions

### Quality Score Conversion
BAM files store quality scores as integers; FASTQ uses ASCII characters. Use `qual_to_fastq_numpy()` for efficient conversion (adds 33 to each quality score).

### Multi-Sample Workflow
1. Parse barcode-to-donor TSV mapping
2. Create temp directory for intermediate files
3. Process each input file, writing to temp files per donor
4. Track all temp files in defaultdict
5. Concatenate temp files into final donor-specific outputs
6. Temp directory auto-cleaned via context manager

## Testing

### Static Test Inputs (important)

Tests read their inputs from **committed static files** in `tests/data/`. They do
*not* fabricate inputs per-test. Consequences:

- `tests/data/` is **tracked in git**, not ignored. The root `.gitignore` entry is
  anchored as `/data/` precisely so it does not also match `tests/data/`.
- `tests/generate_test_data.py` is the single source of truth for every file in
  `tests/data/`, built from the shared builders in `tests/conftest.py`.
- **After changing the generator or any builder, run `python tests/generate_test_data.py`.**
  `test_committed_data_matches_generator` regenerates into a temp dir and compares
  against `tests/data/`, failing if they diverge. It compares `.gz` files
  *decompressed* (gzip embeds an mtime, so raw bytes differ between runs);
  everything else, including BAM/BGZF (which uses mtime 0), is byte-compared.
- Every file in `tests/data/` must have an entry in `generate_test_data._MANIFEST`
  (enforced by `test_generate_writes_expected_files`).
- Tests **read** static inputs and **write** outputs only to pytest's `tmp_path`.
  Never write into `tests/data/`.

### tests/conftest.py

- `data_path(name)`: absolute path to a file in `tests/data/`
- Fixtures returning static paths: `fastq_file`, `fastq_gz_file`, `bam_file`,
  `barcoded_bam_file`, `barcode_donor_tsv`
- `make_args`: factory building an argparse-like namespace for calling command
  functions directly, bypassing the CLI
- Shared constants: `DONOR1_BARCODE`, `DONOR2_BARCODE`, `HTO_BARCODE`, `HTO_UMI`,
  `HTO_HASHTAG`, `HTO_LEFT_BUFFER` (kept in sync with the generated files)
- Deterministic builders (`make_sequence`, `write_fastq`, `write_bam`,
  `make_segment`, `default_bam_specs`, ...) used by `generate_test_data.py` and by
  unit tests that construct in-memory objects

### Conventions

- Functions taking **objects** rather than files are tested with constructed
  objects, not static files (e.g. `get_read_barcode(read)`, `calc_n50(list)`,
  `split_read(read, at)`). Only *file* inputs are static.
- The BAM fixture deliberately contains one secondary (flag 256) and one
  supplementary (flag 2048) record so the primary-only contract is falsifiable
  rather than assumed.
- Read-transforming commands are covered by an invariant test asserting total
  bases in == total bases out.
- Coverage deliberately excludes plumbing that carries no logic worth asserting:
  progress bars, background counting threads, `except Exception` handlers, and
  `# pragma: no cover` blocks (see `[tool.coverage.report] exclude_also` in
  `pyproject.toml`). Prefer excluding such code over writing brittle
  UI/threading tests for it.

## Version Management

Version is stored in `bamurai/VERSION` file and read via `bamurai/version.py:get_version()`. Update VERSION file for releases.
