"""Smoke tests for the static test-data generator.

These keep generate_test_data.py honest and, crucially, keep the committed
tests/data/ inputs in sync with it: if a builder signature drifts, an output
stops matching what the suite asserts, or a committed data file goes stale,
these fail.
"""

import gzip
import os

from bamurai.core import parse_reads
from bamurai.stats import file_read_stats
from conftest import DATA_DIR
import generate_test_data as gen


# The complete manifest of static inputs the suite consumes. Kept explicit
# (rather than derived) so it independently pins what the generator must emit.
EXPECTED_FILES = {
    # valid core inputs
    "reads.fastq", "reads.fastq.gz", "reads.bam", "reads.sam",
    "barcoded.bam", "barcoded.sam", "barcoded.fastq",
    # barcode->donor mapping TSVs (incl. column-detection variants)
    "mapping.tsv", "cell_donor.tsv", "both_columns.tsv",
    "missing_barcode.tsv", "missing_donor.tsv", "custom_columns.tsv",
    "numeric_barcodes.tsv",
    # malformed / edge-case FASTQ inputs (+ unsupported extension)
    "bad_header.fastq", "bad_separator.fastq", "length_mismatch.fastq",
    "invalid_chars.fastq", "truncated.fastq", "empty.fastq",
    "iupac_ok.fastq", "unsupported.txt",
    # HTO read pairs (plain, gzipped, multi-record)
    "hto_R1.fastq", "hto_R2.fastq", "hto_R1.fastq.gz", "hto_R2.fastq.gz",
    "hto_multi_R1.fastq", "hto_multi_R2.fastq",
    # crafted BAMs (tag variants, concat inputs, incomplete records)
    "rx.bam", "both_tags.bam", "concat_in_0.bam", "concat_in_1.bam",
    "noseq.bam", "noqual.bam",
    # manifest
    "MANIFEST.txt",
}


def _content(path):
    """Comparable content for a data file.

    ``.gz`` files embed an mtime in their gzip header, so raw bytes differ
    between runs; compare the decompressed payload instead. Everything else
    (including BAM/BGZF, which uses mtime 0) is byte-stable, so compare raw.
    """
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as f:
            return f.read()
    with open(path, "rb") as f:
        return f.read()


def test_generate_writes_expected_files(tmp_path):
    gen.generate(str(tmp_path))
    produced = set(os.listdir(str(tmp_path)))
    assert produced == EXPECTED_FILES
    # Every produced file except MANIFEST.txt must be documented in the manifest.
    assert set(gen._MANIFEST) | {"MANIFEST.txt"} == EXPECTED_FILES


def test_committed_data_matches_generator(tmp_path):
    """The committed tests/data/ files must match generate_test_data output.

    The suite reads its inputs from committed static files, so a change to
    generate_test_data.py (or a hand-edited data file) without regenerating
    would silently desync the two. Regenerate with
    ``python tests/generate_test_data.py`` to fix a failure here.
    """
    gen.generate(str(tmp_path))
    fresh = set(os.listdir(str(tmp_path)))
    committed = set(os.listdir(DATA_DIR))
    assert committed == fresh, (
        "tests/data/ file set is out of sync with the generator "
        "(run: python tests/generate_test_data.py) -- "
        f"missing={sorted(fresh - committed)}, extra={sorted(committed - fresh)}"
    )
    stale = [
        name for name in sorted(committed)
        if _content(os.path.join(DATA_DIR, name))
        != _content(os.path.join(str(tmp_path), name))
    ]
    assert not stale, (
        "committed tests/data/ files are stale "
        f"(run: python tests/generate_test_data.py): {stale}"
    )


def test_generated_bam_matches_fixture_semantics(tmp_path):
    gen.generate(str(tmp_path))
    # BAM must yield only the 5 primary reads (secondary/supplementary skipped).
    reads = list(parse_reads(str(tmp_path / "reads.bam")))
    assert [len(r) for r in reads] == [50, 120, 250, 80, 300]


def test_generated_fastq_and_bam_agree(tmp_path):
    gen.generate(str(tmp_path))
    fq_stats = file_read_stats(str(tmp_path / "reads.fastq"))
    bam_stats = file_read_stats(str(tmp_path / "reads.bam"))
    assert fq_stats == bam_stats
    assert fq_stats["total_reads"] == 5


def test_sam_twin_records_readable(tmp_path):
    import pysam
    gen.generate(str(tmp_path))
    with pysam.AlignmentFile(str(tmp_path / "reads.sam"), "r", check_sq=False) as sam:
        records = list(sam)
    # 5 primary records plus exactly one secondary and one supplementary.
    primary = [r for r in records if not (r.is_secondary or r.is_supplementary)]
    assert len(primary) == 5
    assert sum(r.is_secondary for r in records) == 1
    assert sum(r.is_supplementary for r in records) == 1
