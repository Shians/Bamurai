"""Tests for bamurai.core: the Read dataclass and parsing primitives."""

import pytest

from bamurai.core import (
    Read,
    qual_to_fastq_numpy,
    parse_reads,
    split_read,
)
from conftest import make_sequence, make_qualities, qual_ints_to_ascii


# ---------------------------------------------------------------------------
# Read dataclass
# ---------------------------------------------------------------------------

class TestRead:
    def test_valid_read(self):
        read = Read("r1", "ACGT", "IIII")
        assert read.read_id == "r1"
        assert read.sequence == "ACGT"
        assert read.quality == "IIII"

    def test_len(self):
        assert len(Read("r1", "ACGTACGT", "IIIIIIII")) == 8

    def test_empty_read_is_valid(self):
        read = Read("r1", "", "")
        assert len(read) == 0
        assert read.is_valid()

    def test_mismatched_lengths_raise(self):
        with pytest.raises(ValueError, match="equal length"):
            Read("bad", "ACGT", "II")

    def test_is_valid_true(self):
        assert Read("r1", "ACGT", "IIII").is_valid()

    def test_to_fastq_roundtrip_format(self):
        read = Read("r1", "ACGT", "IIII")
        assert read.to_fastq() == "@r1\nACGT\n+\nIIII"


# ---------------------------------------------------------------------------
# qual_to_fastq_numpy
# ---------------------------------------------------------------------------

class TestQualConversion:
    def test_known_values(self):
        # Phred 0 -> '!' (33), Phred 30 -> '?' (63), Phred 40 -> 'I' (73)
        assert qual_to_fastq_numpy([0, 30, 40]) == "!?I"

    def test_empty(self):
        assert qual_to_fastq_numpy([]) == ""

    def test_matches_manual_offset(self):
        quals = make_qualities(20, value=25)
        assert qual_to_fastq_numpy(quals) == qual_ints_to_ascii(quals)


# ---------------------------------------------------------------------------
# split_read
# ---------------------------------------------------------------------------

class TestSplitRead:
    def test_no_split_positions_appends_index(self):
        read = Read("r1", "ACGTACGT", "IIIIIIII")
        result = split_read(read, at=[])
        assert len(result) == 1
        assert result[0].read_id == "r1_0"
        assert result[0].sequence == "ACGTACGT"

    def test_no_split_does_not_mutate_input(self):
        # split_read must not rename the caller's read in place.
        read = Read("r1", "ACGTACGT", "IIIIIIII")
        split_read(read, at=[])
        assert read.read_id == "r1"

    def test_split_does_not_mutate_input(self):
        read = Read("r1", "AAAACCCC", "IIIIJJJJ")
        split_read(read, at=[4])
        assert read.read_id == "r1"
        assert read.sequence == "AAAACCCC"

    def test_single_split(self):
        read = Read("r1", "AAAACCCC", "IIIIJJJJ")
        result = split_read(read, at=[4])
        assert [r.sequence for r in result] == ["AAAA", "CCCC"]
        assert [r.quality for r in result] == ["IIII", "JJJJ"]
        assert [r.read_id for r in result] == ["r1_0", "r1_1"]

    def test_multiple_splits(self):
        read = Read("r1", "AABBCCDD", "12345678")
        result = split_read(read, at=[2, 4, 6])
        assert [r.sequence for r in result] == ["AA", "BB", "CC", "DD"]
        assert [r.read_id for r in result] == ["r1_0", "r1_1", "r1_2", "r1_3"]

    def test_fragments_cover_full_read(self):
        read = Read("r1", make_sequence(100), qual_ints_to_ascii(make_qualities(100)))
        result = split_read(read, at=[25, 50, 75])
        assert "".join(r.sequence for r in result) == read.sequence
        assert sum(len(r) for r in result) == 100


# ---------------------------------------------------------------------------
# parse_reads
# ---------------------------------------------------------------------------

class TestParseReadsFastq:
    def test_plain_fastq(self, fastq_file):
        reads = list(parse_reads(fastq_file))
        assert len(reads) == 5
        assert [len(r) for r in reads] == [50, 120, 250, 80, 300]
        assert reads[0].read_id == "read_0"

    def test_gzipped_fastq(self, fastq_gz_file):
        reads = list(parse_reads(fastq_gz_file))
        assert len(reads) == 5
        assert [len(r) for r in reads] == [50, 120, 250, 80, 300]

    def test_fastq_strips_at_symbol(self, fastq_file):
        reads = list(parse_reads(fastq_file))
        assert not reads[0].read_id.startswith("@")

    def test_fastq_all_valid(self, fastq_file):
        assert all(r.is_valid() for r in parse_reads(fastq_file))


class TestParseReadsBam:
    def test_primary_only(self, bam_file):
        """Secondary and supplementary alignments must be filtered out."""
        reads = list(parse_reads(bam_file))
        assert len(reads) == 5
        assert [len(r) for r in reads] == [50, 120, 250, 80, 300]

    def test_read_ids(self, bam_file):
        reads = list(parse_reads(bam_file))
        assert [r.read_id for r in reads] == [
            "read_0", "read_1", "read_2", "read_3", "read_4",
        ]

    def test_qualities_offset_by_33(self, bam_file):
        # First read built with quality value 30 -> '?'
        reads = list(parse_reads(bam_file))
        assert reads[0].quality[0] == "?"
