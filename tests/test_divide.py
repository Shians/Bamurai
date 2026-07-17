"""Tests for bamurai.divide: fixed-fragment-count splitting."""

from bamurai.core import Read, parse_reads
from bamurai.divide import calculate_split_pieces, divide_reads
from conftest import make_sequence, make_qualities, qual_ints_to_ascii


def _read(length):
    return Read("r", make_sequence(length), qual_ints_to_ascii(make_qualities(length)))


class TestCalculateSplitPieces:
    def test_two_pieces(self):
        assert calculate_split_pieces(_read(100), num_pieces=2) == [50]

    def test_four_pieces(self):
        assert calculate_split_pieces(_read(100), num_pieces=4) == [25, 50, 75]

    def test_below_min_length_not_split(self):
        # 100 / 4 = 25 < min_length 50 -> no split
        assert calculate_split_pieces(_read(100), num_pieces=4, min_length=50) == []

    def test_at_min_length_boundary_splits(self):
        # 100 / 2 = 50, not < 50, so it splits
        assert calculate_split_pieces(_read(100), num_pieces=2, min_length=50) == [50]

    def test_uneven_division_uses_floor(self):
        # 101 // 3 = 33 -> [33, 66], last fragment absorbs remainder
        assert calculate_split_pieces(_read(101), num_pieces=3) == [33, 66]


class TestDivideReadsIntegration:
    def test_each_read_divided(self, fastq_file, tmp_path, make_args):
        out = tmp_path / "divide.fastq"
        args = make_args(reads=fastq_file, num_fragments=2, min_length=0,
                         output=str(out))
        divide_reads(args)
        # 5 reads, each into 2 pieces -> 10 output reads.
        assert len(list(parse_reads(str(out)))) == 10

    def test_min_length_prevents_short_splits(self, fastq_file, tmp_path, make_args):
        out = tmp_path / "divide.fastq"
        # num_pieces=2, min_length=100: per-fragment lengths are
        # [25,60,125,40,150] for reads [50,120,250,80,300]. Only the 250 and
        # 300 reads clear the threshold and get divided; the rest stay whole.
        args = make_args(reads=fastq_file, num_fragments=2, min_length=100,
                         output=str(out))
        divide_reads(args)
        # 3 unsplit (1 each) + 2 split (2 each) = 7.
        assert len(list(parse_reads(str(out)))) == 7

    def test_bases_preserved(self, fastq_file, tmp_path, make_args):
        out = tmp_path / "divide.fastq"
        in_bases = sum(len(r) for r in parse_reads(fastq_file))
        args = make_args(reads=fastq_file, num_fragments=3, min_length=0,
                         output=str(out))
        divide_reads(args)
        out_bases = sum(len(r) for r in parse_reads(str(out)))
        assert in_bases == out_bases
