"""Tests for bamurai.split: target-length splitting."""

from bamurai.core import Read, parse_reads
from bamurai.split import calculate_split_len, split_reads
from conftest import make_sequence, make_qualities, qual_ints_to_ascii


def _read(length):
    return Read("r", make_sequence(length), qual_ints_to_ascii(make_qualities(length)))


class TestCalculateSplitLen:
    def test_shorter_than_target_not_split(self):
        assert calculate_split_len(_read(50), target_len=100) == []

    def test_exactly_target_not_split(self):
        # round(100/100) == 1 -> range(1, 1) is empty
        assert calculate_split_len(_read(100), target_len=100) == []

    def test_two_fragments(self):
        # 250 / 100 -> round(2.5) == 2 splits, split_size = 125
        assert calculate_split_len(_read(250), target_len=100) == [125]

    def test_three_fragments(self):
        # 300 / 100 -> 3 pieces, split_size = 100
        assert calculate_split_len(_read(300), target_len=100) == [100, 200]

    def test_fragments_near_target_length(self):
        read = _read(1000)
        locs = calculate_split_len(read, target_len=200)
        fragments = [locs[0]] + [locs[i] - locs[i - 1] for i in range(1, len(locs))]
        fragments.append(len(read) - locs[-1])
        # Every fragment should be within a reasonable band around the target.
        assert all(150 <= f <= 250 for f in fragments)


class TestSplitReadsIntegration:
    def test_output_file_written(self, fastq_file, tmp_path, make_args):
        out = tmp_path / "split.fastq"
        args = make_args(reads=fastq_file, len_target=100, output=str(out))
        split_reads(args)
        assert out.exists()
        result = list(parse_reads(str(out)))
        # Reads of length [50,120,250,80,300] at target 100 ->
        # [1, 1, 2, 1, 3] = 8 fragments.
        assert len(result) == 8

    def test_total_bases_preserved(self, fastq_file, tmp_path, make_args):
        out = tmp_path / "split.fastq"
        in_bases = sum(len(r) for r in parse_reads(fastq_file))
        args = make_args(reads=fastq_file, len_target=100, output=str(out))
        split_reads(args)
        out_bases = sum(len(r) for r in parse_reads(str(out)))
        assert in_bases == out_bases

    def test_works_on_bam_input(self, bam_file, tmp_path, make_args):
        out = tmp_path / "split.fastq"
        args = make_args(reads=bam_file, len_target=100, output=str(out))
        split_reads(args)
        assert len(list(parse_reads(str(out)))) == 8
