"""Tests for bamurai.stats: N50 and file statistics."""

from bamurai.stats import calc_n50, file_read_stats, file_stats
from conftest import data_path


class TestCalcN50:
    def test_empty(self):
        assert calc_n50([]) is None

    def test_single_value(self):
        assert calc_n50([100]) == 100

    def test_known_value(self):
        # Lengths 2,3,4,5,6 -> total 20, half 10. Sorted desc: 6,5,4,3,2;
        # cumulative 6,11 -> reaches half at 5.
        assert calc_n50([2, 3, 4, 5, 6]) == 5

    def test_uniform_lengths(self):
        assert calc_n50([100, 100, 100, 100]) == 100

    def test_does_not_mutate_input(self):
        # calc_n50 must not reorder the caller's list.
        lengths = [2, 5, 3, 6, 4]
        calc_n50(lengths)
        assert lengths == [2, 5, 3, 6, 4]


class TestFileReadStats:
    def test_fastq_stats(self, fastq_file):
        stats = file_read_stats(fastq_file)
        # Lengths [50,120,250,80,300]: total 800, avg 160.
        assert stats["total_reads"] == 5
        assert stats["throughput"] == 800
        assert stats["avg_read_len"] == 160
        # Sorted desc 300,250,...; cumulative 300, 550 crosses half (400) at 250.
        assert stats["n50"] == 250

    def test_bam_stats_primary_only(self, bam_file):
        stats = file_read_stats(bam_file)
        # Secondary/supplementary must not inflate the count.
        assert stats["total_reads"] == 5
        assert stats["throughput"] == 800

    def test_fastq_and_bam_agree(self, fastq_file, bam_file):
        assert file_read_stats(fastq_file) == file_read_stats(bam_file)

    def test_empty_file_stats(self):
        # An empty FASTQ yields zeroed stats (n50 coerced to 0 for display,
        # distinct from calc_n50([]) which returns None).
        stats = file_read_stats(data_path("empty.fastq"))
        assert stats == {
            "total_reads": 0,
            "avg_read_len": 0,
            "throughput": 0,
            "n50": 0,
        }


class TestFileStatsOutput:
    def test_human_readable(self, fastq_file, make_args, capsys):
        file_stats(make_args(reads=fastq_file, tsv=False))
        out = capsys.readouterr().out
        assert "Total reads: 5" in out
        assert "N50: 250" in out

    def test_human_readable_throughput_gb(self, fastq_file, make_args, capsys):
        # 800 bases -> 800 / 1e9 rounded to 2 dp -> 0.0 Gb.
        file_stats(make_args(reads=fastq_file, tsv=False))
        out = capsys.readouterr().out
        assert "Throughput (Gb): 0.0" in out

    def test_tsv_output(self, fastq_file, make_args, capsys):
        file_stats(make_args(reads=fastq_file, tsv=True))
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert lines[0] == "file_name\ttotal_reads\tavg_read_len\tthroughput\tn50"
        fields = lines[1].split("\t")
        assert fields[1] == "5"  # total_reads
        assert fields[3] == "800"  # throughput
