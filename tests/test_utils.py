"""Tests for bamurai.utils: general helpers."""

import gzip

import pytest

from bamurai.utils import (
    is_fastq,
    smart_open,
    calculate_percentage,
    print_elapsed_time_pretty,
)


class TestIsFastq:
    @pytest.mark.parametrize("name", [
        "x.fastq", "x.fq", "x.fastq.gz", "x.fq.gz",
        "X.FASTQ", "path/to/reads.FQ.GZ",
    ])
    def test_true(self, name):
        assert is_fastq(name)

    @pytest.mark.parametrize("name", ["x.bam", "x.sam", "x.cram", "x.txt"])
    def test_false(self, name):
        assert not is_fastq(name)


class TestSmartOpen:
    def test_plain_roundtrip(self, tmp_path):
        path = tmp_path / "f.txt"
        with smart_open(str(path), "wt", encoding="utf-8") as f:
            f.write("hello")
        with smart_open(str(path), "rt", encoding="utf-8") as f:
            assert f.read() == "hello"

    def test_gzip_roundtrip(self, tmp_path):
        path = tmp_path / "f.txt.gz"
        with smart_open(str(path), "wt", encoding="utf-8") as f:
            f.write("hello gzip")
        # Confirm it is genuinely gzip-compressed on disk.
        with gzip.open(str(path), "rt") as f:
            assert f.read() == "hello gzip"

    def test_missing_file_raises_ioerror(self, tmp_path):
        with pytest.raises(IOError):
            smart_open(str(tmp_path / "nope.txt"), "rt")


class TestCalculatePercentage:
    def test_normal(self):
        assert calculate_percentage(1, 4) == 25.0

    def test_zero_total_is_safe(self):
        assert calculate_percentage(5, 0) == 0

    def test_full(self):
        assert calculate_percentage(10, 10) == 100.0


class TestPrintElapsedTime:
    def test_seconds_branch(self, caplog):
        import time
        with caplog.at_level("INFO"):
            print_elapsed_time_pretty(time.time() - 5)
        assert "Time taken" in caplog.text

    def test_minutes_branch(self, caplog):
        import time
        # ~400s elapsed -> the "Xm Ys" branch (300 <= elapsed < 3600).
        with caplog.at_level("INFO"):
            print_elapsed_time_pretty(time.time() - 400)
        assert "Time taken" in caplog.text
        assert "6m" in caplog.text

    def test_hours_branch(self, caplog):
        import time
        # ~4000s elapsed -> the "Xh Ym Zs" branch (elapsed >= 3600).
        with caplog.at_level("INFO"):
            print_elapsed_time_pretty(time.time() - 4000)
        assert "Time elapsed" in caplog.text
        assert "1h" in caplog.text
