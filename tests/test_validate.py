"""Tests for bamurai.validate: FASTQ and BAM integrity checks.

Inputs are static objects under tests/data/ (see generate_test_data.py).
"""

from bamurai.validate import validate_fastq, validate_bam, validate_file
from conftest import data_path


class TestValidateFastq:
    def test_valid_file(self, fastq_file):
        assert validate_fastq(fastq_file) is True

    def test_valid_gzipped(self, fastq_gz_file):
        assert validate_fastq(fastq_gz_file) is True

    def test_bad_header(self):
        # missing leading '@'
        assert validate_fastq(data_path("bad_header.fastq")) is False

    def test_bad_separator(self):
        # '-' instead of '+'
        assert validate_fastq(data_path("bad_separator.fastq")) is False

    def test_length_mismatch(self):
        # qual shorter than seq
        assert validate_fastq(data_path("length_mismatch.fastq")) is False

    def test_invalid_characters(self):
        # 'Z' not IUPAC
        assert validate_fastq(data_path("invalid_chars.fastq")) is False

    def test_iupac_ambiguity_codes_allowed(self):
        assert validate_fastq(data_path("iupac_ok.fastq")) is True

    def test_truncated_record_missing_lines(self):
        # Header present but sequence/separator/quality lines missing.
        assert validate_fastq(data_path("truncated.fastq")) is False

    def test_empty_file_is_valid_with_zero_records(self):
        assert validate_fastq(data_path("empty.fastq")) is True


class TestValidateBam:
    def test_valid_file(self, bam_file):
        assert validate_bam(bam_file) is True

    def test_missing_file_returns_false(self, tmp_path):
        assert validate_bam(str(tmp_path / "nope.bam")) is False

    def test_missing_sequence_returns_false(self, capsys):
        # A record with a name but no sequence must fail validation.
        assert validate_bam(data_path("noseq.bam")) is False
        assert "Missing query sequence" in capsys.readouterr().out

    def test_missing_qualities_returns_false(self, capsys):
        # A record with a sequence but absent qualities must fail validation.
        assert validate_bam(data_path("noqual.bam")) is False
        assert "Missing query qualities" in capsys.readouterr().out


class TestValidateFileDispatch:
    def test_dispatch_fastq(self, fastq_file, make_args, capsys):
        validate_file(make_args(reads=fastq_file))
        assert "valid FASTQ file" in capsys.readouterr().out

    def test_dispatch_bam(self, bam_file, make_args, capsys):
        validate_file(make_args(reads=bam_file))
        assert "valid BAM file" in capsys.readouterr().out

    def test_unsupported_extension(self, make_args, capsys):
        validate_file(make_args(reads=data_path("unsupported.txt")))
        assert "BAM or FASTQ" in capsys.readouterr().out
