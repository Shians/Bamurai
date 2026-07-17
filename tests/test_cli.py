"""End-to-end tests that drive the argparse entry point in bamurai.cli.

These exercise the CLI wiring itself (argument parsing and dispatch to the
command functions) rather than the command internals, which are covered
elsewhere. Each test patches ``sys.argv`` and calls ``main()``.
"""

import sys

import pytest

from bamurai.cli import main
from conftest import data_path


def _run(monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", ["bamurai"] + argv)
    main()


class TestCliDispatch:
    def test_stats(self, monkeypatch, capsys, fastq_file):
        _run(monkeypatch, ["stats", fastq_file])
        out = capsys.readouterr().out
        assert "Total reads: 5" in out

    def test_stats_tsv(self, monkeypatch, capsys, bam_file):
        _run(monkeypatch, ["stats", bam_file, "--tsv"])
        out = capsys.readouterr().out
        assert out.startswith("file_name\t")

    def test_split_to_output(self, monkeypatch, tmp_path, fastq_file):
        out = tmp_path / "split.fastq"
        _run(monkeypatch, ["split", fastq_file, "-l", "100", "-o", str(out)])
        assert out.exists()

    def test_divide_to_output(self, monkeypatch, tmp_path, fastq_file):
        out = tmp_path / "divide.fastq"
        _run(monkeypatch, ["divide", fastq_file, "-n", "2", "-o", str(out)])
        assert out.exists()

    def test_validate(self, monkeypatch, capsys, fastq_file):
        _run(monkeypatch, ["validate", fastq_file])
        assert "valid FASTQ file" in capsys.readouterr().out

    def test_chunk(self, monkeypatch, tmp_path, fastq_file):
        prefix = str(tmp_path / "chunk")
        _run(monkeypatch, ["chunk", fastq_file, "-s", "1K", "-p", prefix])
        assert (tmp_path / "chunk_1.fastq").exists()


class TestCliSampleCommands:
    """Smoke tests for the multi-sample and HTO subcommands via argparse.

    These verify the argparse wiring (flag names, nargs, dest mapping) for the
    subcommands that were previously only invoked as direct function calls.
    """

    def test_split_samples(self, monkeypatch, tmp_path, barcoded_bam_file,
                           barcode_donor_tsv):
        out_dir = tmp_path / "split_out"
        _run(monkeypatch, [
            "split_samples", "--input", barcoded_bam_file,
            "--tsv", barcode_donor_tsv, "--output-dir", str(out_dir),
        ])
        assert (out_dir / "donor1.bam").exists()

    def test_extract_sample(self, monkeypatch, tmp_path, barcoded_bam_file,
                            barcode_donor_tsv):
        out = tmp_path / "donor1.bam"
        _run(monkeypatch, [
            "extract_sample", "--bam", barcoded_bam_file,
            "--tsv", barcode_donor_tsv, "--donor-id", "donor1",
            "--output", str(out),
        ])
        assert out.exists()

    def test_assign_samples(self, monkeypatch, tmp_path, barcoded_bam_file,
                            barcode_donor_tsv):
        out = tmp_path / "assigned.bam"
        _run(monkeypatch, [
            "assign_samples", "--bam", barcoded_bam_file,
            "--tsv", barcode_donor_tsv, "--output", str(out),
        ])
        assert out.exists()

    def test_get_hto(self, monkeypatch, tmp_path):
        out = tmp_path / "hto.tsv"
        _run(monkeypatch, [
            "get_hto",
            "--r1", data_path("hto_R1.fastq"),
            "--r2", data_path("hto_R2.fastq"),
            "--bc-len", "16", "--umi-len", "12", "--output", str(out),
        ])
        assert out.exists()


class TestCliMeta:
    def test_version(self, monkeypatch, capsys):
        with pytest.raises(SystemExit) as exc:
            _run(monkeypatch, ["--version"])
        assert exc.value.code == 0
        assert capsys.readouterr().out.strip()  # prints a version string

    def test_no_command_prints_help(self, monkeypatch, capsys):
        _run(monkeypatch, [])
        out = capsys.readouterr().out
        assert "Bamurai" in out
        assert "usage" in out.lower()
