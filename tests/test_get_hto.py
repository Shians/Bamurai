"""Tests for bamurai.get_hto: HTO extraction from 10x FASTQ pairs.

Inputs are static objects under tests/data/ (see generate_test_data.py).
"""

import csv

from bamurai.get_hto import get_hto
from conftest import data_path, HTO_BARCODE, HTO_UMI, HTO_HASHTAG


def _hto_args(make_args, r1, r2, out):
    return make_args(
        r1=data_path(r1), r2=data_path(r2), bc_len=16, umi_len=12,
        output=str(out), hashtag_len=15, hashtag_left_buffer=10,
    )


class TestGetHto:
    def test_extracts_barcode_umi_hto(self, tmp_path, make_args):
        out = tmp_path / "hto.tsv"
        get_hto(_hto_args(make_args, "hto_R1.fastq", "hto_R2.fastq", out))

        with open(str(out)) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))

        assert len(rows) == 1
        row = rows[0]
        assert row["read_name"] == "read1"
        assert row["cell_barcode"] == HTO_BARCODE
        assert row["umi"] == HTO_UMI
        assert row["hto"] == HTO_HASHTAG

    def test_quality_values_computed(self, tmp_path, make_args):
        # All qualities are 'I' (Phred 40), so every segment average is 40.0.
        out = tmp_path / "hto.tsv"
        get_hto(_hto_args(make_args, "hto_R1.fastq", "hto_R2.fastq", out))

        with open(str(out)) as f:
            row = next(csv.DictReader(f, delimiter="\t"))
        assert float(row["bc_qual"]) == 40.0
        assert float(row["umi_qual"]) == 40.0
        assert float(row["hto_qual"]) == 40.0

    def test_empty_segment_quality_is_zero(self, tmp_path, make_args):
        # A left buffer past the end of R2 (only ~30bp) yields an empty HTO
        # segment; avg_qual("") must return 0 rather than divide by zero.
        out = tmp_path / "hto.tsv"
        args = make_args(
            r1=data_path("hto_R1.fastq"), r2=data_path("hto_R2.fastq"),
            bc_len=16, umi_len=12, output=str(out),
            hashtag_len=15, hashtag_left_buffer=40,
        )
        get_hto(args)

        with open(str(out)) as f:
            row = next(csv.DictReader(f, delimiter="\t"))
        assert row["hto"] == ""
        assert float(row["hto_qual"]) == 0.0
        # The R1-derived segments are unaffected by the R2 buffer.
        assert float(row["bc_qual"]) == 40.0
        assert float(row["umi_qual"]) == 40.0

    def test_gzipped_input(self, tmp_path, make_args):
        # get_hto must transparently read gzipped R1/R2 files.
        out = tmp_path / "hto.tsv"
        get_hto(_hto_args(make_args, "hto_R1.fastq.gz", "hto_R2.fastq.gz", out))

        with open(str(out)) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        assert len(rows) == 1
        assert rows[0]["cell_barcode"] == HTO_BARCODE
        assert rows[0]["hto"] == HTO_HASHTAG

    def test_multiple_records(self, tmp_path, make_args):
        # Two read pairs must yield two output rows with distinct barcodes.
        out = tmp_path / "hto.tsv"
        get_hto(_hto_args(make_args, "hto_multi_R1.fastq", "hto_multi_R2.fastq", out))

        with open(str(out)) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        assert [r["read_name"] for r in rows] == ["readA", "readB"]
        assert [r["cell_barcode"] for r in rows] == ["A" * 16, "C" * 16]

    def test_quality_columns_present(self, tmp_path, make_args):
        out = tmp_path / "hto.tsv"
        get_hto(_hto_args(make_args, "hto_R1.fastq", "hto_R2.fastq", out))

        with open(str(out)) as f:
            header = f.readline().strip().split("\t")
        assert header == [
            "read_name", "cell_barcode", "umi", "hto",
            "bc_qual", "umi_qual", "hto_qual",
        ]
