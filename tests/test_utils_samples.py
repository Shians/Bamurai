"""Tests for bamurai.utils_samples: barcode/donor mapping helpers.

Inputs are static objects under tests/data/ (see generate_test_data.py).
"""

import os

import pysam
import pytest

from bamurai.utils_samples import (
    parse_barcode_donor_mapping,
    get_barcodes_for_donor,
    ensure_directory_exists,
    get_read_barcode,
    concatenate_bam_files,
)
from conftest import (
    data_path,
    default_bam_header,
    make_segment,
    DONOR1_BARCODE,
    DONOR2_BARCODE,
)


class TestParseBarcodeDonorMapping:
    def test_barcode_column(self):
        assert parse_barcode_donor_mapping(data_path("mapping.tsv")) == {
            DONOR1_BARCODE: "donor1", DONOR2_BARCODE: "donor2"}

    def test_cell_column_autodetected(self):
        assert parse_barcode_donor_mapping(data_path("cell_donor.tsv")) == {"AAA": "d1"}

    def test_both_barcode_and_cell_raises(self):
        with pytest.raises(ValueError, match="specify --barcode-column"):
            parse_barcode_donor_mapping(data_path("both_columns.tsv"))

    def test_missing_barcode_column_raises(self):
        with pytest.raises(ValueError, match="No 'barcode' or 'cell'"):
            parse_barcode_donor_mapping(data_path("missing_barcode.tsv"))

    def test_missing_donor_column_raises(self):
        with pytest.raises(ValueError, match="No 'donor_id'"):
            parse_barcode_donor_mapping(data_path("missing_donor.tsv"))

    def test_custom_columns(self):
        result = parse_barcode_donor_mapping(
            data_path("custom_columns.tsv"), barcode_column="bc",
            donor_id_column="sample")
        assert result == {"AAA": "d1"}

    def test_custom_column_not_found_raises(self):
        with pytest.raises(ValueError, match="not found in TSV"):
            parse_barcode_donor_mapping(
                data_path("mapping.tsv"), barcode_column="missing")


class TestGetBarcodesForDonor:
    def test_filters_by_donor(self):
        mapping = {"AAA": "d1", "CCC": "d2", "GGG": "d1"}
        assert get_barcodes_for_donor(mapping, "d1") == {"AAA", "GGG"}

    def test_unknown_donor_empty(self):
        assert get_barcodes_for_donor({"AAA": "d1"}, "d9") == set()


class TestEnsureDirectoryExists:
    def test_creates_nested_dirs(self, tmp_path):
        target = tmp_path / "a" / "b" / "file.bam"
        ensure_directory_exists(str(target))
        assert os.path.isdir(str(tmp_path / "a" / "b"))

    def test_no_directory_component_is_noop(self):
        ensure_directory_exists("file.bam")  # should not raise


class TestGetReadBarcode:
    def _read(self, tag=None):
        header = default_bam_header()
        tags = [(tag, "AAACCC", "Z")] if tag else None
        return make_segment(header, "r", "ACGT", [30, 30, 30, 30], tags=tags)

    @pytest.mark.parametrize("tag", ["CB", "XC", "BC"])
    def test_reads_known_tags(self, tag):
        assert get_read_barcode(self._read(tag)) == "AAACCC"

    def test_no_barcode_returns_none(self):
        assert get_read_barcode(self._read(None)) is None

    def test_tag_precedence_cb_over_bc(self):
        # CB is checked before BC, so CB wins when both are present.
        header = default_bam_header()
        read = make_segment(
            header, "r", "ACGT", [30, 30, 30, 30],
            tags=[("CB", "FROMCB", "Z"), ("BC", "FROMBC", "Z")],
        )
        assert get_read_barcode(read) == "FROMCB"


class TestConcatenateBamFiles:
    def test_merges_records(self, tmp_path, capsys):
        paths = [data_path("concat_in_0.bam"), data_path("concat_in_1.bam")]
        out = tmp_path / "merged.bam"
        concatenate_bam_files(paths, str(out))
        with pysam.AlignmentFile(str(out), "rb") as bam:
            names = [r.query_name for r in bam]
        assert names == ["r0", "r1"]

    def test_empty_list_is_noop(self, tmp_path):
        out = tmp_path / "merged.bam"
        concatenate_bam_files([], str(out))
        assert not out.exists()
