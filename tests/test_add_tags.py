"""Tests for the add_tags command.

Uses the shared `bam_file` fixture (5 primary reads read_0..read_4, plus a
secondary alignment named read_0 and a supplementary alignment named read_2)
to exercise tagging, unmatched/blank-cell handling, and the primary-only rule.
"""

import pysam
import pytest

from bamurai.add_tags import add_tags, add_tags_file


def _write_tsv(path, headers, rows):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(headers) + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")


def _reads_by_name(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return list(bam)


class TestAddTags:
    def test_tags_matching_reads(self, bam_file, tmp_path, make_args):
        tsv = tmp_path / "tags.tsv"
        out = tmp_path / "out.bam"
        _write_tsv(tsv, ["read_id", "XX", "YY"], [["read_0", "hello", "1"], ["read_1", "world", "2"]])

        add_tags(make_args(bam=bam_file, tsv=str(tsv), output=str(out)))

        reads = {r.query_name: r for r in _reads_by_name(str(out)) if not (r.is_secondary or r.is_supplementary)}
        assert reads["read_0"].get_tag("XX") == "hello"
        assert reads["read_0"].get_tag("YY") == 1
        assert reads["read_1"].get_tag("XX") == "world"
        assert reads["read_1"].get_tag("YY") == 2

    def test_unmatched_read_passes_through_untagged(self, bam_file, tmp_path, make_args):
        tsv = tmp_path / "tags.tsv"
        out = tmp_path / "out.bam"
        _write_tsv(tsv, ["read_id", "XX"], [["read_0", "hello"]])

        add_tags(make_args(bam=bam_file, tsv=str(tsv), output=str(out)))

        reads = {r.query_name: r for r in _reads_by_name(str(out)) if not (r.is_secondary or r.is_supplementary)}
        assert reads["read_0"].get_tag("XX") == "hello"
        assert not reads["read_1"].has_tag("XX")

    def test_blank_cell_skips_only_that_tag(self, bam_file, tmp_path, make_args):
        tsv = tmp_path / "tags.tsv"
        out = tmp_path / "out.bam"
        _write_tsv(tsv, ["read_id", "XX", "YY"], [["read_0", "", "42"]])

        add_tags(make_args(bam=bam_file, tsv=str(tsv), output=str(out)))

        reads = {r.query_name: r for r in _reads_by_name(str(out)) if not (r.is_secondary or r.is_supplementary)}
        assert not reads["read_0"].has_tag("XX")
        assert reads["read_0"].get_tag("YY") == 42

    def test_secondary_alignment_not_tagged(self, bam_file, tmp_path, make_args):
        tsv = tmp_path / "tags.tsv"
        out = tmp_path / "out.bam"
        _write_tsv(tsv, ["read_id", "XX"], [["read_0", "hello"]])

        add_tags(make_args(bam=bam_file, tsv=str(tsv), output=str(out)))

        all_reads = _reads_by_name(str(out))
        primary = next(r for r in all_reads if r.query_name == "read_0" and not r.is_secondary)
        secondary = next(r for r in all_reads if r.query_name == "read_0" and r.is_secondary)
        assert primary.get_tag("XX") == "hello"
        assert not secondary.has_tag("XX")

    def test_invalid_tag_column_length_raises(self, bam_file, tmp_path, make_args):
        tsv = tmp_path / "tags.tsv"
        out = tmp_path / "out.bam"
        _write_tsv(tsv, ["read_id", "XXX"], [["read_0", "hello"]])

        with pytest.raises(ValueError):
            add_tags(make_args(bam=bam_file, tsv=str(tsv), output=str(out)))

    def test_missing_read_id_column_raises(self, bam_file, tmp_path, make_args):
        tsv = tmp_path / "tags.tsv"
        out = tmp_path / "out.bam"
        _write_tsv(tsv, ["id", "XX"], [["read_0", "hello"]])

        with pytest.raises(ValueError):
            add_tags(make_args(bam=bam_file, tsv=str(tsv), output=str(out)))

    def test_no_tag_columns_raises(self, bam_file, tmp_path, make_args):
        tsv = tmp_path / "tags.tsv"
        out = tmp_path / "out.bam"
        _write_tsv(tsv, ["read_id"], [["read_0"]])

        with pytest.raises(ValueError):
            add_tags(make_args(bam=bam_file, tsv=str(tsv), output=str(out)))

    def test_add_tags_file_helper_used_directly(self, bam_file, tmp_path):
        tsv = tmp_path / "tags.tsv"
        out = tmp_path / "out.bam"
        _write_tsv(tsv, ["read_id", "XX"], [["read_0", "hello"]])

        add_tags_file(bam_file, str(out), str(tsv))

        reads = {r.query_name: r for r in _reads_by_name(str(out)) if not (r.is_secondary or r.is_supplementary)}
        assert reads["read_0"].get_tag("XX") == "hello"
