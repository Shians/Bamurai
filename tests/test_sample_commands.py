"""End-to-end tests for the multi-sample commands.

Covers split_samples, extract_sample and assign_samples using the barcoded
BAM fixture and its matching barcode->donor TSV. Recall that the fixture has:
  * 2 reads with barcode AAA... -> donor1
  * 2 reads with barcode CCC... -> donor2
  * 1 read with an unknown barcode -> unmapped
"""

import pysam
import pytest

from bamurai.split_samples import split_samples
from bamurai.extract_sample import extract_sample, extract_reads_from_bam
from bamurai.assign_samples import assign_samples
from bamurai.utils_samples import parse_barcode_donor_mapping, get_barcodes_for_donor
from conftest import data_path


def _count(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return sum(1 for _ in bam)


def _fastq_read_names(path):
    """Return the read names (first whitespace/tab token) in a FASTQ file."""
    names = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                names.append(line[1:].split()[0].split("\t")[0])
    return names


class TestSplitSamples:
    def test_bam_split_by_donor(self, barcoded_bam_file, barcode_donor_tsv,
                                tmp_path, make_args):
        out_dir = tmp_path / "split_out"
        args = make_args(
            input=[barcoded_bam_file],
            tsv=barcode_donor_tsv,
            output_dir=str(out_dir),
            barcode_column=None,
            donor_id_column=None,
        )
        split_samples(args)

        assert _count(str(out_dir / "donor1.bam")) == 2
        assert _count(str(out_dir / "donor2.bam")) == 2
        assert _count(str(out_dir / "unmapped.bam")) == 1

    def test_fastq_split_by_donor(self, tmp_path, barcode_donor_tsv, make_args):
        # FASTQ split reads the barcode from a BC:Z: field in the read name.
        out_dir = tmp_path / "split_out"
        args = make_args(
            input=[data_path("barcoded.fastq")],
            tsv=barcode_donor_tsv,
            output_dir=str(out_dir),
            barcode_column=None,
            donor_id_column=None,
        )
        split_samples(args)

        # Every donor (and the unmapped bucket) must receive exactly its read.
        assert _fastq_read_names(str(out_dir / "donor1.fastq")) == ["r0"]
        assert _fastq_read_names(str(out_dir / "donor2.fastq")) == ["r1"]
        assert _fastq_read_names(str(out_dir / "unmapped.fastq")) == ["r2"]

    def test_bam_split_multiple_inputs(self, barcoded_bam_file,
                                       barcode_donor_tsv, tmp_path, make_args):
        # Two input BAMs exercise the per-donor concatenation path (each donor
        # gets >1 temp file). The same barcoded BAM twice doubles every bucket.
        out_dir = tmp_path / "split_out"
        args = make_args(
            input=[barcoded_bam_file, barcoded_bam_file],
            tsv=barcode_donor_tsv,
            output_dir=str(out_dir),
            barcode_column=None,
            donor_id_column=None,
        )
        split_samples(args)

        assert _count(str(out_dir / "donor1.bam")) == 4
        assert _count(str(out_dir / "donor2.bam")) == 4
        assert _count(str(out_dir / "unmapped.bam")) == 2


class TestExtractSample:
    def test_extract_single_donor(self, barcoded_bam_file, barcode_donor_tsv,
                                   tmp_path, make_args):
        out = tmp_path / "donor1.bam"
        args = make_args(
            bam=[barcoded_bam_file],
            tsv=barcode_donor_tsv,
            donor_id="donor1",
            output=str(out),
            barcode_column=None,
            donor_id_column=None,
        )
        extract_sample(args)
        assert _count(str(out)) == 2

    def test_extract_reads_from_bam_helper(self, barcoded_bam_file,
                                           barcode_donor_tsv, tmp_path):
        mapping = parse_barcode_donor_mapping(barcode_donor_tsv)
        barcodes = get_barcodes_for_donor(mapping, "donor2")
        out = tmp_path / "d2.bam"
        count, path = extract_reads_from_bam(
            barcoded_bam_file, barcodes, output_file=str(out))
        assert count == 2
        assert _count(path) == 2

    def test_multiple_bam_inputs(self, barcoded_bam_file, barcode_donor_tsv,
                                 tmp_path, make_args):
        out = tmp_path / "donor1.bam"
        args = make_args(
            bam=[barcoded_bam_file, barcoded_bam_file],
            tsv=barcode_donor_tsv,
            donor_id="donor1",
            output=str(out),
            barcode_column=None,
            donor_id_column=None,
        )
        extract_sample(args)
        # Same donor across two identical inputs -> 2 + 2 reads.
        assert _count(str(out)) == 4

    def test_no_barcodes_for_donor(self, barcoded_bam_file, barcode_donor_tsv,
                                   tmp_path, make_args, capsys):
        # A donor absent from the mapping yields no barcodes: extract_sample
        # must report it and return without writing an output file.
        out = tmp_path / "missing.bam"
        args = make_args(
            bam=[barcoded_bam_file],
            tsv=barcode_donor_tsv,
            donor_id="donorX",  # not present in the mapping
            output=str(out),
            barcode_column=None,
            donor_id_column=None,
        )
        extract_sample(args)
        assert "No barcodes found for donor_id 'donorX'" in capsys.readouterr().out
        assert not out.exists()


class TestAssignSamples:
    def test_rg_tags_assigned(self, barcoded_bam_file, barcode_donor_tsv,
                              tmp_path, make_args):
        out = tmp_path / "assigned.bam"
        args = make_args(
            bam=barcoded_bam_file,
            tsv=barcode_donor_tsv,
            output=str(out),
            barcode_column=None,
            donor_id_column=None,
        )
        assign_samples(args)

        rg_by_name = {}
        with pysam.AlignmentFile(str(out), "rb") as bam:
            # RG lines for both donors should be present in the header.
            rg_ids = {rg["ID"] for rg in bam.header.to_dict().get("RG", [])}
            assert {"donor1", "donor2"}.issubset(rg_ids)
            for read in bam:
                rg_by_name[read.query_name] = (
                    read.get_tag("RG") if read.has_tag("RG") else None)

        assert rg_by_name["bc_read_0"] == "donor1"
        assert rg_by_name["bc_read_2"] == "donor2"
        assert rg_by_name["bc_read_4"] is None  # unknown barcode, unassigned

    def test_rx_tag_does_not_assign_a_donor(self, tmp_path, barcode_donor_tsv,
                                            make_args):
        # RX holds the UMI, not a cell barcode, so it cannot name a donor even
        # when its value happens to match one. rx.bam holds a single read
        # tagged RX=<donor1 barcode>; it must come back unassigned.
        out = tmp_path / "assigned.bam"
        args = make_args(
            bam=data_path("rx.bam"), tsv=barcode_donor_tsv, output=str(out),
            barcode_column=None, donor_id_column=None,
        )
        assign_samples(args)

        with pysam.AlignmentFile(str(out), "rb") as bam:
            read = next(iter(bam))
            assert not read.has_tag("RG")

    def test_rx_ignored_when_cb_present(self, tmp_path, barcode_donor_tsv,
                                        make_args):
        # both_tags.bam holds a read with CB=donor1 barcode and RX=donor2
        # barcode. CB names the cell, so donor1 wins and RX is never consulted.
        out = tmp_path / "assigned.bam"
        args = make_args(
            bam=data_path("both_tags.bam"), tsv=barcode_donor_tsv,
            output=str(out), barcode_column=None, donor_id_column=None,
        )
        assign_samples(args)

        with pysam.AlignmentFile(str(out), "rb") as bam:
            read = next(iter(bam))
            assert read.get_tag("RG") == "donor1"

    def test_rg_header_order_is_deterministic(self, barcoded_bam_file,
                                              barcode_donor_tsv, tmp_path,
                                              make_args):
        # RG records are emitted in sorted donor order, not set-iteration order,
        # so the same input always yields the same header.
        out = tmp_path / "assigned.bam"
        args = make_args(
            bam=barcoded_bam_file, tsv=barcode_donor_tsv, output=str(out),
            barcode_column=None, donor_id_column=None,
        )
        assign_samples(args)

        with pysam.AlignmentFile(str(out), "rb") as bam:
            rg_ids = [rg["ID"] for rg in bam.header.to_dict()["RG"]]
        assert rg_ids == sorted(rg_ids)

    def test_rerun_does_not_duplicate_rg_records(self, barcoded_bam_file,
                                                 barcode_donor_tsv, tmp_path,
                                                 make_args):
        # Feeding an already-tagged BAM back in must not append a second copy of
        # each RG record: @RG.ID has to stay unique.
        first = tmp_path / "pass1.bam"
        second = tmp_path / "pass2.bam"
        assign_samples(make_args(
            bam=barcoded_bam_file, tsv=barcode_donor_tsv, output=str(first),
            barcode_column=None, donor_id_column=None,
        ))
        assign_samples(make_args(
            bam=str(first), tsv=barcode_donor_tsv, output=str(second),
            barcode_column=None, donor_id_column=None,
        ))

        with pysam.AlignmentFile(str(second), "rb") as bam:
            rg_ids = [rg["ID"] for rg in bam.header.to_dict()["RG"]]
        assert rg_ids == sorted(set(rg_ids))

    def test_missing_donor_column_raises(self, barcoded_bam_file, tmp_path,
                                         make_args):
        # missing_donor.tsv has 'barcode' and 'sample' but no 'donor_id' column.
        args = make_args(
            bam=barcoded_bam_file, tsv=data_path("missing_donor.tsv"),
            output=str(tmp_path / "out.bam"),
            barcode_column=None, donor_id_column=None,
        )
        with pytest.raises(ValueError, match="No 'donor_id'"):
            assign_samples(args)
