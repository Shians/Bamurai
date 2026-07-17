#!/usr/bin/env python
"""
Materialise every static input object the test suite consumes into ``tests/data/``.

The suite reads its inputs from committed static files rather than fabricating
them per-test. This script is the single source of truth for those files: it
uses the shared builders in ``conftest.py`` so the static objects are exactly
what the fixtures reference.

Every BAM is accompanied by a plain-text ``.sam`` twin (BAM is binary) so the
records can be read directly, e.g. with ``cat`` or ``less`` and no samtools.

Usage:
    python tests/generate_test_data.py            # writes to tests/data/
    python tests/generate_test_data.py --out DIR  # writes to DIR
"""

import argparse
import gzip
import os
import sys

import pysam

# Ensure the sibling conftest module is importable when run as a script.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from conftest import (  # noqa: E402
    default_fastq_records,
    default_bam_header,
    default_bam_specs,
    barcoded_bam_specs,
    make_segment,
    write_fastq,
    write_bam,
    make_sequence,
    make_qualities,
    qual_ints_to_ascii,
    HTO_BARCODE,
    HTO_UMI,
    HTO_HASHTAG,
    HTO_LEFT_BUFFER,
    DONOR1_BARCODE,
    DONOR2_BARCODE,
)

DEFAULT_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _bam_to_sam(bam_path, sam_path):
    """Write a human-readable SAM text twin of a BAM file."""
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        with pysam.AlignmentFile(sam_path, "w", template=bam) as sam:
            for read in bam:
                sam.write(read)


def _write_bam_with_sam(path, specs, header=None):
    write_bam(path, specs, header=header)
    _bam_to_sam(path, path[:-4] + ".sam")


# ---------------------------------------------------------------------------
# Barcode / donor mapping TSVs
# ---------------------------------------------------------------------------

def write_barcode_donor_tsv(path):
    """The canonical barcode->donor mapping matching :func:`barcoded_bam_specs`."""
    with open(path, "w") as f:
        f.write(
            "barcode\tdonor_id\n"
            f"{DONOR1_BARCODE}\tdonor1\n"
            f"{DONOR2_BARCODE}\tdonor2\n"
        )


def write_mapping_variants(out_dir):
    """TSV variants exercising the flexible column-detection logic."""
    def p(name):
        return os.path.join(out_dir, name)

    # 'cell' column instead of 'barcode' (auto-detected).
    with open(p("cell_donor.tsv"), "w") as f:
        f.write("cell\tdonor_id\nAAA\td1\n")

    # Both 'barcode' and 'cell' present -> ambiguous, must raise.
    with open(p("both_columns.tsv"), "w") as f:
        f.write("barcode\tcell\tdonor_id\nAAA\tXXX\td1\n")

    # No 'barcode'/'cell' column -> must raise.
    with open(p("missing_barcode.tsv"), "w") as f:
        f.write("something\tdonor_id\nAAA\td1\n")

    # No 'donor_id' column -> must raise (also used by assign_samples).
    with open(p("missing_donor.tsv"), "w") as f:
        f.write(f"barcode\tsample\n{DONOR1_BARCODE}\tdonor1\n")

    # Custom column names for explicit --barcode-column/--donor-id-column.
    with open(p("custom_columns.tsv"), "w") as f:
        f.write("bc\tsample\nAAA\td1\n")

    # Purely numeric barcodes and donor IDs. These must survive parsing as
    # strings; inferred as integers they would never match a BAM tag value.
    with open(p("numeric_barcodes.tsv"), "w") as f:
        f.write("barcode\tdonor_id\n1234\t1\n5678\t2\n")


# ---------------------------------------------------------------------------
# FASTQ inputs (valid, malformed, and the unsupported-extension case)
# ---------------------------------------------------------------------------

def write_barcoded_fastq(path):
    """
    A FASTQ whose read names carry a ``BC:Z:`` barcode field.

    Mirrors the input shape ``split_samples`` expects for FASTQ: the barcode is
    parsed from a tab-separated ``BC:Z:`` token in the read name.
    """
    records = [
        (f"r0\tBC:Z:{DONOR1_BARCODE}", "ACGT", "IIII"),  # -> donor1
        (f"r1\tBC:Z:{DONOR2_BARCODE}", "TTTT", "IIII"),  # -> donor2
        ("r2\tBC:Z:UNKNOWNBARCODE0", "GGGG", "IIII"),    # -> unmapped
    ]
    with open(path, "w") as f:
        for name, seq, qual in records:
            f.write(f"@{name}\n{seq}\n+\n{qual}\n")


def write_malformed_fastqs(out_dir):
    """The deliberately broken FASTQ inputs the validator must reject."""
    def p(name):
        return os.path.join(out_dir, name)

    with open(p("bad_header.fastq"), "w") as f:
        f.write("read1\nACGT\n+\nIIII\n")            # missing leading '@'
    with open(p("bad_separator.fastq"), "w") as f:
        f.write("@read1\nACGT\n-\nIIII\n")           # '-' instead of '+'
    with open(p("length_mismatch.fastq"), "w") as f:
        f.write("@read1\nACGT\n+\nII\n")             # qual shorter than seq
    with open(p("invalid_chars.fastq"), "w") as f:
        f.write("@read1\nACGZ\n+\nIIII\n")           # 'Z' not IUPAC
    with open(p("truncated.fastq"), "w") as f:
        f.write("@read1\n")                          # header only
    with open(p("empty.fastq"), "w") as f:
        f.write("")                                  # zero records
    # A valid FASTQ using IUPAC ambiguity codes (must pass validation).
    write_fastq(p("iupac_ok.fastq"), [("r1", "ACGTNRYK", "IIIIIIII")])
    # A non-BAM/FASTQ extension for the dispatcher's unsupported branch.
    with open(p("unsupported.txt"), "w") as f:
        f.write("hello")


# ---------------------------------------------------------------------------
# HTO read pairs for get_hto
# ---------------------------------------------------------------------------

def _hto_pair_records():
    """Return the (R1, R2) sequences for a single 10x-style HTO pair."""
    r1_seq = HTO_BARCODE + HTO_UMI
    r2_seq = HTO_LEFT_BUFFER + HTO_HASHTAG + "AAAAA"
    return r1_seq, r2_seq


def write_hto_pair(r1_path, r2_path, gzipped=False):
    """A single 10x-style HTO read pair (optionally gzipped)."""
    r1_seq, r2_seq = _hto_pair_records()
    opener = gzip.open if gzipped else open
    with opener(r1_path, "wt") as f:
        f.write(f"@read1 1:N:0\n{r1_seq}\n+\n{'I' * len(r1_seq)}\n")
    with opener(r2_path, "wt") as f:
        f.write(f"@read1 2:N:0\n{r2_seq}\n+\n{'I' * len(r2_seq)}\n")


def write_hto_multi(r1_path, r2_path):
    """Two HTO read pairs (readA, readB) with distinct cell barcodes."""
    bc_a, bc_b = "A" * 16, "C" * 16
    r2_body = HTO_LEFT_BUFFER + HTO_HASHTAG + "AAAAA"
    with open(r1_path, "w") as f:
        f.write(f"@readA 1:N:0\n{bc_a + HTO_UMI}\n+\n{'I' * 28}\n")
        f.write(f"@readB 1:N:0\n{bc_b + HTO_UMI}\n+\n{'I' * 28}\n")
    with open(r2_path, "w") as f:
        f.write(f"@readA 2:N:0\n{r2_body}\n+\n{'I' * len(r2_body)}\n")
        f.write(f"@readB 2:N:0\n{r2_body}\n+\n{'I' * len(r2_body)}\n")


# ---------------------------------------------------------------------------
# Crafted BAMs (barcode-tag variants and deliberately incomplete records)
# ---------------------------------------------------------------------------

def write_crafted_bams(out_dir):
    """BAM inputs for the multi-sample and BAM-validation tests."""
    def p(name):
        return os.path.join(out_dir, name)

    # A read carrying only an RX tag. RX is the UMI, not a cell barcode, so this
    # read must stay unassigned even though its value matches donor1's barcode.
    write_bam(p("rx.bam"), [{
        "name": "rx_read_0",
        "sequence": make_sequence(80, seed=0),
        "quals": make_qualities(80),
        "flag": 4,
        "tags": [("RX", DONOR1_BARCODE, "Z")],
    }])

    # A read carrying both CB (donor1) and RX (donor2): CB names the cell, so
    # donor1 must win and the RX value must be ignored.
    write_bam(p("both_tags.bam"), [{
        "name": "both_tags_read",
        "sequence": make_sequence(80, seed=0),
        "quals": make_qualities(80),
        "flag": 4,
        "tags": [("CB", DONOR1_BARCODE, "Z"), ("RX", DONOR2_BARCODE, "Z")],
    }])

    # Two single-read BAMs for the concatenate_bam_files test.
    for i in range(2):
        write_bam(p(f"concat_in_{i}.bam"), [{
            "name": f"r{i}",
            "sequence": make_sequence(20, seed=i),
            "quals": make_qualities(20),
            "flag": 4,
        }])

    # A record with a name but no sequence (validation must reject).
    header = default_bam_header()
    with pysam.AlignmentFile(p("noseq.bam"), "wb", header=header) as bam:
        seg = pysam.AlignedSegment(header)
        seg.query_name = "r0"
        seg.flag = 4
        bam.write(seg)

    # A record with a sequence but absent qualities (validation must reject).
    with pysam.AlignmentFile(p("noqual.bam"), "wb", header=header) as bam:
        seg = pysam.AlignedSegment(header)
        seg.query_name = "r0"
        seg.query_sequence = "ACGTACGT"  # qualities left unset -> stored missing
        seg.flag = 4
        bam.write(seg)


# Human-readable notes written to MANIFEST.txt, keyed by filename.
_MANIFEST = {
    "reads.fastq":
        "5 FASTQ reads, lengths [50, 120, 250, 80, 300].",
    "reads.fastq.gz":
        "Gzip of reads.fastq (identical records).",
    "reads.bam":
        "5 primary reads (same as reads.fastq) PLUS 1 secondary (flag 256) "
        "and 1 supplementary (flag 2048) record that parse_reads must ignore.",
    "reads.sam":
        "Plain-text SAM twin of reads.bam for inspection.",
    "barcoded.bam":
        "5 primary reads carrying CB barcode tags: 2xdonor1, 2xdonor2, "
        "1 unknown barcode (-> unmapped).",
    "barcoded.sam":
        "Plain-text SAM twin of barcoded.bam.",
    "barcoded.fastq":
        "3 FASTQ reads with BC:Z: barcodes in the read name (donor1, donor2, "
        "unknown) for split_samples FASTQ input.",
    "mapping.tsv":
        "Barcode->donor TSV matching barcoded.bam / barcoded.fastq.",
    "cell_donor.tsv":
        "Mapping TSV using a 'cell' column instead of 'barcode'.",
    "both_columns.tsv":
        "Mapping TSV with both 'barcode' and 'cell' (ambiguous -> error).",
    "missing_barcode.tsv":
        "Mapping TSV with no 'barcode'/'cell' column (-> error).",
    "missing_donor.tsv":
        "Mapping TSV with no 'donor_id' column (-> error).",
    "custom_columns.tsv":
        "Mapping TSV with custom 'bc'/'sample' column names.",
    "numeric_barcodes.tsv":
        "Mapping TSV whose barcodes and donor IDs are all digits (must parse "
        "as strings, not integers).",
    "bad_header.fastq":
        "Malformed FASTQ: header missing leading '@'.",
    "bad_separator.fastq":
        "Malformed FASTQ: separator line is '-' not '+'.",
    "length_mismatch.fastq":
        "Malformed FASTQ: quality shorter than sequence.",
    "invalid_chars.fastq":
        "Malformed FASTQ: non-IUPAC character 'Z' in sequence.",
    "truncated.fastq":
        "Malformed FASTQ: header line only, record truncated.",
    "empty.fastq":
        "Empty FASTQ (zero records) - valid input.",
    "iupac_ok.fastq":
        "Valid FASTQ exercising IUPAC ambiguity codes (N,R,Y,K).",
    "unsupported.txt":
        "Non-BAM/FASTQ file for the validate dispatcher's unsupported branch.",
    "hto_R1.fastq":
        "10x R1: 16bp cell barcode + 12bp UMI.",
    "hto_R2.fastq":
        "10x R2: 10bp left buffer + 15bp hashtag oligo.",
    "hto_R1.fastq.gz":
        "Gzip of hto_R1.fastq.",
    "hto_R2.fastq.gz":
        "Gzip of hto_R2.fastq.",
    "hto_multi_R1.fastq":
        "10x R1 with two records (readA, readB) with distinct barcodes.",
    "hto_multi_R2.fastq":
        "10x R2 with two records matching hto_multi_R1.fastq.",
    "rx.bam":
        "Single read carrying only an RX barcode tag (-> donor1).",
    "both_tags.bam":
        "Single read with CB (donor1) and RX (donor2) tags; CB takes precedence.",
    "concat_in_0.bam":
        "Single-read BAM (r0) for the concatenate_bam_files test.",
    "concat_in_1.bam":
        "Single-read BAM (r1) for the concatenate_bam_files test.",
    "noseq.bam":
        "BAM record with a name but no sequence (validation must reject).",
    "noqual.bam":
        "BAM record with a sequence but absent qualities (validation rejects).",
}


def generate(out_dir):
    os.makedirs(out_dir, exist_ok=True)

    def p(name):
        return os.path.join(out_dir, name)

    # Core valid objects, built from the shared conftest builders.
    write_fastq(p("reads.fastq"), default_fastq_records())
    write_fastq(p("reads.fastq.gz"), default_fastq_records(), gzipped=True)
    _write_bam_with_sam(p("reads.bam"), default_bam_specs())
    _write_bam_with_sam(p("barcoded.bam"), barcoded_bam_specs())

    # Mapping TSVs.
    write_barcode_donor_tsv(p("mapping.tsv"))
    write_mapping_variants(out_dir)

    # FASTQ command inputs (valid and malformed).
    write_barcoded_fastq(p("barcoded.fastq"))
    write_malformed_fastqs(out_dir)

    # HTO read pairs.
    write_hto_pair(p("hto_R1.fastq"), p("hto_R2.fastq"))
    write_hto_pair(p("hto_R1.fastq.gz"), p("hto_R2.fastq.gz"), gzipped=True)
    write_hto_multi(p("hto_multi_R1.fastq"), p("hto_multi_R2.fastq"))

    # Crafted BAMs (tag variants and incomplete records).
    write_crafted_bams(out_dir)

    # Manifest describing each artifact.
    with open(p("MANIFEST.txt"), "w") as f:
        f.write("Static input objects for the Bamurai test suite.\n")
        f.write("Regenerate with: python tests/generate_test_data.py\n\n")
        for name in sorted(_MANIFEST):
            f.write(f"{name}\n    {_MANIFEST[name]}\n")

    written = sorted(os.listdir(out_dir))
    print(f"Wrote {len(written)} files to {out_dir}:")
    for name in written:
        size = os.path.getsize(p(name))
        print(f"  {name} ({size} bytes)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", default=DEFAULT_OUT,
        help="Output directory (default: tests/data/)",
    )
    args = parser.parse_args()
    generate(args.out)


if __name__ == "__main__":
    main()
