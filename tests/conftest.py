"""
Shared pytest fixtures and synthetic data builders for the Bamurai test suite.

The suite reads its inputs from committed **static** objects under
``tests/data/`` rather than fabricating them per-test. The builder functions
here are the shared machinery ``tests/generate_test_data.py`` uses to write
those static files, so what the fixtures point at is exactly what the builders
produce. Everything is deterministic so assertions can be exact.

Key facts about the data model that the builders honour:
  * ``parse_reads`` only yields *primary* alignments for BAM files, so the
    builders deliberately mix in secondary/supplementary records that the
    tests then assert are skipped.
  * pysam stores qualities as integer arrays; FASTQ stores them as ASCII
    (Phred + 33). ``qual_to_fastq_numpy`` bridges the two.
"""

import gzip
import os
from array import array
from types import SimpleNamespace

import pysam
import pytest


# ---------------------------------------------------------------------------
# Static data location and shared constants
# ---------------------------------------------------------------------------

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def data_path(name):
    """Absolute path to a static input object in ``tests/data/``."""
    return os.path.join(DATA_DIR, name)


# Known barcode/HTO constants shared by the static-data generator and the tests
# that assert against those files. Keep in sync with generate_test_data.py.
DONOR1_BARCODE = "AAACCCAAAGGGTTT"
DONOR2_BARCODE = "CCCGGGTTTAAACCC"

HTO_BARCODE = "AAACCCAAAGGGTTTC"   # 16bp cell barcode
HTO_UMI = "ACGTACGTACGT"           # 12bp UMI
HTO_HASHTAG = "GATCGATCGATCGAT"    # 15bp hashtag oligo
HTO_LEFT_BUFFER = "N" * 10         # 10bp left buffer in R2


# ---------------------------------------------------------------------------
# Low level sequence helpers
# ---------------------------------------------------------------------------

_BASES = "ACGT"


def make_sequence(length, seed=0):
    """Return a deterministic DNA sequence of the requested length."""
    return "".join(_BASES[(i * 7 + seed * 3) % 4] for i in range(length))


def make_qualities(length, value=30):
    """Return a list of integer Phred quality scores of the requested length."""
    return [value] * length


def qual_ints_to_ascii(quals):
    """Convert a list of integer Phred scores to a FASTQ quality string."""
    return "".join(chr(q + 33) for q in quals)


# ---------------------------------------------------------------------------
# FASTQ builders
# ---------------------------------------------------------------------------

def write_fastq(path, records, gzipped=False):
    """
    Write ``records`` to ``path`` as FASTQ.

    Each record is a ``(read_id, sequence, quality_string)`` tuple. ``read_id``
    is written verbatim after the leading ``@`` (no extra whitespace added).
    """
    opener = gzip.open if gzipped else open
    with opener(path, "wt") as handle:
        for read_id, sequence, quality in records:
            handle.write(f"@{read_id}\n{sequence}\n+\n{quality}\n")
    return str(path)


def default_fastq_records():
    """A small, known set of FASTQ records used across several tests."""
    records = []
    for i, length in enumerate([50, 120, 250, 80, 300]):
        seq = make_sequence(length, seed=i)
        qual = qual_ints_to_ascii(make_qualities(length, value=30 + i))
        records.append((f"read_{i}", seq, qual))
    return records


# ---------------------------------------------------------------------------
# BAM builders
# ---------------------------------------------------------------------------

def default_bam_header():
    """A minimal but complete BAM header with a single reference."""
    return pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "unsorted"},
            "SQ": [{"SN": "chr1", "LN": 100000}],
        }
    )


def make_segment(header, name, sequence, quals, flag=4, tags=None):
    """
    Build a single ``pysam.AlignedSegment``.

    Reads are unmapped by default (``flag=4``) which keeps the builders simple
    while still exercising the sequence/quality/tag paths Bamurai cares about.
    Use ``flag`` to set the secondary (256) or supplementary (2048) bits.
    """
    segment = pysam.AlignedSegment(header)
    segment.query_name = name
    segment.query_sequence = sequence
    # query_qualities must be assigned *after* query_sequence, which resets it.
    segment.query_qualities = array("B", quals)
    segment.flag = flag
    segment.reference_id = -1
    segment.reference_start = -1
    if tags:
        for tag, value, value_type in tags:
            segment.set_tag(tag, value, value_type=value_type)
    return segment


def write_bam(path, segments_spec, header=None):
    """
    Write a BAM file from a list of segment specs.

    Each spec is a dict accepted by :func:`make_segment` (minus ``header``):
    ``{"name", "sequence", "quals", "flag"(optional), "tags"(optional)}``.
    """
    header = header or default_bam_header()
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for spec in segments_spec:
            bam.write(
                make_segment(
                    header,
                    spec["name"],
                    spec["sequence"],
                    spec["quals"],
                    flag=spec.get("flag", 4),
                    tags=spec.get("tags"),
                )
            )
    return str(path)


def default_bam_specs():
    """
    Five primary reads plus one secondary and one supplementary record.

    The secondary/supplementary records share names and sequences with primary
    reads so tests can confirm they are filtered rather than counted.
    """
    specs = []
    for i, length in enumerate([50, 120, 250, 80, 300]):
        seq = make_sequence(length, seed=i)
        specs.append(
            {
                "name": f"read_{i}",
                "sequence": seq,
                "quals": make_qualities(length, value=30 + i),
                "flag": 4,  # unmapped primary
            }
        )
    # A secondary alignment (flag bit 256) that must be ignored.
    specs.append(
        {
            "name": "read_0",
            "sequence": make_sequence(50, seed=0),
            "quals": make_qualities(50),
            "flag": 256,
        }
    )
    # A supplementary alignment (flag bit 2048) that must be ignored.
    specs.append(
        {
            "name": "read_2",
            "sequence": make_sequence(100, seed=2),
            "quals": make_qualities(100),
            "flag": 2048,
        }
    )
    return specs


def barcoded_bam_specs():
    """
    Primary reads carrying CB barcode tags, for multi-sample command tests.

    Barcodes map to donors as: AAA*->donor1, CCC*->donor2, and one read with an
    unknown barcode that should land in the ``unmapped`` bucket.
    """
    barcodes = [
        DONOR1_BARCODE,
        DONOR1_BARCODE,
        DONOR2_BARCODE,
        DONOR2_BARCODE,
        "GGGTTTAAACCCGGG",  # unknown barcode -> unmapped
    ]
    specs = []
    for i, barcode in enumerate(barcodes):
        length = 100 + i * 10
        specs.append(
            {
                "name": f"bc_read_{i}",
                "sequence": make_sequence(length, seed=i),
                "quals": make_qualities(length),
                "flag": 4,
                "tags": [("CB", barcode, "Z")],
            }
        )
    return specs


# ---------------------------------------------------------------------------
# Fixtures: concrete files
# ---------------------------------------------------------------------------

@pytest.fixture
def fastq_file():
    """A plain FASTQ file with the default record set."""
    return data_path("reads.fastq")


@pytest.fixture
def fastq_gz_file():
    """A gzipped FASTQ file with the default record set."""
    return data_path("reads.fastq.gz")


@pytest.fixture
def bam_file():
    """A BAM file with 5 primary + 1 secondary + 1 supplementary records."""
    return data_path("reads.bam")


@pytest.fixture
def barcoded_bam_file():
    """A BAM file whose reads carry CB barcode tags."""
    return data_path("barcoded.bam")


@pytest.fixture
def barcode_donor_tsv():
    """A barcode->donor mapping TSV matching :func:`barcoded_bam_specs`."""
    return data_path("mapping.tsv")


# ---------------------------------------------------------------------------
# Fixtures: helpers
# ---------------------------------------------------------------------------

@pytest.fixture
def make_args():
    """Factory building an argparse-like namespace from keyword arguments."""
    def _make(**kwargs):
        return SimpleNamespace(**kwargs)

    return _make
