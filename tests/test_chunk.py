"""Tests for bamurai.chunk: size-based file chunking."""

import glob
import os

from bamurai.core import parse_reads
from bamurai.chunk import chunk_reads, do_chunk


class TestChunkReads:
    def test_size_suffix_parsing(self, fastq_file, tmp_path, make_args):
        # 1K minimum chunk; exercises the 'K' unit path. The default reads are
        # ~1.6 KB on disk, so this rolls into a second chunk after the threshold.
        prefix = str(tmp_path / "chunk")
        chunk_reads(make_args(reads=fastq_file, size="1K", prefix=prefix))
        assert os.path.exists(f"{prefix}_1.fastq")

    def test_megabyte_suffix_single_chunk(self, fastq_file, tmp_path, make_args):
        # 1M dwarfs the ~1.6 KB input -> a single chunk. Exercises the 'M' unit.
        prefix = str(tmp_path / "chunkM")
        chunk_reads(make_args(reads=fastq_file, size="1M", prefix=prefix))
        assert os.path.exists(f"{prefix}_1.fastq")
        assert not os.path.exists(f"{prefix}_2.fastq")

    def test_gigabyte_suffix_single_chunk(self, fastq_file, tmp_path, make_args):
        # 1G dwarfs the input -> a single chunk. Exercises the 'G' unit.
        prefix = str(tmp_path / "chunkG")
        chunk_reads(make_args(reads=fastq_file, size="1G", prefix=prefix))
        assert os.path.exists(f"{prefix}_1.fastq")
        assert not os.path.exists(f"{prefix}_2.fastq")

    def test_bare_integer_size(self, fastq_file, tmp_path, make_args):
        # No unit suffix -> the size is bytes. 100 bytes forces multiple chunks.
        prefix = str(tmp_path / "chunkN")
        chunk_reads(make_args(reads=fastq_file, size="100", prefix=prefix))
        chunks = sorted(glob.glob(f"{prefix}_*.fastq"))
        assert len(chunks) > 1
        total = sum(len(list(parse_reads(c))) for c in chunks)
        assert total == 5

    def test_multiple_chunks_created(self, fastq_file, tmp_path):
        prefix = str(tmp_path / "chunk")
        # Very small chunk size forces a new file after (nearly) every read.
        do_chunk(fastq_file, chunk_size=100, output_prefix=prefix)
        chunks = sorted(glob.glob(f"{prefix}_*.fastq"))
        assert len(chunks) > 1

    def test_all_reads_preserved_across_chunks(self, fastq_file, tmp_path):
        prefix = str(tmp_path / "chunk")
        do_chunk(fastq_file, chunk_size=100, output_prefix=prefix)
        chunks = sorted(glob.glob(f"{prefix}_*.fastq"))
        total = sum(len(list(parse_reads(c))) for c in chunks)
        assert total == 5

    def test_bam_input_chunked_to_fastq(self, bam_file, tmp_path):
        prefix = str(tmp_path / "chunk")
        do_chunk(bam_file, chunk_size=100, output_prefix=prefix)
        chunks = sorted(glob.glob(f"{prefix}_*.fastq"))
        total = sum(len(list(parse_reads(c))) for c in chunks)
        assert total == 5  # primary reads only
