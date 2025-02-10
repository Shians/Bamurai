import pysam
import time
import numpy as np

bam_path = "data/nsc_1_mat.bam"

with pysam.AlignmentFile(bam_path, "rb") as bam:
    start = time.time()
    for read in bam.fetch():
        data = read.query_name  # Minimal access
    print("Fetch qname:", time.time() - start, "seconds")

with pysam.AlignmentFile(bam_path, "rb") as bam:
    start = time.time()
    for read in bam.fetch():
        data = read.query_sequence
    print("Fetch qseq:", time.time() - start, "seconds")

def qual_to_fastq(qualities):
    """Convert query_qualities to a FASTQ QUAL string using bytearray (FAST)."""
    return bytearray(q + 33 for q in qualities).decode()

def qual_to_fastq_numpy(qualities):
    """Convert query_qualities to FASTQ QUAL using NumPy (Best for Large Arrays)."""
    return (np.array(qualities, dtype=np.uint8) + 33).tobytes().decode()

with pysam.AlignmentFile(bam_path, "rb") as bam:
    start = time.time()
    for read in bam.fetch():
        data = qual_to_fastq_numpy(read.query_qualities)
    print("Fetch qual:", time.time() - start, "seconds")
