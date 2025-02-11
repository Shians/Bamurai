import os
import pysam
import numpy as np
from dataclasses import dataclass

@dataclass
class Read:
    """Class to represent a FASTQ read."""
    read_id: str
    sequence: str
    quality: str

    def __init__(self, read_id, sequence, quality):
        self.read_id = read_id
        self.sequence = sequence
        self.quality = quality
        self.validate()

    def __len__(self):
        return len(self.sequence)
    
    # validate the read
    def validate(self):
        if len(self.sequence) != len(self.quality):
            raise ValueError(f"Sequence and quality strings must be of equal length. Offending read: {self.read_id}")
        
    def is_valid(self):
        try:
            self.validate()
            return True
        except ValueError:
            return False
        
    def to_fastq(self):
        return f"@{self.read_id}\n{self.sequence}\n+\n{self.quality}"

def qual_to_fastq_numpy(qualities):
    """Convert query_qualities to FASTQ QUAL using NumPy (Best for Large Arrays)."""
    return (np.array(qualities, dtype=np.uint8) + 33).tobytes().decode()

def parse_reads(bam_file):
    """Parse reads from a BAM file."""
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            yield Read(read.query_name, read.query_sequence, qual_to_fastq_numpy(read.query_qualities))

def keep_n_bases(read, n, on = "left"):
    """Trim n bases of a read."""
    read_len = len(read)

    if n > read_len:
        return read
    elif on == "left":
        seq = read.sequence[:n]
        qual = read.quality[:n]
    elif on == "right":
        seq = read.sequence[(read_len - n):]
        qual = read.quality[(read_len - n):]
    
    return Read(read.read_id, seq, qual)
    
def split_read(read, at: list[int]):
    """Split a read at a given positions."""
    reads = []

    if len(at) == 0:
        read.read_id = f'{read.read_id}_0'
        return [read]
    
    count = 0
    start = 0

    for pos in at:
        reads.append(Read(f'{read.read_id}_{count}', read.sequence[start:pos], read.quality[start:pos]))
        start = pos
        count += 1
    
    reads.append(Read(f'{read.read_id}_{count}', read.sequence[start:], read.quality[start:]))

    return reads

def read_version():
    version_file = os.path.join(os.path.dirname(__file__), "..", "VERSION")
    with open(version_file) as vf:
        return vf.read().strip()