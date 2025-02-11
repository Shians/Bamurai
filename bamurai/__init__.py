
from .core import (
    Read,
    parse_reads,
    split_read,
    calculate_split,
    read_version
)

from .stats import (
    file_stats,
    is_fastq,
    calc_n50,
    bam_file_stats,
    fastq_file_stats
)

__version__ = read_version()