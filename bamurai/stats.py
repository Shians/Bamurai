import pysam
import gzip
from bamurai.utils import is_fastq

def calc_n50(read_lengths):
    """Calculate the N50 statistic for a list of read lengths."""
    read_lengths.sort(reverse=True)
    total_bp = sum(read_lengths)
    half_bp = total_bp / 2
    bp_sum = 0
    for read_len in read_lengths:
        bp_sum += read_len
        if bp_sum >= half_bp:
            return read_len

def bam_file_stats(bam_file):
    """Calculate statistics for a BAM file."""
    read_lengths = []
    total_reads = 0
    for read in pysam.AlignmentFile(bam_file, "rb", check_sq=False):
        read_lengths.append(read.query_length)
        total_reads += 1

    throughput = sum(read_lengths)
    avg_read_len = round(throughput / len(read_lengths))
    n50 = calc_n50(read_lengths)
    return {
        "total_reads": total_reads,
        "avg_read_len": avg_read_len,
        "throughput": throughput,
        "n50": n50
    }

def fastq_file_stats(fastq_file):
    """Calculate statistics for a FASTQ file."""
    read_lengths = []
    total_reads = 0

    if (fastq_file.endswith(".gz")):
        f = gzip.open(fastq_file, "rt")
    else:
        f = open(fastq_file, "r")

    while True:
        read_id = f.readline().strip()
        if not read_id:
            break
        sequence = f.readline().strip()
        f.readline()  # skip the "+" line
        quality = f.readline().strip()
        read_lengths.append(len(sequence))
        total_reads += 1

    f.close()

    throughput = sum(read_lengths)
    avg_read_len = round(throughput / len(read_lengths))
    n50 = calc_n50(read_lengths)
    return {
        "total_reads": total_reads,
        "avg_read_len": avg_read_len,
        "throughput": throughput,
        "n50": n50
    }

def file_stats(args):
    if args.reads.endswith(".bam"):
        stats = bam_file_stats(args.reads)
    elif is_fastq(args.reads):
        stats = fastq_file_stats(args.reads)

    if args.tsv:
        # print in tsv style
        print(f"file_name\ttotal_reads\tavg_read_len\tthroughput\tn50")
        print(f"{args.reads}\t{stats['total_reads']}\t{stats['avg_read_len']}\t{stats['throughput']}\t{stats['n50']}")
    else:
        print(f"Statistics for {args.reads}:")
        print(f"  Total reads: {stats['total_reads']}")
        print(f"  Average read length: {stats['avg_read_len']}")
        print(f"  Throughput (Gb): {round(stats['throughput'] / 1e9, 2)}")
        print(f"  N50: {stats['n50']}")
