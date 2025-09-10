from bamurai.utils import is_fastq, create_progress_bar_for_file, count_reads_async_generic
from bamurai.core import parse_reads

def calc_n50(read_lengths):
    """Calculate the N50 statistic for a list of read lengths."""
    if not read_lengths:
        return None
    read_lengths.sort(reverse=True)
    total_bp = sum(read_lengths)
    half_bp = total_bp / 2
    bp_sum = 0
    for read_len in read_lengths:
        bp_sum += read_len
        if bp_sum >= half_bp:
            return read_len

def file_read_stats(read_file):
    """Calculate statistics for a BAM or FASTQ file using parse_reads."""
    read_lengths = []
    total_reads = 0
    
    # Create progress bar
    pbar = create_progress_bar_for_file(read_file, "Calculating statistics")
    count_thread = count_reads_async_generic(read_file, pbar)
    
    for read in parse_reads(read_file):
        read_lengths.append(len(read))
        total_reads += 1
        pbar.update(1)

    pbar.close()

    if not read_lengths:
        return {
            "total_reads": 0,
            "avg_read_len": 0,
            "throughput": 0,
            "n50": 0
        }

    throughput = sum(read_lengths)
    avg_read_len = round(throughput / len(read_lengths))
    n50 = calc_n50(read_lengths)
    return {
        "total_reads": total_reads,
        "avg_read_len": avg_read_len,
        "throughput": throughput,
        "n50": n50
    }

    # fastq_file_stats is now handled by file_read_stats

def file_stats(args):
    stats = file_read_stats(args.reads)

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
