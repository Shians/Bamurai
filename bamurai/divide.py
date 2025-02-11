import sys
import time
from bamurai.core import *
from bamurai.utils import print_elapsed_time_pretty

def calculate_split_pieces(read, num_pieces: int, min_length: int = 0):
    """Calculate split locations for a read given a number of pieces"""
    if len(read)/num_pieces < min_length:
        return []

    # find the number of splits
    split_size = len(read) // num_pieces

    return [i * split_size for i in range(1, num_pieces)]

def divide_reads(args):
    print("Running Bamurai divide...", file=sys.stderr)
    start_time = time.time()

    total_input_reads = 0
    total_output_reads = 0
    total_unsplit_reads = 0
    # pretty print the arguments in arg: value format
    arg_desc_dict = {
        "reads": "Input file",
        "num_pieces": "Number of pieces",
        "min_length": "Minimum length",
        "output": "Output file"
    }
    print("Arguments:", file=sys.stderr)
    for arg, value in vars(args).items():
        if arg in arg_desc_dict:
            print(f"  {arg_desc_dict[arg]}: {value}", file=sys.stderr)

    # Read the input reads file
    read_lens = []

    # clear the output file
    if args.output:
        f = open(args.output, "w")

    for read in parse_reads(args.reads):
        total_input_reads += 1
        split_locs = calculate_split_pieces(read, num_pieces=args.num_fragments, min_length=args.min_length)
        split = split_read(read, at = split_locs)

        if len(split) == 1:
            total_unsplit_reads += 1

        for read in split:
            total_output_reads += 1
            read_lens.append(len(read))

            if args.output:
                f.write(read.to_fastq())
                f.write("\n")
            else:
                print(read.to_fastq())

    if args.output:
        f.close()

    avg_read_len = round(sum(read_lens) / len(read_lens))
    print(f"Total input reads: {total_input_reads}", file=sys.stderr)
    print(f"Total output reads: {total_output_reads}", file=sys.stderr)
    print(f"Total unsplit reads: {total_unsplit_reads}", file=sys.stderr)
    print(f"Average split read length: {avg_read_len}", file=sys.stderr)
    print_elapsed_time_pretty(start_time)
