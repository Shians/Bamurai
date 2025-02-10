import argparse
import sys
import time
from bamurai.core import *

def split_reads(args):
    # record runtime
    print("Running Bamurai split...", file=sys.stderr)
    start_time = time.time()

    # pretty print the arguments in arg: value format
    arg_desc_dict = {
        "reads": "Input file",
        "len_target": "Target length",
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
        split_locs = calculate_split(read, target_len = args.len_target)
        split = split_read(read, at = split_locs)
        for read in split:
            read_lens.append(len(read))
            
            if args.output:
                f.write(read.to_fastq())
                f.write("\n")
            else:
                print(read.to_fastq())

    if args.output:
        f.close()

    avg_read_len = round(sum(read_lens) / len(read_lens))
    print(f"Average split read length: {avg_read_len}", file=sys.stderr)
    print(f"Time taken: {round(time.time() - start_time, 2)} seconds", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description="Bamurai: A tool for processing BAM files")
    subparsers = parser.add_subparsers(dest="command")

    # Subparser for the "split" command
    parser_split = subparsers.add_parser("split", help="Split reads in a BAM file to a target length")
    parser_split.add_argument("reads", type=str, help="Input reads file (BAM)")
    parser_split.add_argument("-l", "--len_target", type=int, help="Target length for splitting reads")
    parser_split.add_argument("-o", "--output", type=str, nargs='?', help="Output file")
    parser_split.set_defaults(func=split_reads)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()