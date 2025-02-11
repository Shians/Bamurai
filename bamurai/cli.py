import argparse
import sys
import time
from bamurai.core import *
from bamurai.stats import *
from bamurai.split import *
from bamurai.divide import *
from bamurai import __version__

def main():
    parser = argparse.ArgumentParser(description="Bamurai: A tool for processing BAM files")
    subparsers = parser.add_subparsers(dest="command")

    # Subparser for the "split" command
    parser_split = subparsers.add_parser("split", help="Split reads in a BAM file to a target length")
    parser_split.add_argument("reads", type=str, help="Input reads file (BAM)")
    parser_split.add_argument("-l", "--len_target", type=int, help="Target length for splitting reads")
    parser_split.add_argument("-o", "--output", type=str, nargs='?', help="Output file")
    parser_split.set_defaults(func=split_reads)

    # Subparser for the "stats" command
    parser_stat = subparsers.add_parser("stats", help="Calculate statistics for a BAM or FASTQ(.gz) file")
    parser_stat.add_argument("reads", type=str, help="Input reads file (BAM)")
    parser_stat.add_argument("--tsv", action="store_true", help="Output in TSV format", default=False)
    parser_stat.set_defaults(func=file_stats)

    # Print version if "--version" is passed
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    args = parser.parse_args()

    if args.command:
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()