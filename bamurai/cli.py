import argparse
import textwrap
from bamurai.core import *
from bamurai.stats import *
from bamurai.split import *
from bamurai.divide import *
from bamurai import __version__

def main():
    parser = argparse.ArgumentParser(description="Bamurai: A tool for processing BAM/FASTQ files")
    subparsers = parser.add_subparsers(dest="command")

    class CustomFormatter(argparse.RawDescriptionHelpFormatter):
        """Custom formatter to wrap text and remove common leading whitespace"""
        def __init__(self, prog, indent_increment=2, max_help_position=24, width=80):
            super().__init__(prog, indent_increment, max_help_position, width)

        def _fill_text(self, text, width, indent):
            # Remove common leading whitespace using textwrap.dedent
            text = textwrap.dedent(text)

            # Split text into lines and wrap each line individually
            lines = []
            for line in text.split('\n'):
                # Preserve empty lines and list items
                if not line.strip() or line.lstrip().startswith('-'):
                    lines.append(line)
                else:
                    # Wrap normal text to width
                    wrapped = textwrap.fill(line, width=width-len(indent))
                    lines.append(wrapped)

            return indent + '\n'.join(lines)

    input_read_arg_description = "Input reads file (BAM/FASTQ)"
    output_file_arg_description = "Output file (FASTQ)"

    # Subparser for the "split" command
    parser_split = subparsers.add_parser(
        "split",
        help="Split reads in a BAM/FASTQ file to a target length",
        description = """
        Split reads in a BAM/FASTQ file to a target length. Each read will be split into fragments as close to the target length as possible. The output will be in FASTQ format written to the output file specified. If no output file is defined then the otuput is written to stdout. Reads that are shorter than the target length are not split.
        """,
        formatter_class=CustomFormatter
    )
    parser_split.add_argument("reads", type=str, help=input_read_arg_description)
    parser_split.add_argument("-l", "--len_target", type=int, help="Target length for splitting reads")
    parser_split.add_argument("-o", "--output", type=str, nargs='?', help=output_file_arg_description)
    parser_split.set_defaults(func=split_reads)

    # Subparser for the "stats" command
    parser_stat = subparsers.add_parser(
        "stats",
        help="Calculate statistics for a BAM or FASTQ(.gz) file",
        description = """
        Calculate statistics for a BAM or FASTQ(.gz) file. The statistics include:

        - Total number of reads
        - Average read length
        - Total throughput (in gigabases)
        - N50 read length
        """,
        formatter_class=CustomFormatter
    )
    parser_stat.add_argument("reads", type=str, help=input_read_arg_description)
    parser_stat.add_argument("--tsv", action="store_true", help="Output in TSV format", default=False)
    parser_stat.set_defaults(func=file_stats)

    # Subparser for the "divide" command
    parser_divide = subparsers.add_parser(
        "divide",
        help="Divide reads in a BAM/FASTQ into fixed number of fragments",
        description = """
        Divide reads in a BAM/FASTQ file into a fixed number of fragments. The output will be in FASTQ format written to the output file specified. If no output file is defined then the otuput is written to stdout. Reads that are shorter than the minimum length are not divided.
        """,
        formatter_class=CustomFormatter
    )
    parser_divide.add_argument("reads", type=str, help=input_read_arg_description)
    parser_divide.add_argument("-n", "--num_fragments", type=int, help="Number of fragments to divide reads into (default = 2)", default=2)
    parser_divide.add_argument("-m", "--min_length", type=int, help="Minimum length for a fragment, reads will not be divided if resultant length is less than this (default = 100)", default=100)
    parser_divide.add_argument("-o", "--output", type=str, nargs='?', help=output_file_arg_description)
    parser_divide.set_defaults(func=divide_reads)

    # Print version if "--version" is passed
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    args = parser.parse_args()

    if args.command:
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
