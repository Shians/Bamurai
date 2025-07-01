import argparse
import textwrap

from bamurai.core import *
from bamurai.stats import *
from bamurai.split import *
from bamurai.divide import *
from bamurai.validate import *
from bamurai.chunk import *
from bamurai.split_samples import *
from bamurai.extract_sample import *
from bamurai.assign_samples import *
from bamurai.get_hto import *
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
        Divide reads in a BAM/FASTQ file into a fixed number of fragments. The output will be in FASTQ format written to the output file specified. If no output file is defined then the output is written to stdout. Reads that are shorter than the minimum length are not divided.
        """,
        formatter_class=CustomFormatter
    )
    parser_divide.add_argument("reads", type=str, help=input_read_arg_description)
    parser_divide.add_argument("-n", "--num_fragments", type=int, help="Number of fragments to divide reads into (default = 2)", default=2)
    parser_divide.add_argument("-m", "--min_length", type=int, help="Minimum length for a fragment, reads will not be divided if resultant length is less than this (default = 100)", default=100)
    parser_divide.add_argument("-o", "--output", type=str, nargs='?', help=output_file_arg_description)
    parser_divide.set_defaults(func=divide_reads)

    # Subparser for the "validate" command
    parser_validate = subparsers.add_parser(
        "validate",
        help="Validate a BAM or FASTQ file",
        description = """
        Validate a BAM or FASTQ file to ensure it is correctly formatted. The validation checks include:

        For BAM files:
        - Check that the header is present and correctly formatted
        - Check that each record has a query name, sequence, and quality scores
        - Check that the sequence and quality lengths are equal

        For FASTQ files:
        - Check that the header starts with '@'
        - Check that the separator line starts with '+'
        - Check that the sequence and quality lengths are equal
        - Check that the sequence contains only valid IUPAC characters
        """,
        formatter_class=CustomFormatter
    )
    parser_validate.add_argument("reads", type=str, help=input_read_arg_description)
    parser_validate.set_defaults(func=validate_file)

    # Subparser for the "chunk" command
    parser_chunk = subparsers.add_parser(
        "chunk",
        help="Split BAM/FASTQ file into chunks of specified size",
        description="""
        Split BAM/FASTQ file into chunks of at least the specified size. Output files
        will be named <prefix>_1.fastq, <prefix>_2.fastq, etc. Each chunk will be at
        least as large as the specified size, but may be larger to avoid splitting
        individual reads.
        """,
        formatter_class=CustomFormatter
    )
    parser_chunk.add_argument("reads", type=str, help=input_read_arg_description)
    parser_chunk.add_argument(
        "-s", "--size",
        type=str,
        required=True,
        help="Minimum chunk size (e.g. 1G, 100M, 1000K)"
    )
    parser_chunk.add_argument(
        "-p", "--prefix",
        type=str,
        default="chunk",
        help="Output file prefix (default: 'chunk')"
    )
    parser_chunk.set_defaults(func=chunk_reads)

    # Subparser for the "split_samples" command
    split_parser = subparsers.add_parser(
        "split_samples",
        help="Split BAM or FASTQ file by donor ID"
    )
    split_parser.add_argument("--input", required=True, help="Input BAM or FASTQ file(s)", nargs='+')
    split_parser.add_argument("--tsv", required=True, help="two-column TSV file mapping barcodes to donor IDs, with headers 'barcode' and 'donor_id'")
    split_parser.add_argument("--output-dir", default="output", help="Output directory for split files (default: 'output')")
    split_parser.set_defaults(func=split_samples)

    # Subparser for the "extract_sample" command
    extract_parser = subparsers.add_parser(
        "extract_sample",
        help="Extract all reads for a given donor_id from BAM file(s) using a barcode-to-donor mapping"
    )
    extract_parser.add_argument("--bam", required=True, help="Input BAM file(s)", nargs='+')
    extract_parser.add_argument("--tsv", required=True, help="TSV file mapping barcodes to donor IDs, with headers 'barcode' and 'donor_id'")
    extract_parser.add_argument("--donor-id", required=True, help="Donor ID to extract")
    extract_parser.add_argument("--output", required=True, help="Output BAM file for extracted reads")
    extract_parser.set_defaults(func=extract_sample)

    # Subparser for the "assign_sample" command
    assign_parser = subparsers.add_parser(
        "assign_samples",
        help="Assign donor_id to RG tag in BAM file using barcode-to-donor mapping"
    )
    assign_parser.add_argument("--bam", required=True, help="Input BAM file")
    assign_parser.add_argument("--tsv", required=True, help="TSV file mapping barcodes to donor IDs, with headers 'barcode' and 'donor_id'")
    assign_parser.add_argument("--output", required=True, help="Output BAM file with RG tags assigned")
    assign_parser.set_defaults(func=assign_samples)

    # Subparser for the "get_hto" command
    parser_get_hto = subparsers.add_parser(
        "get_hto",
        help="Extract HTO information from 10x FASTQ files. ",
        description="""
        Extract HTO (Hashtag Oligo) information from 10x FASTQ files. Assumes that
        the first read (R1) contains the cell barcode and UMI, and the second read
        (R2) contains the HTO sequence. The output will be in a tab-separated
        format with columns: read_name, cell_barcode, umi, hto.
        """,
        formatter_class=CustomFormatter
    )
    parser_get_hto.add_argument("--r1", required=True, help="FASTQ R1 file.")
    parser_get_hto.add_argument("--r2", required=True, help="FASTQ R2 file.")
    parser_get_hto.add_argument("--bc-len", type=int, required=True, help="Cell barcode length.")
    parser_get_hto.add_argument("--umi-len", type=int, required=True, help="UMI length.")
    parser_get_hto.add_argument("--output", type=str, required=True, help="Output file for HTO information.")
    parser_get_hto.add_argument("--hashtag-len", type=int, default=15, help="Hashtag length (default: 15).")
    parser_get_hto.add_argument("--hashtag-left-buffer", type=int, default=10, help="Hashtag left buffer (default: 10).")
    parser_get_hto.set_defaults(func=get_hto)

    # Print version if "--version" is passed
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    args = parser.parse_args()

    if args.command:
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
