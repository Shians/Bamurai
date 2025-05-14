import pysam
from typing import Dict, Set
from bamurai.utils_samples import (
    parse_barcode_donor_mapping,
    get_barcodes_for_donor,
    ensure_directory_exists,
    get_read_barcode
)

def extract_sample(args):
    """
    Extract reads from a BAM file for a specific donor ID

    Args:
        args: Command-line arguments containing:
            - bam: Path to input BAM file
            - tsv: Path to TSV file mapping barcodes to donor IDs
            - donor_id: Donor ID to extract reads for
            - output: Path to output BAM file
    """
    # Parse barcode-to-donor mapping
    barcode_donor_map = parse_barcode_donor_mapping(args.tsv)

    # Get barcodes for the specified donor
    donor_barcodes = get_barcodes_for_donor(barcode_donor_map, args.donor_id)

    if not donor_barcodes:
        print(f"No barcodes found for donor_id '{args.donor_id}'.")
        return

    print(f"Found {len(donor_barcodes)} barcodes for donor '{args.donor_id}'")

    # Create the output directory if it doesn't exist
    ensure_directory_exists(args.output)

    # Process the BAM file and extract reads for the specified donor
    read_count = 0
    with pysam.AlignmentFile(args.bam, "rb") as input_file:
        # Create output BAM file using the template of the input
        with pysam.AlignmentFile(args.output, "wb", template=input_file) as output_file:
            # Process each read
            for read in input_file:
                # Extract barcode from read
                barcode = get_read_barcode(read)

                # Write the read to output if its barcode matches the donor
                if barcode and barcode in donor_barcodes:
                    output_file.write(read)
                    read_count += 1

    print(f"Extracted {read_count} reads for donor '{args.donor_id}' to {args.output}")
