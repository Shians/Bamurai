import pysam
import os
import tempfile
from typing import Dict, Set, List
from bamurai.utils_samples import (
    parse_barcode_donor_mapping,
    get_barcodes_for_donor,
    ensure_directory_exists,
    get_read_barcode,
    concatenate_bam_files
)

def extract_reads_from_bam(
    input_bam: str,
    donor_barcodes: Set[str],
    output_file: str = None,
    temp_dir: str = None
) -> tuple[int, str]:
    """
    Extract reads from a BAM file for a specific donor's barcodes

    Args:
        input_bam: Path to input BAM file
        donor_barcodes: Set of barcodes associated with the donor
        output_file: Path to output BAM file (if None, create a temp file)
        temp_dir: Directory for temporary files

    Returns:
        tuple containing:
        - Count of extracted reads
        - Path to output file
    """
    if output_file is None:
        # Create a temporary output file in the temp directory
        output_file = os.path.join(temp_dir, f"{os.path.basename(input_bam)}.temp")

    read_count = 0
    with pysam.AlignmentFile(input_bam, "rb") as input_file:
        # Create output BAM file using the template of the input
        with pysam.AlignmentFile(output_file, "wb", template=input_file) as output_file_handle:
            # Process each read
            for read in input_file:
                # Extract barcode from read
                barcode = get_read_barcode(read)

                # Write the read to output if its barcode matches the donor
                if barcode and barcode in donor_barcodes:
                    output_file_handle.write(read)
                    read_count += 1

    return read_count, output_file

def extract_sample(args):
    """
    Extract reads from BAM file(s) for a specific donor ID

    Args:
        args: Command-line arguments containing:
            - bam: Path(s) to input BAM file(s)
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

    # Handle whether we received a list of BAM files or just one
    bam_files = args.bam if isinstance(args.bam, list) else [args.bam]

    # Create the output directory if it doesn't exist
    ensure_directory_exists(args.output)

    # Process single BAM file directly to output
    if len(bam_files) == 1:
        read_count, _ = extract_reads_from_bam(
            bam_files[0],
            donor_barcodes,
            output_file=args.output
        )
        print(f"Extracted {read_count} reads for donor '{args.donor_id}' to {args.output}")
        return

    # Process multiple BAM files using temp files
    temp_files = []
    total_read_count = 0

    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Processing {len(bam_files)} BAM files...")

        # Extract reads from each BAM file to temp files
        for bam_file in bam_files:
            read_count, temp_file = extract_reads_from_bam(
                bam_file,
                donor_barcodes,
                temp_dir=temp_dir
            )
            temp_files.append(temp_file)
            total_read_count += read_count
            print(f"Extracted {read_count} reads from {bam_file}")

        # Concatenate all temp files into final output
        concatenate_bam_files(temp_files, args.output)

    print(f"Total: Extracted {total_read_count} reads for donor '{args.donor_id}' to {args.output}")
