import pysam
import pandas as pd
import os
import argparse
from collections import defaultdict
from typing import Dict, List, Set

def parse_barcode_donor_mapping(csv_file: str) -> Dict[str, str]:
    """
    Parse a CSV file with barcode and donor_id columns

    Args:
        csv_file: Path to CSV file containing barcode to donor mapping

    Returns:
        Dictionary mapping barcodes to donor IDs
    """
    df = pd.read_csv(csv_file)
    if 'barcode' not in df.columns or 'donor_id' not in df.columns:
        raise ValueError("CSV file must contain 'barcode' and 'donor_id' columns")

    return dict(zip(df.barcode, df.donor_id))

def split_bam_by_donor(
    input_bam: str,
    barcode_donor_map: Dict[str, str],
    output_dir: str
) -> None:
    """
    Split a BAM file by donor ID

    Args:
        input_bam: Path to input BAM file
        barcode_donor_map: Dictionary mapping barcodes to donor IDs
        output_dir: Directory to write output BAM files
    """

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Identify all unique donor IDs
    unique_donors = set(barcode_donor_map.values())

    # Open input BAM file
    with pysam.AlignmentFile(input_bam, "rb") as input_file:
        # Create a dictionary to store output files, keyed by donor ID
        output_files = {}

        # Create an output BAM file for each donor
        for donor_id in unique_donors:
            output_path = os.path.join(output_dir, f"{donor_id}.bam")
            output_files[donor_id] = pysam.AlignmentFile(
                output_path, "wb", template=input_file
            )

        # Create an output file for unmapped reads
        unmapped_path = os.path.join(output_dir, "unmapped.bam")
        output_files["unmapped"] = pysam.AlignmentFile(
            unmapped_path, "wb", template=input_file
        )

        # Process reads in the input BAM file
        for read in input_file:
            # Extract barcode from read (assuming it's stored in a tag, e.g., 'CB')
            barcode = None
            if read.has_tag('CB'):
                barcode = read.get_tag('CB')

            # Determine which donor the read belongs to
            donor_id = barcode_donor_map.get(barcode, "unmapped") if barcode else "unmapped"

            # Write the read to the appropriate file
            output_files[donor_id].write(read)

        # Close all output files
        for out_file in output_files.values():
            out_file.close()

    print(f"Split BAM into {len(unique_donors)} donor files, plus unmapped reads")

def split_samples(args):
    # Parse barcode-to-donor mapping
    barcode_donor_map = parse_barcode_donor_mapping(args.csv)

    # Split BAM file by donor
    split_bam_by_donor(
        args.bam,
        barcode_donor_map,
        args.output_dir
    )
