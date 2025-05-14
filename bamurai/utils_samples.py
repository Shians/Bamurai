import pysam
import pandas as pd
import os
from typing import Dict, Set, List

def parse_barcode_donor_mapping(tsv_file: str) -> Dict[str, str]:
    """
    Parse a TSV file with barcode and donor_id columns

    Args:
        tsv_file: Path to TSV file containing barcode to donor mapping

    Returns:
        Dictionary mapping barcodes to donor IDs
    """
    df = pd.read_csv(tsv_file, sep='\t')
    if 'barcode' not in df.columns or 'donor_id' not in df.columns:
        raise ValueError("TSV file must contain 'barcode' and 'donor_id' columns")

    return dict(zip(df.barcode, df.donor_id))

def get_barcodes_for_donor(barcode_donor_map: Dict[str, str], donor_id: str) -> Set[str]:
    """
    Get the set of barcodes associated with a specific donor

    Args:
        barcode_donor_map: Dictionary mapping barcodes to donor IDs
        donor_id: The donor ID to extract barcodes for

    Returns:
        Set of barcodes associated with the donor
    """
    return {barcode for barcode, donor in barcode_donor_map.items() if donor == donor_id}

def ensure_directory_exists(file_path: str) -> None:
    """
    Ensures the directory for the given file path exists

    Args:
        file_path: Path to a file
    """
    directory = os.path.dirname(file_path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)

def get_read_barcode(read) -> str:
    """
    Extract barcode from a read by checking common barcode tags

    Args:
        read: A pysam aligned segment (read)

    Returns:
        The barcode string or None if no barcode is found
    """
    for tag in ('CB', 'XC', 'BC'):  # Common barcode tags
        if read.has_tag(tag):
            return read.get_tag(tag)
    return None

def concatenate_bam_files(file_list: List[str], output_path: str) -> None:
    """
    Concatenate multiple BAM files into a single file

    Args:
        file_list: List of BAM files to concatenate
        output_path: Path to write the concatenated BAM file
    """
    if not file_list:
        return

    # Open the first file to use as a template
    with pysam.AlignmentFile(file_list[0], "rb") as template:
        # Create the output file using the template
        with pysam.AlignmentFile(output_path, "wb", template=template) as outfile:
            # Iterate through all input files
            for bam_file in file_list:
                with pysam.AlignmentFile(bam_file, "rb") as infile:
                    # Copy all reads from the input file to the output file
                    for read in infile:
                        outfile.write(read)

    print(f"Concatenated {len(file_list)} files into {output_path}")
