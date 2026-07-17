import pysam
import pandas as pd
import os
from typing import Dict, Set, List

def parse_barcode_donor_mapping(tsv_file: str, barcode_column: str = None, donor_id_column: str = None) -> Dict[str, str]:
    """
    Parse a TSV file with barcode and donor_id columns

    Args:
        tsv_file: Path to TSV file containing barcode to donor mapping
        barcode_column: Optional column name for barcodes (auto-detects 'barcode' or 'cell' if not provided)
        donor_id_column: Optional column name for donor IDs (defaults to 'donor_id')

    Returns:
        Dictionary mapping barcodes to donor IDs
    """
    df = pd.read_csv(tsv_file, sep='\t')
    columns = df.columns.tolist()

    # Auto-detect barcode column
    if barcode_column is None:
        if 'barcode' in columns and 'cell' in columns:
            raise ValueError("Both 'barcode' and 'cell' columns are present in the TSV. Please specify --barcode-column explicitly.")
        elif 'barcode' in columns:
            barcode_column = 'barcode'
        elif 'cell' in columns:
            barcode_column = 'cell'
        else:
            raise ValueError("No 'barcode' or 'cell' column found in the TSV. Please specify --barcode-column.")
    else:
        if barcode_column not in columns:
            raise ValueError(f"Specified barcode column '{barcode_column}' not found in TSV. Available columns: {', '.join(columns)}")

    # Auto-detect donor_id column
    if donor_id_column is None:
        if 'donor_id' in columns:
            donor_id_column = 'donor_id'
        else:
            raise ValueError("No 'donor_id' column found in the TSV. Please specify --donor-id-column.")
    else:
        if donor_id_column not in columns:
            raise ValueError(f"Specified donor_id column '{donor_id_column}' not found in TSV. Available columns: {', '.join(columns)}")

    # Both are compared against, or written into, string BAM tags. Without the
    # cast pandas infers numeric-looking barcodes and donor IDs as integers,
    # which then never match the tag values read off a read.
    return dict(zip(df[barcode_column].astype(str), df[donor_id_column].astype(str)))

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

# Cell barcode tags in order of preference. Per the SAM tag spec CB is the
# optionally-corrected cell barcode and CR its uncorrected form. XC is not in the
# spec (X* is reserved for local use) but is the Drop-seq/early-10x convention.
# BC and RX are deliberately absent: BC identifies a sample and RX a molecule
# (the UMI), so neither names the cell a donor is assigned to.
CELL_BARCODE_TAGS = ('CB', 'CR', 'XC')

def get_read_barcode(read) -> str | None:
    """
    Extract the cell barcode from a read, checking CB, then CR, then XC

    Args:
        read: A pysam aligned segment (read)

    Returns:
        The barcode string or None if no barcode is found
    """
    for tag in CELL_BARCODE_TAGS:
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
