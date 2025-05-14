import pysam
import os
import tempfile
import shutil
import uuid
from collections import defaultdict
from typing import Dict, List, Set, Tuple
from bamurai.utils_samples import (
    parse_barcode_donor_mapping,
    get_read_barcode,
    concatenate_bam_files
)

def split_bam_by_donor(
    input_bam: str,
    barcode_donor_map: Dict[str, str],
    temp_dir: str
) -> Tuple[Set[str], Dict[str, str]]:
    """
    Split a BAM file by donor ID and write to temporary files

    Args:
        input_bam: Path to input BAM file
        barcode_donor_map: Dictionary mapping barcodes to donor IDs
        temp_dir: Directory to write temporary BAM files

    Returns:
        Tuple containing:
        - Set of unique donor IDs
        - Dictionary mapping donor IDs to temporary file paths
    """
    # Create a unique identifier for this BAM file's outputs
    bam_uuid = str(uuid.uuid4())[:8]
    bam_basename = os.path.basename(input_bam).split('.')[0]

    # Identify all unique donor IDs
    unique_donors = set(barcode_donor_map.values())

    # Dictionary to store temp file paths for each donor
    temp_files = {}

    # Open input BAM file
    with pysam.AlignmentFile(input_bam, "rb") as input_file:
        # Create a dictionary to store output files, keyed by donor ID
        output_files = {}

        # Create a temporary output BAM file for each donor
        for donor_id in unique_donors:
            temp_path = os.path.join(temp_dir, f"{bam_basename}_{donor_id}_{bam_uuid}.bam")
            temp_files[donor_id] = temp_path
            output_files[donor_id] = pysam.AlignmentFile(
                temp_path, "wb", template=input_file
            )

        # Create a temp output file for unmapped reads
        unmapped_temp_path = os.path.join(temp_dir, f"{bam_basename}_unmapped_{bam_uuid}.bam")
        temp_files["unmapped"] = unmapped_temp_path
        output_files["unmapped"] = pysam.AlignmentFile(
            unmapped_temp_path, "wb", template=input_file
        )

        # Process reads in the input BAM file
        for read in input_file:
            # Extract barcode from read
            barcode = get_read_barcode(read)

            # Determine which donor the read belongs to
            donor_id = barcode_donor_map.get(barcode, "unmapped") if barcode else "unmapped"

            # Write the read to the appropriate file
            output_files[donor_id].write(read)

        # Close all output files
        for out_file in output_files.values():
            out_file.close()

    print(f"Split {input_bam} into {len(unique_donors)} temporary donor files, plus unmapped reads")
    return unique_donors, temp_files

def split_samples(args):
    # Parse barcode-to-donor mapping
    barcode_donor_map = parse_barcode_donor_mapping(args.tsv)

    # Handle whether we received a list of BAM files or just one
    bam_files = args.bam if isinstance(args.bam, list) else [args.bam]

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Create a temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Created temporary directory: {temp_dir}")

        # Keep track of all donor IDs and their temporary files
        all_donors = set()
        donor_to_temp_files = defaultdict(list)

        # Process each BAM file
        for bam_file in bam_files:
            # Split BAM file by donor, creating temporary files
            donors, temp_files = split_bam_by_donor(
                bam_file,
                barcode_donor_map,
                temp_dir
            )

            # Update our tracking of donors and temp files
            all_donors.update(donors)
            for donor_id, temp_file in temp_files.items():
                donor_to_temp_files[donor_id].append(temp_file)

        # Now concatenate the temporary files for each donor
        for donor_id in all_donors.union({"unmapped"}):
            temp_files_for_donor = donor_to_temp_files[donor_id]
            final_output_path = os.path.join(args.output_dir, f"{donor_id}.bam")

            if len(temp_files_for_donor) == 1:
                # If there's only one file, just copy it
                shutil.copy2(temp_files_for_donor[0], final_output_path)
                print(f"Copied {donor_id} file to {final_output_path}")
            else:
                # Otherwise concatenate the files
                concatenate_bam_files(temp_files_for_donor, final_output_path)

        print("All processing complete. Temporary files will be deleted.")
