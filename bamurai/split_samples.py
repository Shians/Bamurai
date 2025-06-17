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
from bamurai.utils import is_fastq

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

def split_fastq_by_donor(
    input_fastq: str,
    barcode_donor_map: Dict[str, str],
    temp_dir: str
) -> Tuple[Set[str], Dict[str, str]]:
    """
    Split a FASTQ file by donor ID and write to temporary files.

    Args:
        input_fastq: Path to input FASTQ file
        barcode_donor_map: Dictionary mapping barcodes to donor IDs
        temp_dir: Directory to write temporary FASTQ files

    Returns:
        Tuple containing:
        - Set of unique donor IDs
        - Dictionary mapping donor IDs to temporary file paths
    """
    import gzip
    # Create a unique identifier for this FASTQ file's outputs
    fastq_uuid = str(uuid.uuid4())[:8]
    fastq_basename = os.path.basename(input_fastq).split('.')[0]

    # Identify all unique donor IDs
    unique_donors = set(barcode_donor_map.values())

    # Dictionary to store temp file paths for each donor
    temp_files = {}
    output_files = {}

    # Create a temporary output FASTQ file for each donor
    for donor_id in unique_donors:
        temp_path = os.path.join(temp_dir, f"{fastq_basename}_{donor_id}_{fastq_uuid}.fastq")
        temp_files[donor_id] = temp_path
        output_files[donor_id] = open(temp_path, "wt")

    # Create a temp output file for unmapped reads
    unmapped_temp_path = os.path.join(temp_dir, f"{fastq_basename}_unmapped_{fastq_uuid}.fastq")
    temp_files["unmapped"] = unmapped_temp_path
    output_files["unmapped"] = open(unmapped_temp_path, "wt")

    # Open input FASTQ file (support gzipped files)
    open_func = gzip.open if input_fastq.endswith('.gz') else open
    with open_func(input_fastq, 'rt') as infile:
        while True:
            # Read 4 lines for each FASTQ record
            lines = [infile.readline() for _ in range(4)]
            if not lines[0]:
                break  # EOF
            # Extract barcode from read name (first line)
            readname = lines[0].strip()
            barcode = None
            # Look for tab-separated fields in the readname
            if '\t' in readname:
                fields = readname.split('\t')
                for field in fields:
                    if field.startswith('BC:Z:'):
                        barcode = field[5:]
                        break
            # Determine which donor the read belongs to
            donor_id = barcode_donor_map.get(barcode, "unmapped") if barcode else "unmapped"
            # Write the read to the appropriate file
            output_files[donor_id].writelines(lines)

    # Close all output files
    for out_file in output_files.values():
        out_file.close()

    print(f"Split {input_fastq} into {len(unique_donors)} temporary donor files, plus unmapped reads")
    return unique_donors, temp_files

def split_samples(args):
    # Parse barcode-to-donor mapping
    barcode_donor_map = parse_barcode_donor_mapping(args.tsv)

    # Handle whether we received a list of files or just one
    input_files = args.input if isinstance(args.input, list) else [args.input]

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Create a temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Created temporary directory: {temp_dir}")

        # Keep track of all donor IDs and their temporary files
        all_donors = set()
        donor_to_temp_files = defaultdict(list)
        filetype = None

        # Process each input file
        for input_file in input_files:
            if input_file.endswith('.bam'):
                filetype = 'bam'
                donors, temp_files = split_bam_by_donor(
                    input_file,
                    barcode_donor_map,
                    temp_dir
                )
            elif is_fastq(input_file):
                filetype = 'fastq'
                donors, temp_files = split_fastq_by_donor(
                    input_file,
                    barcode_donor_map,
                    temp_dir
                )
            else:
                print(f"Unsupported file type: {input_file}")
                continue

            all_donors.update(donors)
            for donor_id, temp_file in temp_files.items():
                donor_to_temp_files[donor_id].append(temp_file)

        # Now concatenate the temporary files for each donor
        for donor_id in all_donors.union({"unmapped"}):
            temp_files_for_donor = donor_to_temp_files[donor_id]
            if filetype == 'bam':
                final_output_path = os.path.join(args.output_dir, f"{donor_id}.bam")
                if len(temp_files_for_donor) == 1:
                    shutil.copy2(temp_files_for_donor[0], final_output_path)
                    print(f"Copied {donor_id} file to {final_output_path}")
                elif len(temp_files_for_donor) > 1:
                    concatenate_bam_files(temp_files_for_donor, final_output_path)
            elif filetype == 'fastq':
                final_output_path = os.path.join(args.output_dir, f"{donor_id}.fastq")
                with open(final_output_path, 'wt') as outfile:
                    for temp_file in temp_files_for_donor:
                        with open(temp_file, 'rt') as infile:
                            shutil.copyfileobj(infile, outfile)
                print(f"Wrote {donor_id} FASTQ to {final_output_path}")

        print("All processing complete. Temporary files will be deleted.")
