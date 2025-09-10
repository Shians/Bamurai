import pysam
import csv
import tempfile
import shutil
import os
import threading
from tqdm import tqdm
from bamurai.utils import calculate_percentage

def count_reads_async(bam_path, progress_bar):
    """Count total reads in BAM file asynchronously and update progress bar total."""
    try:
        with pysam.AlignmentFile(bam_path, "rb") as infile:
            total = sum(1 for _ in infile)
            progress_bar.total = total
            progress_bar.refresh()
    except Exception:
        # If counting fails, leave progress bar without total
        pass

def assign_samples(args):
    """
    Assign donor_id to RG tag in BAM file using barcode-to-donor mapping TSV.
    Args:
        args: argparse.Namespace with .bam (list), .tsv (str), .output (str), .barcode_column (str, optional), .donor_id_column (str, optional)
    """
    # Load barcode-to-donor mapping with flexible column names
    with open(args.tsv, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        columns = reader.fieldnames
        if not columns:
            raise ValueError("The TSV file is empty or improperly formatted.")
        # Determine barcode column
        barcode_col = getattr(args, 'barcode_column', None)
        donor_id_col = getattr(args, 'donor_id_column', None)
        # Barcode column detection logic
        if not barcode_col:
            if 'cell' in columns and 'barcode' in columns:
                raise ValueError("Both 'cell' and 'barcode' columns are present in the TSV. Please specify --barcode-column explicitly.")
            elif 'barcode' in columns:
                barcode_col = 'barcode'
            elif 'cell' in columns:
                barcode_col = 'cell'
            else:
                raise ValueError("No 'barcode' or 'cell' column found in the TSV. Please specify --barcode-column.")
        if not donor_id_col:
            if 'donor_id' in columns:
                donor_id_col = 'donor_id'
            else:
                raise ValueError("No 'donor_id' column found in the TSV. Please specify --donor-id-column.")
        barcode_to_donor = {}
        for row in reader:
            barcode = row[barcode_col]
            donor_id = row[donor_id_col]
            barcode_to_donor[barcode] = donor_id

    # Open output BAM for writing to a temporary file
    with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as tmpfile:
        tmp_output = tmpfile.name
    try:
        with pysam.AlignmentFile(args.bam, "rb") as infile:
            header = infile.header.to_dict()
            # Add RG records for each donor_id
            donor_ids = set(barcode_to_donor.values())
            header.setdefault('RG', [])
            for donor_id in donor_ids:
                header['RG'].append({'ID': donor_id, 'SM': donor_id})
            with pysam.AlignmentFile(tmp_output, "wb", header=header) as outfile:
                no_match_count = 0
                total_reads = 0
                donor_counts = {donor_id: 0 for donor_id in donor_ids}

                # Initialize progress bar with unknown total
                pbar = tqdm(desc="Processing reads", unit="reads")

                # Start background thread to count total reads
                count_thread = threading.Thread(target=count_reads_async, args=(args.bam, pbar))
                count_thread.daemon = True
                count_thread.start()

                for read in infile:
                    total_reads += 1
                    pbar.update(1)

                    # Extract barcode (assume in CB tag or RX tag)
                    barcode = None

                    if read.has_tag('CB'):
                        barcode = read.get_tag('CB')
                    elif read.has_tag('RX'):
                        barcode = read.get_tag('RX')

                    if barcode and (barcode in barcode_to_donor):
                        donor_id = barcode_to_donor[barcode]
                        read.set_tag('RG', donor_id, value_type='Z')
                        donor_counts[donor_id] += 1
                    else:
                        no_match_count += 1

                    outfile.write(read)

                pbar.close()
        shutil.move(tmp_output, args.output)
        print(f"RG tags assigned and written to {args.output}")
        print(f"Total reads processed: {total_reads}")
        print(f"Reads assigned to donors:")
        for donor_id in sorted(donor_counts.keys()):
            percentage = calculate_percentage(donor_counts[donor_id], total_reads)
            print(f"  {donor_id}: {donor_counts[donor_id]} [{percentage:.1f}%]")
        unassigned_percentage = calculate_percentage(no_match_count, total_reads)
        print(f"Unassigned reads: {no_match_count} [{unassigned_percentage:.1f}%]")
    finally:
        # Clean up temp file if something went wrong
        if os.path.exists(tmp_output):
            os.remove(tmp_output)
