import pysam
import tempfile
import shutil
import os
from bamurai.utils import calculate_percentage, create_progress_bar_for_file, count_reads_async_generic
from bamurai.utils_samples import get_read_barcode, parse_barcode_donor_mapping

def assign_samples(args):
    """
    Assign donor_id to RG tag in BAM file using barcode-to-donor mapping TSV.
    Args:
        args: argparse.Namespace with .bam (str), .tsv (str), .output (str), .barcode_column (str, optional), .donor_id_column (str, optional)
    """
    # Load barcode-to-donor mapping with flexible column names
    barcode_to_donor = parse_barcode_donor_mapping(
        args.tsv,
        getattr(args, 'barcode_column', None),
        getattr(args, 'donor_id_column', None),
    )

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

                # Initialize progress bar using generic utility
                pbar = create_progress_bar_for_file(args.bam, "Processing reads")

                # Start background thread to count total reads
                count_thread = count_reads_async_generic(args.bam, pbar)

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
