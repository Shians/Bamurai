import pysam
import csv
import tempfile
import shutil
import os

def assign_samples(args):
    """
    Assign donor_id to RG tag in BAM file using barcode-to-donor mapping TSV.
    Args:
        args: argparse.Namespace with .bam (list), .tsv (str), .output (str)
    """
    # Load barcode-to-donor mapping
    barcode_to_donor = {}
    with open(args.tsv, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            barcode = row['barcode']
            donor_id = row['donor_id']
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
                for read in infile:
                    # Extract barcode (assume in CB tag or RX tag)
                    barcode = None
                    if read.has_tag('CB'):
                        barcode = read.get_tag('CB')
                    elif read.has_tag('RX'):
                        barcode = read.get_tag('RX')
                    if barcode and barcode in barcode_to_donor:
                        donor_id = barcode_to_donor[barcode]
                        read.set_tag('RG', donor_id, value_type='Z')
                    else:
                        no_match_count += 1
                    outfile.write(read)
        shutil.move(tmp_output, args.output)
        print(f"RG tags assigned and written to {args.output}")
        if no_match_count > 0:
            print(f"{no_match_count} reads did not have a matching barcode in the mapping and were left unchanged.")
    finally:
        # Clean up temp file if something went wrong
        if os.path.exists(tmp_output):
            os.remove(tmp_output)
