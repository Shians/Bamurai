import pysam
import gzip
from bamurai.utils import create_progress_bar_for_file, count_reads_async_generic

def validate_file(args):
    """Validate a file to ensure it is correctly formatted."""
    file_path = args.reads
    if file_path.endswith('.bam'):
        validate_bam(file_path)
    elif file_path.endswith('.fastq') or file_path.endswith('.fq') or file_path.endswith('.gz'):
        validate_fastq(file_path)
    else:
        print("File must be in BAM or FASTQ format.")

def validate_fastq(file_path):
    """Validate a FASTQ file to ensure it is correctly formatted."""
    if file_path.endswith('.gz'):
        f = gzip.open(file_path, 'rt')
    elif file_path.endswith('.fastq') or file_path.endswith('.fq'):
        f = open(file_path, 'r', encoding='utf-8')
    else:
        raise ValueError("File must be in FASTQ format.")

    # Create progress bar
    pbar = create_progress_bar_for_file(file_path, "Validating FASTQ")
    count_thread = count_reads_async_generic(file_path, pbar)

    record = 0
    while True:
        header = f.readline()
        if not header:
            break  # End of file
        header = header.rstrip()
        seq = f.readline().rstrip()
        plus = f.readline().rstrip()
        qual = f.readline().rstrip()
        record += 1
        pbar.update(1)

        # Check all lines are present
        if not header or not seq or not plus or not qual:
            if not header:
                print(f"Error at record {record}: Missing header line")
            if not seq:
                print(f"Error at record {record}: Missing sequence line")
            if not plus:
                print(f"Error at record {record}: Missing separator line")
            if not qual:
                print(f"Error at record {record}: Missing quality line")
            pbar.close()
            f.close()
            return False

        # Check
        if not str(header).startswith('@'):
            print(f"Error at record {record}: Header does not start with '@'")
            pbar.close()
            f.close()
            return False
        if not str(plus).startswith('+'):
            print(f"Error at record {record}: Separator line does not start with '+'")
            pbar.close()
            f.close()
            return False
        if len(seq) != len(qual):
            print(f"Error at record {record}: Sequence and quality lengths differ")
            pbar.close()
            f.close()
            return False

        # Check that sequence contains only valid IUPAC characters
        valid_chars = 'ACGTURYKMSWBDHVNacgturykmswbhdvn'
        if not all([str(c) in valid_chars for c in seq]):
            print(f"Error at record {record}: Invalid sequence characters")
            all_invalid_pos = [i for i, c in enumerate(seq) if str(c) not in valid_chars]
            all_invalid_char = [seq[i] for i in all_invalid_pos]

            invalid_pos_str = ', '.join([str(i) for i in all_invalid_pos])
            invalid_char_str = ', '.join([str(c) for c in all_invalid_char])

            print(f"Offending character at positions {invalid_pos_str}, characters: {invalid_char_str}")
            pbar.close()
            f.close()
            return False

    pbar.close()
    f.close()

    print(f"{file_path} is a valid FASTQ file with {record} records.")
    return True

def validate_bam(bam_file):
    """Validate a BAM file to ensure it is correctly formatted."""
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        print("Error opening BAM file:", e)
        return False

    # Check header integrity
    if bam.header is None:
        print("Missing or invalid header in BAM file.")
        return False

    # Create progress bar
    pbar = create_progress_bar_for_file(bam_file, "Validating BAM")
    count_thread = count_reads_async_generic(bam_file, pbar)
    record = 0

    try:
        # Iterate through records to ensure they can be read without error
        for read in bam:
            record += 1
            pbar.update(1)
            
            # Minimal check: ensure required fields exist
            if read.query_name is None:
                print(f"Error at record {record}: Missing query name")
                pbar.close()
                return False
            if read.query_sequence is None:
                print(f"Error at record {record}: Missing query sequence")
                pbar.close()
                return False
            if read.query_qualities is None:
                print(f"Error at record {record}: Missing query qualities")
                pbar.close()
                return False
            # check that the sequence and quality lengths are equal
            if len(read.query_sequence) != len(read.query_qualities):
                print(f"Error at record {record}: Sequence and quality lengths differ")
                pbar.close()
                return False

    except Exception as e:
        print(f"Error reading BAM file at record {record}:", e)
        pbar.close()
        return False

    pbar.close()
    print(f"{bam_file} is a valid BAM file with {record} records.")
    return True
