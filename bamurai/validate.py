import pysam
import gzip

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
            return False

        # Check
        if not header.startswith('@'):
            print(f"Error at record {record}: Header does not start with '@'")
            return False
        if not plus.startswith('+'):
            print(f"Error at record {record}: Separator line does not start with '+'")
            return False
        if len(seq) != len(qual):
            print(f"Error at record {record}: Sequence and quality lengths differ")
            return False

        # Check that sequence contains only valid IUPAC characters
        valid_chars = 'ACGTURYKMSWBDHVNacgturykmswbhdvn'
        if not all([c in valid_chars for c in seq]):
            print(f"Error at record {record}: Invalid sequence characters")
            all_invalid_pos = [i for i, c in enumerate(seq) if c not in valid_chars]
            all_invalid_char = [seq[i] for i in all_invalid_pos]

            invalid_pos_str = ', '.join([str(i) for i in all_invalid_pos])
            invalid_char_str = ', '.join(all_invalid_char)

            print(f"Offending character at positions {invalid_pos_str}, characters: {invalid_char_str}")
            return False


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

    try:
        # Iterate through records to ensure they can be read without error
        record = 0
        for read in bam:
            record += 1
            # Minimal check: ensure required fields exist
            if read.query_name is None:
                print(f"Error at record {record}: Missing query name")
                return False
            if read.query_sequence is None:
                print(f"Error at record {record}: Missing query sequence")
                return False
            if read.query_qualities is None:
                print(f"Error at record {record}: Missing query qualities")
                return False
            # check that the sequence and quality lengths are equal
            if len(read.query_sequence) != len(read.query_qualities):
                print(f"Error at record {record}: Sequence and quality lengths differ")
                return False

    except Exception as e:
        print(f"Error reading BAM file at record {record}:", e)
        return False

    print(f"{bam_file} is a valid BAM file with {record} records.")
    return True
