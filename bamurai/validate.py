import pysam
import gzip

def validate_fastq(file_path):
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as f:
            record = 0
            while True:
                header = f.readline().rstrip()
                if not header:
                    break
    elif file_path.endswith('.fastq') or file_path.endswith('.fq'):
        with open(file_path, 'r', encoding='utf-8') as f:
            record = 0
            while True:
                header = f.readline().rstrip()
                if not header:
                    break  # End of file
                seq = f.readline().rstrip()
                plus = f.readline().rstrip()
                qual = f.readline().rstrip()
                record += 1

                if not header.startswith('@'):
                    print(f"Error at record {record}: Header does not start with '@'")
                    return False
                if not plus.startswith('+'):
                    print(f"Error at record {record}: Separator line does not start with '+'")
                    return False
                if len(seq) != len(qual):
                    print(f"Error at record {record}: Sequence and quality lengths differ")
                    return False
    else:
        print("Error: File must be in FASTQ format.")

    print("FASTQ file is valid.")
    return True
