import gzip
from itertools import zip_longest

def get_hto(args):
    """
    Extracts HTO information from 10x FASTQ files.
    """
    r1_file = args.r1
    r2_file = args.r2
    barcode_len = args.bc_len
    umi_len = args.umi_len
    output_file = args.output
    hashtag_len = args.hashtag_len
    hashtag_left_buffer = args.hashtag_left_buffer

    r1_open = gzip.open if r1_file.endswith('.gz') else open
    r2_open = gzip.open if r2_file.endswith('.gz') else open

    with r1_open(r1_file, 'rt') as f_r1, r2_open(r2_file, 'rt') as f_r2:
        with open(output_file, 'w', encoding='utf-8') as out_f:
            out_f.write("read_name\tcell_barcode\tumi\thto\n")
            while True:
                r1_lines = [f_r1.readline() for _ in range(4)]
                r2_lines = [f_r2.readline() for _ in range(4)]
                if not r1_lines[0] or not r2_lines[0]:
                    break

                read_name = r1_lines[0].strip().split()[0][1:]  # Remove '@' and take first word
                r1_seq = r1_lines[1].strip()
                r2_seq = r2_lines[1].strip()

                cell_barcode = r1_seq[:barcode_len]
                umi = r1_seq[barcode_len:(barcode_len + umi_len)]
                hto = r2_seq[hashtag_left_buffer:(hashtag_left_buffer + hashtag_len)]

                out_f.write(f"{read_name}\t{cell_barcode}\t{umi}\t{hto}\n")
