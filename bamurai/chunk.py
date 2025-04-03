from bamurai.core import parse_reads

def chunk_reads(args):
    # Convert size string to bytes
    units = {'K': 1024, 'M': 1024**2, 'G': 1024**3}
    size = args.size.upper()
    unit = size[-1]
    if unit in units:
        chunk_size = int(size[:-1]) * units[unit]
    else:
        chunk_size = int(size)

    # Do the chunking
    do_chunk(args.reads, chunk_size, args.prefix)

def do_chunk(input_file, chunk_size, output_prefix):
    """Split input file into chunks of at least chunk_size bytes"""
    current_size = 0
    current_chunk = 1
    current_out = None

    for read in parse_reads(input_file):
        # Open new file if needed
        if current_out is None:
            current_out = open(f"{output_prefix}_{current_chunk}.fastq", 'w')
            current_size = 0

        # Write read
        read_str = read.to_fastq() + "\n"
        current_out.write(read_str)
        current_size += len(read_str.encode('utf-8'))

        # Check if chunk is big enough
        if current_size >= chunk_size:
            current_out.close()
            current_out = None
            current_chunk += 1

    if current_out:
        current_out.close()
