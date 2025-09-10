import sys
import time
import logging
from bamurai.core import *
from bamurai.utils import print_elapsed_time_pretty, smart_open, create_progress_bar_for_file, count_reads_async_generic
from bamurai.logging_config import configure_logging

def calculate_split_len(read, target_len: int):
    """Calculate split locations for a read given a target length"""
    if len(read) < target_len:
        return []

    # find the number of splits
    split_loc = round(len(read) / target_len)
    split_size = len(read) // split_loc

    return [i * split_size for i in range(1, split_loc)]

def split_reads(args):
    configure_logging()
    logger = logging.getLogger("bamurai.split")
    logger.info("Running Bamurai split...")
    start_time = time.time()

    total_input_reads = 0
    total_output_reads = 0
    total_unsplit_reads = 0
    # pretty print the arguments in arg: value format
    arg_desc_dict = {
        "reads": "Input file",
        "len_target": "Target length",
        "output": "Output file"
    }
    logger.info("Arguments:")
    for arg, value in vars(args).items():
        if arg in arg_desc_dict:
            logger.info("  %s: %s", arg_desc_dict[arg], value)

    # Read the input reads file
    read_lens = []

    # clear the output file
    if args.output:
        f = smart_open(args.output, "wt", encoding="utf-8")
    else:
        f = None

    # Create progress bar
    pbar = create_progress_bar_for_file(args.reads, "Splitting reads")
    count_thread = count_reads_async_generic(args.reads, pbar)

    for read in parse_reads(args.reads):
        total_input_reads += 1
        pbar.update(1)
        
        split_locs = calculate_split_len(read, target_len = args.len_target)
        split = split_read(read, at = split_locs)

        if len(split) == 1:
            total_unsplit_reads += 1

        for read in split:
            total_output_reads += 1
            read_lens.append(len(read))

            if args.output:
                assert f is not None
                f.write(read.to_fastq())
                f.write("\n")
            else:
                print(read.to_fastq())

    pbar.close()

    if args.output:
        assert f is not None
        f.close()

    avg_read_len = round(sum(read_lens) / len(read_lens)) if read_lens else 0
    logger.info("Total input reads: %d", total_input_reads)
    logger.info("Total output reads: %d", total_output_reads)
    logger.info("Total unsplit reads: %d", total_unsplit_reads)
    logger.info("Average split read length: %d", avg_read_len)
    logger.info("Time taken: %.2f seconds", time.time() - start_time)
