import gzip
import time
import logging
from bamurai.logging_config import LOGGING_FORMAT, LOGGING_DATEFMT

def print_elapsed_time_pretty(start_time, logger=None):
    """Log elapsed time in a pretty format."""
    elapsed_time = time.time() - start_time

    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    if minutes < 1:
        seconds = round(elapsed_time, 2)
    else:
        seconds = int(elapsed_time % 60)

    if logger is None:
        logger = logging.getLogger("bamurai.utils")

    # below 5 minutes log in seconds
    if elapsed_time < 300:
        logger.info("Time taken: %ss", seconds)
    # below 1 hour log in minutes and seconds
    elif elapsed_time < 3600:
        logger.info("Time taken: %dm %ds", minutes, seconds)
    # above 1 hour log in hours, minutes and seconds
    else:
        logger.info("Time elapsed: %dh %dm %ds", hours, minutes, seconds)

def is_fastq(path):
    """Check if a file is a FASTQ file."""
    path = path.lower()
    return path.endswith(".fastq") or \
        path.endswith(".fq") or \
        path.endswith(".fastq.gz") or \
        path.endswith(".fq.gz")

def smart_open(filename, mode="rt", encoding=None):
    """Open a file normally or with gzip based on file extension. Supports text and binary modes."""
    if filename.endswith('.gz'):
        if 't' in mode and encoding is not None:
            return gzip.open(filename, mode, encoding=encoding)
        else:
            return gzip.open(filename, mode)
    else:
        if 't' in mode and encoding is not None:
            return open(filename, mode, encoding=encoding)
        else:
            return open(filename, mode)

def calculate_percentage(count, total):
    """Calculate percentage with safety check for division by zero."""
    return (count / total * 100) if total > 0 else 0

# Progress bar utilities
def create_progress_bar_for_file(filepath, desc="Processing", unit="reads", mininterval=0.2):
    """Create appropriate progress bar based on file type using tqdm's built-in time buffering."""
    import os
    from tqdm import tqdm

    if filepath.endswith(('.bam', '.sam', '.cram')):
        # BAM files: Start with unknown total, will be updated async
        return tqdm(desc=desc, unit=unit, mininterval=mininterval)
    elif filepath.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        # FASTQ files: Use file size for progress
        file_size = os.path.getsize(filepath)
        return tqdm(total=file_size, desc=desc, unit="B", unit_scale=True, mininterval=mininterval)
    else:
        # Unknown format: Use unknown total
        return tqdm(desc=desc, unit=unit, mininterval=mininterval)

def count_reads_async_generic(filepath, progress_bar):
    """Count reads asynchronously for different file types."""
    import threading
    import pysam

    def count_bam_reads():
        try:
            with pysam.AlignmentFile(filepath, "rb") as infile:
                total = sum(1 for read in infile if not (read.is_secondary or read.is_supplementary))
                progress_bar.total = total
                progress_bar.refresh()
        except Exception:
            pass

    def count_fastq_reads():
        try:
            line_count = 0
            with smart_open(filepath, "r") as f:
                for _ in f:
                    line_count += 1
            total_reads = line_count // 4
            progress_bar.total = total_reads
            progress_bar.refresh()
        except Exception:
            pass

    if filepath.endswith(('.bam', '.sam', '.cram')):
        count_func = count_bam_reads
    elif filepath.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        count_func = count_fastq_reads
    else:
        return  # Unknown format

    count_thread = threading.Thread(target=count_func)
    count_thread.daemon = True
    count_thread.start()
    return count_thread

def create_multi_file_progress_bar(filepaths, desc="Processing files", mininterval=0.2):
    """Create progress bar for multiple files using tqdm's built-in time buffering."""
    from tqdm import tqdm
    import threading
    import pysam

    total_pbar = tqdm(desc=desc, unit="reads", mininterval=mininterval)

    def count_all_reads():
        try:
            total = 0
            for filepath in filepaths:
                if filepath.endswith(('.bam', '.sam', '.cram')):
                    with pysam.AlignmentFile(filepath, "rb") as infile:
                        total += sum(1 for read in infile if not (read.is_secondary or read.is_supplementary))
                elif filepath.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
                    with smart_open(filepath, "r") as f:
                        line_count = sum(1 for _ in f)
                    total += line_count // 4
            total_pbar.total = total
            total_pbar.refresh()
        except Exception:
            pass

    count_thread = threading.Thread(target=count_all_reads)
    count_thread.daemon = True
    count_thread.start()

    return total_pbar, count_thread
