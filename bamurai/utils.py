import time
import sys

def print_elapsed_time_pretty(start_time):
    """Print elapsed time in a pretty format."""
    elapsed_time = time.time() - start_time

    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    if minutes < 1:
        seconds = round(elapsed_time, 2)
    else:
        seconds = int(elapsed_time % 60)

    # below 5 minutes print in seconds
    if elapsed_time < 300:
        print(f"Time taken: {seconds}s", file=sys.stderr)
    # below 1 hour print in minutes and seconds
    elif elapsed_time < 3600:
        print(f"Time taken: {minutes}m {seconds}s", file=sys.stderr)
    # above 1 hour print in hours, minutes and seconds
    else:
        print(f"Time elapsed: {hours}h {minutes}m {seconds}s", file=sys.stderr)

def is_fastq(path):
    """Check if a file is a FASTQ file."""
    path = path.lower()
    return path.endswith(".fastq") or \
        path.endswith(".fq") or \
        path.endswith(".fastq.gz") or \
        path.endswith(".fq.gz")
