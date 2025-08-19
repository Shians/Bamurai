import logging

LOGGING_FORMAT = "%(levelname)s [%(asctime)s] %(message)s"
LOGGING_DATEFMT = "%Y-%m-%d %H:%M:%S"

def configure_logging(level=logging.INFO):
    logging.basicConfig(
        format=LOGGING_FORMAT,
        datefmt=LOGGING_DATEFMT,
        level=level
    )
