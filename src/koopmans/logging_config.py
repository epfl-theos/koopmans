"""Logging configuration for `koopmans`."""

import logging
from pathlib import Path


class FormatterWithLineNo(logging.Formatter):
    """Logging formatter that includes the file name and line number."""

    def format(self, record) -> str:
        """Format the log message with file name and line number."""
        original_message = super().format(record)
        return f'{record.levelname} - {original_message} - {record.pathname}:{record.lineno}'


def setup_logging(level=logging.INFO, filename='koopmans.log'):
    """Set up logging."""
    # Remove pre-existing log files
    if Path(filename).exists():
        Path(filename).unlink()

    file_handler = logging.FileHandler(filename, mode='a')
    file_handler.setFormatter(FormatterWithLineNo())

    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # Clear existing handlers (i.e. those added by AiiDA)
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    root_logger.addHandler(file_handler)
