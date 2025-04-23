"""Logging configuration for `koopmans`."""

import inspect
import logging
from pathlib import Path


class FormatterWithLineNo(logging.Formatter):
    """Logging formatter with line number information."""

    def format(self, record):  # noqa: A003
        """Format the log message, adding the file name and line number where the log was called."""
        stack = inspect.stack()
        frame = stack[10]
        file_name = Path(frame.filename).relative_to(Path(__file__).parents[2])
        line_number = frame.lineno

        # Add file name and line number to the log message
        original_message = super().format(record)
        return f'{record.levelname} - {original_message} - {file_name}:{line_number}'


def setup_logging(level=logging.INFO, filename='koopmans.log'):
    """Set up logging."""
    file_handler = logging.FileHandler(filename, mode='w')
    file_handler.setFormatter(FormatterWithLineNo())

    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    root_logger.addHandler(file_handler)
