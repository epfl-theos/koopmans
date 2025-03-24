import inspect
import logging
from pathlib import Path


class FormatterWithLineNo(logging.Formatter):
    def format(self, record):
        stack = inspect.stack()
        frame = stack[10]
        file_name = Path(frame.filename).relative_to(Path(__file__).parents[2])
        line_number = frame.lineno

        # Add file name and line number to the log message
        original_message = super().format(record)
        return f'{record.levelname} - {original_message} - {file_name}:{line_number}'


def setup_logging(level=logging.INFO, filename='koopmans.log'):

    file_handler = logging.FileHandler(filename, mode='w')
    file_handler.setFormatter(FormatterWithLineNo())

    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    root_logger.addHandler(file_handler)

    # # Add a file handler for each module dynamically
    # logger = logging.getLogger()
    # if not logger.handlers:
    #     module_name = logger.name if logger.name else 'default'
    #     log_file = f'logs/{module_name}.log'
    #     module_handler = logging.FileHandler(log_file)
    #     module_handler.setFormatter(FormatterWithLineNo())
    #     logger.addHandler(module_handler)
