"""warnings functions for koopmans.utils."""

import sys
import traceback
import warnings
from typing import Optional, TextIO, Type, Union

from ._io import print_alert


class CalculatorNotConvergedWarning(UserWarning):
    """Warning for when a calculator has not converged."""

    pass


def _warning(message: Union[str, Warning], category: Type[Warning] = UserWarning, filename: str = '',
             lineno: int = -1, file: Optional[TextIO] = None, line: Optional[str] = None) -> None:
    """Print a patched warnings.warn with added aesthetics."""
    print_alert('warning', str(message))


def _warn_with_traceback(message: Union[str, Warning], category: Type[Warning] = UserWarning, filename: str = '',
                         lineno: int = -1, file: Optional[TextIO] = None, line: Optional[str] = None) -> None:
    """Print a warning with traceback."""
    log = file if file and hasattr(file, 'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))


# warnings.showwarning = _warn_with_traceback
warnings.showwarning = _warning


def warn(message: str, cls: Type[Warning] = UserWarning) -> None:
    """Print a warning.

    This is implemented here to allow the monkey-patched warnings.warn to be imported as utils.warn.
    """
    warnings.warn(message, cls)
