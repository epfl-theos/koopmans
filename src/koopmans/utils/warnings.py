"""warnings functions for koopmans.utils."""

import warnings
from typing import Optional, Type, Union

from ._io import create_alert


class IncommensurateProjectionsWarning(UserWarning):
    """Warning for when the projections are incommensurate with the number of occupied bands."""


class CalculatorNotConvergedWarning(UserWarning):
    """Warning for when a calculator has not converged."""


def custom_formatwarning(message: Union[str, Warning], category: Type[Warning] = UserWarning,
                         filename: str = '', lineno: int = -1, line: Optional[str] = None) -> str:
    """Print a patched warnings.warn with added aesthetics."""
    return create_alert("warning", str(message))


def warn(message: str, cls: Type[Warning] = UserWarning) -> None:
    """Print a warning.

    This is implemented here to allow the monkey-patched warnings.warn to be imported as utils.warn.
    """
    warnings.warn(message, cls, stacklevel=2)
    # logger = logging.getLogger(__name__)
    # logger.warning(f'{cls.__name__}: {message}', stacklevel=2)


def configure_warnings() -> None:
    """Configure thw warnings module."""
    warnings.filterwarnings("once", category=IncommensurateProjectionsWarning)
    warnings.formatwarning = custom_formatwarning
