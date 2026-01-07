"""warnings functions for koopmans.utils."""

import sys
import warnings
from typing import Optional, Set, Tuple, Type, Union

from ._io import create_alert


class IncommensurateProjectionsWarning(UserWarning):
    """Warning for when the projections are incommensurate with the number of occupied bands."""


class ImageCorrectionWarning(UserWarning):
    """Warnings for when image charge corrections are/aren't used when appropriate."""


class CalculatorNotConvergedWarning(UserWarning):
    """Warning for when a calculator has not converged."""


# Global set to track which warnings have already been shown
# This prevents duplicates when Pydantic revalidates model instances
_shown_warnings: Set[Tuple[str, str]] = set()


def custom_formatwarning(message: Union[str, Warning], category: Type[Warning] = UserWarning,
                         filename: str = '', lineno: int = -1, line: Optional[str] = None) -> str:
    """Print a patched warnings.warn with added aesthetics."""
    return create_alert("warning", str(message))


def warn(message: str, cls: Type[Warning] = UserWarning) -> None:
    """Print a warning.

    This is implemented here to allow the monkey-patched warnings.warn to be imported as utils.warn.
    """
    # For warnings that should only appear once, check our global registry first
    # For some reason warnings.simplefilter("once") does not always work as expected with Pydantic (hence
    # this implementation)
    if cls in (ImageCorrectionWarning, IncommensurateProjectionsWarning):
        key = (cls.__name__, message)
        if key in _shown_warnings:
            return  # Already shown, skip it
        _shown_warnings.add(key)

    warnings.warn(message, cls, stacklevel=2)
    # import logging
    # logger = logging.getLogger(__name__)
    # logger.warning(f'{cls.__name__}: {message}', stacklevel=2)


def configure_warnings() -> None:
    """Configure the warnings module."""
    warnings.formatwarning = custom_formatwarning
