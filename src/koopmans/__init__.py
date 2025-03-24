'koopmans: a package for performing and automating Koopmans functional calculations'
import sys
from pathlib import Path

from .logging_config import setup_logging

if sys.version_info >= (3, 8):
    from importlib import metadata
else:
    import importlib_metadata as metadata

__version__: str = metadata.version('koopmans')
