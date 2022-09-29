'koopmans: a package for performing and automating Koopmans functional calculations'
import sys
from pathlib import Path

if sys.version_info >= (3, 8):
    from importlib import metadata
else:
    import importlib_metadata as metadata

__version__ = metadata.version('koopmans')
