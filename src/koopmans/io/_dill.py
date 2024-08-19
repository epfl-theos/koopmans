from pathlib import Path
from typing import Any, TextIO

import dill


def read_pkl(filename: Path | str) -> Any:
    with open(filename, 'rb') as fd:
        out = dill.load(fd)
    return out


def write_pkl(obj: Any, filename: Path | str):
    with open(filename, 'wb') as fd:
        dill.dump(obj, fd)
