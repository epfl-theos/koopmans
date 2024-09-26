import os
from pathlib import Path
from typing import Any

import dill


def read_pkl(filename: Path | str) -> Any:
    with open(filename, 'rb') as fd:
        out = dill.load(fd)
    return out


def write_pkl(obj: Any, filename: Path | str):
    filename = Path(filename) if not isinstance(filename, Path) else filename

    with open(filename, 'wb') as fd:
        dill.dump(obj, fd)
