from pathlib import Path
from typing import Any

import dill


class PicklerWithRelativePaths(dill.Pickler):
    def __init__(self, base_directory, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.base_directory = base_directory

    def reducer_override(self, obj):
        if isinstance(obj, Path):
            return Path, pickle_path(obj, self.base_directory)
        else:
            return NotImplemented


def pickle_path(path: Path, base_directory: Path):
    # Convert the path to relative if it's absolute
    if path.is_absolute():
        try:
            rel_path = path.relative_to(base_directory.resolve())
        except ValueError:
            rel_path = path
        return (str(rel_path),)
    else:
        # Return relative paths as-is
        return (str(path),)


def read_pkl(filename: Path | str) -> Any:
    with open(filename, 'rb') as fd:
        out = dill.load(fd)
    return out


def write_pkl(obj: Any, filename: Path | str, base_directory: Path | None = None):
    filename = Path(filename) if not isinstance(filename, Path) else filename

    if base_directory is None:
        base_directory = filename.parent

    with open(filename, 'wb') as fd:
        pickler = PicklerWithRelativePaths(base_directory, fd)
        pickler.dump(obj)
