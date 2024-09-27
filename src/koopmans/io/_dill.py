import os
from pathlib import Path
from typing import Any, Optional

import dill

BASE_PLACEHOLDER = '/BASE_DIRECTORY_PLACEHOLDER/'


class CustomPicklerUsePlaceholder(dill.Pickler):
    def __init__(self, base_directory, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.base_directory = base_directory

    def persistent_id(self, obj):
        """Called when pickling objects. If obj is a Path, return a unique identifier."""
        if isinstance(obj, Path) and obj.is_absolute():
            if self.base_directory == obj:
                return ('Path', BASE_PLACEHOLDER)
            # Return a tuple to identify this as a Path that needs special handling
            tup = ('Path', str((Path(BASE_PLACEHOLDER) / os.path.relpath(obj, self.base_directory))))
            if 'data' in tup[1]:
                assert '..' in tup[1]
            return tup
        return None  # No special handling required for other objects


class CustomUnpicklerReplacePlaceholder(dill.Unpickler):
    def __init__(self, base_directory, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.base_directory = base_directory.resolve()

    def persistent_load(self, pid):
        """Called when unpickling to resolve persistent ids."""
        if isinstance(pid, tuple) and pid[0] == 'Path' and pid[1].startswith(BASE_PLACEHOLDER):
            # Replace the placeholder with the actual current working directory
            path = (self.base_directory / pid[1][len(BASE_PLACEHOLDER):]).resolve()
            if 'data' in pid[1]:
                assert 'tmp' not in str(path)
            return path
        return pid  # If it's not a special persistent id, return as is


class CustomUnpicklerKeepPlaceholder(dill.Unpickler):
    def persistent_load(self, pid):
        if isinstance(pid, tuple) and pid[0] == 'Path':
            return Path(pid[1])


def read_pkl(filename: Path | str, base_directory: Optional[Path] = None) -> Any:
    with open(filename, 'rb') as fd:
        if base_directory is None:
            unpickler = CustomUnpicklerKeepPlaceholder(fd)
        else:
            unpickler = CustomUnpicklerReplacePlaceholder(base_directory, fd)
        out = unpickler.load()

    return out


def write_pkl(obj: Any, filename: Path | str, base_directory: Path | None = None):
    filename = Path(filename) if not isinstance(filename, Path) else filename
    with open(filename, 'wb') as fd:
        if base_directory is None:
            pickler = dill.Pickler(fd)
        else:
            pickler = CustomPicklerUsePlaceholder(base_directory, fd)
        pickler.dump(obj)
