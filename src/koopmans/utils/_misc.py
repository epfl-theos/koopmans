"""Miscellaneous functions for `koopmans.utils`."""

from typing import Any, Generator, Iterable, List, Union

import numpy as np


def flatten(lst: Union[List[Any], Iterable[Any]]) -> Generator[Any, None, None]:
    """Convert a list of any kind of object (numbers, arrays, lists, strings, ecc.) to a generator."""
    for item in lst:
        if isinstance(item, Iterable) and not isinstance(item, str):
            for x in flatten(item):
                yield x
        else:
            yield item


def update_nested_dict(dct_to_update, second_dct):
    """Recursively update a dictionary with another dictionary."""
    for k, v in second_dct.items():
        if k in dct_to_update and isinstance(v, dict):
            update_nested_dict(dct_to_update[k], second_dct[k])
        else:
            dct_to_update[k] = v


def remove_numpy_from_dict(obj: Any) -> Any:
    """
    Recursively convert NumPy arrays to Python lists inside nested structures.
    Supports dict, list, tuple, set, and nested combinations.
    """
    if isinstance(obj, dict):
        return {k: remove_numpy_from_dict(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [remove_numpy_from_dict(v) for v in obj]
    elif isinstance(obj, tuple):
        return tuple(remove_numpy_from_dict(v) for v in obj)
    elif isinstance(obj, set):
        return {remove_numpy_from_dict(v) for v in obj}
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj
