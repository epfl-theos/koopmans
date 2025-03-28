"""Miscellaneous functions for `koopmans.utils`."""


from typing import Any, Generator, Iterable, List, Union


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
