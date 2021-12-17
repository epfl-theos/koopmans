'''

Miscellaneous functions for koopmans.utils

Written by Edward Linscott May 2020

'''


from typing import List, Union, Generator, Iterable
from collections.abc import Iterable


def calc_diff(calcs, silent=False):
    # Returns the differences in the settings of a list of calculators

    # If calcs is a dict, convert it to a list (we only need the values)
    if isinstance(calcs, dict):
        calcs = calcs.values()

    diffs = []

    settings = [c._settings for c in calcs]

    keys = set([k for s in settings for k in s.keys()])
    for key in sorted(keys):
        vals = [s.get(key, None) for s in settings]
        if len(set(vals)) > 1:
            if not silent:
                print(f'{key}: ' + ', '.join(map(str, vals)))
            diffs.append(key)

    return diffs


def list_to_formatted_str(values: List[int]):
    # Converts a list of integers into the format expected by Wannier90
    # e.g. list_to_formatted_str([1, 2, 3, 4, 5, 7]) = "1-5,7"
    if len(values) == 0:
        raise ValueError('list_to_formatted_str() should not be given an empty list')
    assert all(a > b for a, b in zip(values[1:], values[:-1])), 'values must be monotonically increasing'
    indices: List[Union[int, None]] = [None]
    indices += [i + 1 for i in range(len(values) - 1) if values[i + 1] != values[i] + 1]
    indices += [None]
    sectors = [values[slice(a, b)] for a, b in zip(indices[:-1], indices[1:])]
    out = []
    for sector in sectors:
        if len(sector) == 1:
            out.append(str(sector[0]))
        else:
            out.append(f'{sector[0]}-{sector[-1]}')
    return ','.join(out)


def flatten(l: Union[List, Iterable]) -> Generator:
    # Converts a list of any kind of object (numbers, arrays, lists, strings, ecc.)
    # to a generator
    for item in l:
        if isinstance(item, Iterable) and not isinstance(item, str):
            for x in flatten(item):
                yield x
        else:
            yield item
