'''

Miscellaneous functions for koopmans.utils

Written by Edward Linscott May 2020

'''


from typing import List, Union, Generator, Iterable


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


def flatten(l: Union[List, Iterable]) -> Generator:
    # Converts a list of any kind of object (numbers, arrays, lists, strings, ecc.)
    # to a generator
    for item in l:
        if isinstance(item, Iterable) and not isinstance(item, str):
            for x in flatten(item):
                yield x
        else:
            yield item
