#!/usr/bin/env python3

import argparse
import textwrap
from koopmans.io import read_json
from koopmans.workflows.generic import valid_settings

'''
Perform KI/KIPZ calculations
'''


if __name__ == '__main__':

    # Automatically constructing a list of workflow keywords for 'run_koopmans.py --help'
    # from valid_settings
    epilog = ''
    maxlen = max([len(s.name) for s in valid_settings]) + 2
    for s in valid_settings:
        entry = f'  {s.name.ljust(maxlen)}{s.description} ('
        if isinstance(s.type, tuple):
            entry += '/'.join([t.__name__ for t in s.type])
        else:
            entry += s.type.__name__
        entry += f', default {s.default}'
        if s.options is not None and s.type is not bool:
            entry += ', must be ' + '/'.join([str(o) for o in s.options])
        entry += ')'
        for line in textwrap.wrap(entry, subsequent_indent=' ' * (maxlen + 2)):
            epilog += '\n' + line
        epilog += '\n'

    # Construct parser
    parser = argparse.ArgumentParser(
        description='Perform a KI/KIPZ calculation using Quantum Espresso',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="'workflow' block arguments:" + epilog)
    parser.add_argument('json', metavar='system.json', type=str,
                        help='a single JSON file containing the workflow and code settings')

    # Parse arguments
    args = parser.parse_args()

    # Reading in JSON file
    workflow = read_json(args.json)

    # Run workflow
    workflow.run()

