#!/usr/bin/env python3

import argparse
import textwrap
from koopmans_cp import workflow

'''
Perform KI/KIPZ calculations
'''


if __name__ == '__main__':

    # Automatically constructing a list of workflow keywords for 'koopmans.py --help'
    # from workflow.valid_settings
    epilog = ''
    maxlen = max([len(s.name) for s in workflow.valid_settings]) + 2
    for s in workflow.valid_settings:
        entry = f'  {s.name.ljust(maxlen)}{s.description} ({s.type.__name__}, ' \
            f'default {s.default}'
        if s.options is not None and s.type is not bool:
            entry += ', must be ' + '/'.join([str(o) for o in s.options])
        entry += ')'
        for line in textwrap.wrap(entry, subsequent_indent=' '*(maxlen + 2)):
            epilog += '\n' + line
        epilog += '\n'

    # Construct parser
    parser = argparse.ArgumentParser(
        description='Perform a KI/KIPZ calculation using cp.x',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="'calc_param' block arguments:" + epilog)
    parser.add_argument('json', metavar='system.json', type=str,
                        help='a single JSON file containing the workflow and cp.x settings')

    # Parse arguments
    args = parser.parse_args()

    # Run workflow
    workflow.run_from_json(args.json)
