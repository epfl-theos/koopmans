#!/usr/bin/env python3

import argparse
import textwrap
import koopmans.mpl_config
from koopmans.settings import WorkflowSettingsDict
from koopmans.io import read, write


'''
Perform KI/KIPZ calculations
'''


def main():
    # Automatically constructing a list of workflow keywords for 'run_koopmans.py --help'
    # from valid_settings
    epilog = ''
    wf_settings = WorkflowSettingsDict()
    maxlen = max([len(s.name) for s in wf_settings.settings]) + 2
    for s in wf_settings.settings:
        entry = f'  {s.name.ljust(maxlen)}{s.description} ('
        if isinstance(s.kind, tuple):
            entry += '/'.join([t.__name__ for t in s.kind])
        else:
            entry += s.kind.__name__
        if s.default is not None:
            entry += f', default {s.default}'
        if s.options is not None and s.kind is not bool:
            entry += ', must be ' + '/'.join([str(o) for o in s.options])
        entry += ')'
        for line in textwrap.wrap(entry, subsequent_indent=' ' * (maxlen + 2), width=120):
            epilog += '\n' + line
        epilog += '\n'

    # Construct parser
    parser = argparse.ArgumentParser(
        description='Perform a Koopmans spectral functional calculation using Quantum ESPRESSO',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="'workflow' block arguments:" + epilog)
    parser.add_argument('json', metavar='system.json', type=str,
                        help='a single JSON file containing the workflow and code settings')

    # Parse arguments
    args = parser.parse_args()

    # Reading in JSON file
    workflow = read(args.json)

    # Run workflow
    workflow.run()
