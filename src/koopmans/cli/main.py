#!/usr/bin/env python3

import argparse
import re
import sys
import textwrap
import traceback

import koopmans.mpl_config
from koopmans.io import read
from koopmans.utils import indented_print


def _custom_exception_hook(exception_type, exception_value, traceback):
    # Adding spaces to the error name so that it is easier to read
    spaced_text = re.sub(r'([a-z])([A-Z])', r'\1 \2', exception_type.__name__)
    spaced_text = re.sub(r'([A-Z]+)([A-Z][a-z])', r'\1 \2', spaced_text)

    indented_print(f'\n\n  > [!CAUTION] {spaced_text}')
    indented_print(str(exception_value) + '\n', initial_indent='  > ', subsequent_indent='  > ')


def main():
    # Construct parser
    parser = argparse.ArgumentParser(
        description='Perform Koopmans functional calculations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='See https://koopmans-functionals.org for more details')
    parser.add_argument('json', metavar='system.json', type=str,
                        help='a single JSON file containing the workflow and code settings')
    parser.add_argument('-t', '--traceback', action='store_true', help='enable traceback')

    # Parse arguments
    args = parser.parse_args()

    # Reading in JSON file
    workflow = read(args.json)

    # Set traceback behavior
    if not args.traceback:
        sys.tracebacklimit = 0
        sys.excepthook = _custom_exception_hook

    # Run workflow
    workflow.run()
