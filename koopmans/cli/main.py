#!/usr/bin/env python3

import argparse
import textwrap

import koopmans.mpl_config
from koopmans.io import read
from koopmans.settings import WorkflowSettingsDict


def main():
    # Construct parser
    parser = argparse.ArgumentParser(
        description='Perform Koopmans functional calculations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='See https://koopmans-functionals.org for more details')
    parser.add_argument('json', metavar='system.json', type=str,
                        help='a single JSON file containing the workflow and code settings')

    # Parse arguments
    args = parser.parse_args()

    # Reading in JSON file
    workflow = read(args.json)

    # Run workflow
    workflow.run()
