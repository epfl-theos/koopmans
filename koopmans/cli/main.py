#!/usr/bin/env python3

from koopmans.workflows import valid_settings
from koopmans.calculators import qe_bin_directory
from koopmans.io import read, write
from koopmans.utils import chdir
import ase
import os
import argparse
import subprocess
import textwrap
from types import ModuleType
import matplotlib
matplotlib.use('Agg')


'''
Perform KI/KIPZ calculations
'''


def get_version(module):
    if isinstance(module, ModuleType):
        module = module.__path__[0]
    with chdir(module):
        version_label = subprocess.check_output(["git", "describe", "--always", "--tags"]).strip()
    return version_label.decode("utf-8")


def header():

    koopmans_version = get_version(os.path.dirname(__file__))
    ase_version = get_version(ase)
    qe_version = get_version(qe_bin_directory)

    header = [r"  _                                                ",
              r" | | _____   ___  _ __  _ __ ___   __ _ _ __  ___  ",
              r" | |/ / _ \ / _ \| '_ \| '_ ` _ \ / _` | '_ \/ __| ",
              r" |   < (_) | (_) | |_) | | | | | | (_| | | | \__ \ ",
              r" |_|\_\___/ \___/| .__/|_| |_| |_|\__,_|_| |_|___/ ",
              r"                 |_|                               ",
              "",
              " Koopmans spectral functional calculations with Quantum ESPRESSO",
              "",
              " Written by Edward Linscott, Riccardo De Gennaro, and Nicola Colonna",
              "",
              f" using QE version {qe_version}, workflow manager version {koopmans_version}, and ASE version "
              f"{ase_version}"
              ""]
    return '\n'.join(header)


def main():
    # Automatically constructing a list of workflow keywords for 'run_koopmans.py --help'
    # from valid_settings
    epilog = ''
    maxlen = max([len(s.name) for s in valid_settings]) + 2
    for s in valid_settings:
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

    # Write out the header
    print(header())

    # Run workflow
    workflow.run()

    # Save workflow to file
    write(workflow, args.json.replace('.json', '.kwf'))

    # Print farewell message
    print('\n Workflow complete')
