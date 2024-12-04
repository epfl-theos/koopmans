#!/usr/bin/env python3

import argparse
import re
import sys
import textwrap
import traceback

import koopmans.mpl_config
from koopmans.engines import LocalhostEngine
from koopmans.io import read, write
from koopmans.status import Status
from koopmans.utils import print_alert


def _custom_exception_hook(exception_type, exception_value, traceback):
    # Adding spaces to the error name so that it is easier to read
    spaced_text = re.sub(r'([a-z])([A-Z])', r'\1 \2', exception_type.__name__)
    spaced_text = re.sub(r'([A-Z]+)([A-Z][a-z])', r'\1 \2', spaced_text)

    print_alert('caution', str(exception_value), header=spaced_text, indent=1)


def main():
    from koopmans import __version__

    # Construct parser
    parser = argparse.ArgumentParser(
        description='Perform Koopmans functional calculations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='See https://koopmans-functionals.org for more details')
    parser.add_argument('json', metavar='system.json', type=str,
                        help='a single JSON file containing the workflow and code settings')
    parser.add_argument('-t', '--traceback', action='store_true', help='enable traceback')
    parser.add_argument('--engine', choices=['localhost', 'aiida'], default='localhost',
                        help="Specify the execution engine: 'local' or 'aiida' (default: 'local')")
    parser.add_argument('--engine_config', type=str, default='engine.json',
                        help='Specify the engine configuration file (default: engine.json)')
    parser.add_argument('--version', action='version', version=__version__,
                        help="Show the program's version number and exit")

    # Parse arguments
    args = parser.parse_args()

    # Reading in JSON file
    workflow = read(args.json)

    # Create the engine
    if args.engine == 'localhost':
        engine = LocalhostEngine(from_scratch=workflow.parameters.from_scratch)
    elif args.engine == 'aiida':
        from aiida_koopmans.engine.aiida import AiiDAEngine
        if args.engine_config is not None:
            import json
            with open(args.engine_config, 'r') as f:
                engine_config = json.load(f)
        else:
            engine_config = None
        engine = AiiDAEngine(configuration=engine_config)
    else:
        raise NotImplementedError(f"Unknown engine '{args.engine}'")
    
    workflow.engine = engine

    # Set traceback behavior
    if not args.traceback:
        sys.tracebacklimit = 0
        sys.excepthook = _custom_exception_hook

    # Run workflow
    workflow.run_while()

    # Save workflow to file
    write(workflow, workflow.name + '.pkl')

    # Save the ML model to a separate file
    if workflow.ml.train:
        assert workflow.ml_model is not None
        write(workflow.ml_model, workflow.name + '_ml_model.pkl')
