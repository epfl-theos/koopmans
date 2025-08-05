#!/usr/bin/env python3

"""Main Koopmans CLI."""

import argparse
import json
import logging
import os
import re
import sys
import traceback
from pathlib import Path

import koopmans.mpl_config  # noqa: F401
from koopmans.engines import Engine, LocalhostEngine
from koopmans.io import read
from koopmans.logging_config import setup_logging
from koopmans.utils import print_alert
from koopmans.utils.warnings import configure_warnings

# Disabling the auto-loading of the juliacall ipython extension prior to loading ipdb
# isort: off
os.environ['PYTHON_JULIACALL_AUTOLOAD_IPYTHON_EXTENSION'] = 'no'
import ipdb  # noqa: E402
# isort: on


DEFAULT_ENGINE = 'localhost'
AVAILABLE_ENGINES = ['localhost', 'aiida']


def _custom_exception_hook(exception_type, exception_value, traceback):  # noqa: W0613
    # Adding spaces to the error name so that it is easier to read
    spaced_text = re.sub(r'([a-z])([A-Z])', r'\1 \2', exception_type.__name__)
    spaced_text = re.sub(r'([A-Z]+)([A-Z][a-z])', r'\1 \2', spaced_text)

    print_alert('caution', str(exception_value), header=spaced_text, indent=1)


def _pdb_exception_hook(exception_type, exception_value, exception_traceback):
    traceback.print_exception(exception_type, exception_value, exception_traceback)
    ipdb.post_mortem(exception_traceback)


def initialize_engine(engine_arg: str, engine_config: str | None) -> Engine:
    """Initialize the engine based on the command line arguments."""
    if engine_arg == 'localhost':
        engine = LocalhostEngine()
    elif engine_arg == 'aiida':
        # raise NotImplementedError("AiiDA engine is not yet implemented")
        # Uncomment the following lines once aiida-koopmans is available
        from aiida_koopmans.engine.aiida import AiiDAEngine, AiiDAStepsData
        if engine_config is not None:
            with open(engine_config, 'r') as f:
                engine_config = json.load(f)
        else:
            engine_config = None
        engine = AiiDAEngine(step_data=AiiDAStepsData(configuration=engine_config))
    else:
        raise NotImplementedError(f"Unknown engine '{engine_arg}'")
    return engine


def list_pseudo(namespace: argparse.Namespace):
    """List the available pseudopotential libraries."""
    engine = initialize_engine(namespace.engine, getattr(namespace, 'engine_config', None))

    for p in sorted(engine.available_pseudo_libraries()):
        print(p)


def install_pseudo(namespace: argparse.Namespace):
    """Install a pseudopotential file."""
    engine = initialize_engine(namespace.engine, getattr(namespace, 'engine_config', None))

    for f in namespace.files:
        pseudo_file = Path(f).resolve()
        if not pseudo_file.exists():
            raise FileNotFoundError(f"File {pseudo_file} does not exist")
        engine.install_pseudopotential(pseudo_file, library=namespace.library)


class UninstallPseudoAction(argparse.Action):
    """An action to uninstall a pseudopotential library."""

    def __call__(self, parser, namespace, values, option_string=None):
        """Uninstall a pseudopotential library."""
        engine_name = getattr(namespace, 'engine', DEFAULT_ENGINE)
        engine_config = getattr(namespace, 'engine_config', None)
        engine = initialize_engine(engine_name, engine_config)

        for value in values:
            engine.uninstall_pseudopotential_library(value)


def run_workflow(namespace: argparse.Namespace):
    """Run a Koopmans workflow based on the provided JSON file."""
    # Create the engine
    engine = None
    if namespace.engine:
        engine = initialize_engine(namespace.engine, getattr(namespace, 'engine_config', None))

    # If requested, set up logging
    if namespace.log:
        level = logging.DEBUG if namespace.debug else logging.INFO
        setup_logging(level=level)

    # Configure the warnings
    configure_warnings()

    # Reading in JSON file
    workflow = read(namespace.json, engine=engine)

    # Run workflow
    workflow.run()

    print("Workflow complete 🎉")


def main():
    """Run the main Koopmans CLI."""
    from koopmans import __version__

    # Construct parser
    parser = argparse.ArgumentParser(
        description='Perform Koopmans functional calculations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='See https://koopmans-functionals.org for more details')

    # Subcommands
    subparsers = parser.add_subparsers(dest='command', title='subcommands')

    # koopmans --version
    parser.add_argument('--version', action='version', version=__version__,
                        help="show the program's version number and exit")

    def add_engine_flag(p):
        p.add_argument('--engine', choices=AVAILABLE_ENGINES, default=DEFAULT_ENGINE,
                       help="specify the execution engine")

    # koopmans run
    run_parser = subparsers.add_parser("run", help="run a Koopmans spectral functional calculation")
    run_parser.add_argument('json', metavar='system.json', type=str,
                            help='a single JSON file containing the workflow and code settings')
    run_parser.add_argument('-t', '--traceback', action='store_true', help='enable traceback')
    run_parser.add_argument('-l', '--log', action='store_true', help='enable logging')
    run_parser.add_argument('-d', '--debug', action='store_true', help='enable debug logging')
    run_parser.add_argument('--pdb', action='store_true', help='enable interactive debugging')
    add_engine_flag(run_parser)
    run_parser.add_argument('--engine_config', type=str, default='engine.json',
                            help='Specify the engine configuration file (default: engine.json)')
    run_parser.set_defaults(func=run_workflow, log=False, traceback=False)

    # koopmans pseudos
    pseudos_parser = subparsers.add_parser("pseudos", help="list, install, and uninstall pseudopotentials")
    pseudos_subparsers = pseudos_parser.add_subparsers(title='subcommands')

    # koopmans pseudos list
    pseudos_list = pseudos_subparsers.add_parser("list", help="list available pseudopotential libraries")
    add_engine_flag(pseudos_list)
    pseudos_list.set_defaults(func=list_pseudo)

    # koopmans pseudos install
    pseudos_install = pseudos_subparsers.add_parser("install", help="install a local pseudopotential file")
    pseudos_install.add_argument('files', type=str, help="the local .upf file to install", nargs='+', metavar="file")
    pseudos_install.add_argument('--library', type=str, nargs='?',
                                 help="the custom library to put the pseudopotential in", default="CustomPseudos")
    pseudos_install.set_defaults(func=install_pseudo)
    add_engine_flag(pseudos_install)

    # koopmans pseudos uninstall
    pseudos_uninstall = pseudos_subparsers.add_parser("uninstall", help="uninstall a pseudopotential library")
    pseudos_uninstall.add_argument(
        'library', type=str, help="the pseudopotential library to uninstall", nargs='+', action=UninstallPseudoAction)
    add_engine_flag(pseudos_uninstall)

    # Parse arguments
    args = parser.parse_args()

    # Customising traceback
    if getattr(args, 'pdb', False):
        # Use pdb for debugging
        sys.tracebacklimit = None
        sys.excepthook = _pdb_exception_hook
    elif not getattr(args, 'traceback', False):
        # Use custom traceback behavior by default
        sys.tracebacklimit = 0
        sys.excepthook = _custom_exception_hook

    # Call the action
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
