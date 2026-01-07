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

import yaml
from pydantic import ValidationError

import koopmans.mpl_config  # noqa: F401
from koopmans.cli.input_schema.input_file import KoopmansInput, prettify_errors
from koopmans.engines.ase import LocalASEEngine
from koopmans.engines.engine import Engine
from koopmans.io import read
from koopmans.logging_config import setup_logging
from koopmans.utils import indented_print, print_alert
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
        if Path(engine_config).exists():
            with open(engine_config, 'r') as f:
                engine_config = json.load(f)
        else:
            engine_config = {}
        engine = LocalASEEngine.model_validate(engine_config)
    elif engine_arg == 'aiida':
        # raise NotImplementedError("AiiDA engine is not yet implemented")
        # Uncomment the following lines once aiida-koopmans is available
        from aiida_koopmans.engine import AiiDAEngine, AiiDAStepsData
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
    engine = initialize_engine(namespace.engine, namespace.engine_config)

    # If requested, set up logging
    if namespace.log:
        level = logging.DEBUG if namespace.debug else logging.INFO
        setup_logging(level=level)

    # Configure the warnings
    configure_warnings()

    # Print the header
    indented_print(header(), indent=0)

    # Reading in input file
    try:
        koopmans_input = KoopmansInput.from_file(namespace.system)
    except ValidationError as e:
        raise ValueError(
            f"Errors found in the input file: \n\n{prettify_errors(e)}\n\n For more information, see URL_HERE.")
    workflow = koopmans_input.to_workflow(engine=engine)

    # Run workflow
    workflow.run()

    print("Workflow complete 🎉")


def header():
    """Print the output header."""
    from koopmans import __version__

    bf = '**' if sys.stdout.isatty() else ''

    header = ["",
              bf + "koopmans" + bf,
              bf + "========" + bf,
              "",
              "*Koopmans spectral functional calculations with `Quantum ESPRESSO`*",
              "",
              f"📦 **Version:** {__version__}  ",
              "🧑 **Authors:** Edward Linscott, Nicola Colonna, Riccardo De Gennaro, Ngoc Linh Nguyen, Giovanni "
              "Borghi, Andrea Ferretti, Ismaila Dabo, and Nicola Marzari  ",
              "📚 **Documentation:** https://koopmans-functionals.org  ",
              "❓ **Support:** https://groups.google.com/g/koopmans-users  ",
              "🐛 **Report a bug:** https://github.com/epfl-theos/koopmans/issues/new"
              ]

    return '\n'.join(header)


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
    run_parser.add_argument('system', metavar='system.json', type=str,
                            help='a JSON file containing the system details')
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
