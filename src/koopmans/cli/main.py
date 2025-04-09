#!/usr/bin/env python3

"""Main Koopmans CLI."""

import argparse
import re
import sys
from pathlib import Path

import koopmans.mpl_config  # noqa: F401
from koopmans.engines import Engine, LocalhostEngine
from koopmans.io import read
from koopmans.logging_config import setup_logging
from koopmans.utils import print_alert

DEFAULT_ENGINE = 'localhost'


def _custom_exception_hook(exception_type, exception_value, traceback):  # noqa: W0613
    # Adding spaces to the error name so that it is easier to read
    spaced_text = re.sub(r'([a-z])([A-Z])', r'\1 \2', exception_type.__name__)
    spaced_text = re.sub(r'([A-Z]+)([A-Z][a-z])', r'\1 \2', spaced_text)

    print_alert('caution', str(exception_value), header=spaced_text, indent=1)


class ListPseudoAction(argparse.Action):
    """An action to list the available pseudopotential libraries."""

    def __call__(self, parser, namespace, values, option_string=None):
        """List the available pseudopotential libraries."""
        engine_name = getattr(namespace, 'engine', DEFAULT_ENGINE)
        engine_config = getattr(namespace, 'engine_config', None)
        engine = initialize_engine(engine_name, engine_config)

        for p in sorted(engine.available_pseudo_libraries()):
            print(p)


class InstallPseudoAction(argparse.Action):
    """An action to install a pseudopotential file."""

    def __call__(self, parser, namespace, values, option_string=None):
        """Install a pseudopotential file."""
        engine_name = getattr(parser, 'engine', DEFAULT_ENGINE)
        engine_config = getattr(parser, 'engine_config', None)
        engine = initialize_engine(engine_name, engine_config)

        for f in parser.file:
            pseudo_file = Path(f).resolve()
            if not pseudo_file.exists():
                raise FileNotFoundError(f"File {pseudo_file} does not exist")
            engine.install_pseudopotential(pseudo_file, library=parser.library)


class UninstallPseudoAction(argparse.Action):
    """An action to uninstall a pseudopotential library."""

    def __call__(self, parser, namespace, values, option_string=None):
        """Uninstall a pseudopotential library."""
        engine_name = getattr(namespace, 'engine', DEFAULT_ENGINE)
        engine_config = getattr(namespace, 'engine_config', None)
        engine = initialize_engine(engine_name, engine_config)

        for value in values:
            engine.uninstall_pseudopotential_library(value)


def initialize_engine(engine_arg: str, engine_config: str | None) -> Engine:
    """Initialize the engine based on the command line arguments."""
    if engine_arg == 'localhost':
        engine = LocalhostEngine()
    elif engine_arg == 'aiida':
        raise NotImplementedError("AiiDA engine is not yet implemented")
        # Uncomment the following lines once aiida_koopmans is available
        # from aiida_koopmans.engine.aiida import AiiDAEngine
        # if engine_config is not None:
        #     with open(engine_config, 'r') as f:
        #         engine_config = json.load(f)
        # else:
        #     engine_config = None
        # engine = AiiDAEngine(configuration=engine_config)
    else:
        raise NotImplementedError(f"Unknown engine '{engine_arg}'")
    return engine


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
        p.add_argument('--engine', choices=['localhost', 'aiida'], default=DEFAULT_ENGINE,
                       help="specify the execution engine")

    # koopmans run
    run_parser = subparsers.add_parser("run", help="run a Koopmans spectral functional calculation")
    run_parser.add_argument('json', metavar='system.json', type=str,
                            help='a single JSON file containing the workflow and code settings')
    run_parser.add_argument('-t', '--traceback', action='store_true', help='enable traceback')
    run_parser.add_argument('-l', '--log', action='store_true', help='enable logging')
    add_engine_flag(run_parser)
    run_parser.add_argument('--engine_config', type=str, default='engine.json',
                            help='Specify the engine configuration file (default: engine.json)')

    # koopmans pseudos
    pseudos_parser = subparsers.add_parser("pseudos", help="list, install, and uninstall pseudopotentials")
    pseudos_subparsers = pseudos_parser.add_subparsers(title='subcommands')

    # koopmans pseudos list
    pseudos_list = pseudos_subparsers.add_parser("list", help="list available pseudopotential libraries")
    pseudos_list.set_defaults(action=ListPseudoAction)
    add_engine_flag(pseudos_list)

    # koopmans pseudos install
    pseudos_install = pseudos_subparsers.add_parser("install", help="install a local pseudopotential file")
    pseudos_install.add_argument('file', type=str, help="the local .upf file to install", nargs='+')
    pseudos_install.add_argument('--library', type=str, nargs='?',
                                 help="the custom library to put the pseudopotential in", default="CustomPseudos")
    pseudos_install.set_defaults(action=InstallPseudoAction)
    add_engine_flag(pseudos_install)

    # koopmans pseudos uninstall
    pseudos_uninstall = pseudos_subparsers.add_parser("uninstall", help="uninstall a pseudopotential library")
    pseudos_uninstall.add_argument(
        'library', type=str, help="the pseudopotential library to uninstall", nargs='+', action=UninstallPseudoAction)
    add_engine_flag(pseudos_uninstall)

    # Hide traceback
    sys.tracebacklimit = 0
    default_excepthook, sys.excepthook = sys.excepthook, _custom_exception_hook

    # Parse arguments
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        parser.exit()

    # Create the engine
    engine = None
    if getattr(args, 'engine', None):
        engine = initialize_engine(args.engine, getattr(args, 'engine_config', None))

    # For koopmans pseudo list, perform the action and exit
    if args.command == 'pseudos':
        if hasattr(args, 'action'):
            args.action(parser, args, None, None)
        else:
            pseudos_parser.print_help()
        parser.exit()

    # Restore traceback behavior if requested
    if args.traceback:
        sys.tracebacklimit = None
        sys.excepthook = default_excepthook

    # If requested, set up logging
    if args.log:
        setup_logging()

    # Reading in JSON file
    workflow = read(args.json, engine=engine)

    # Run workflow
    workflow.run()
