#!/usr/bin/env python3

import argparse
from ase.io import espresso_cp as cp_io
from koopmans_cp.io import write_json, read_json
from koopmans_cp.defaults import defaults
from koopmans_cp.workflow import keywords_altered_during_workflow
import copy
import textwrap


def cpi_to_json(cpi, json, to_exclude=['nat', 'ntyp', 'pseudo_dir'], calc_param={}):
    '''

    Converts a cpi file to a json file

    Arguments
    ---------
        cpi: the name of the .cpi file to read in
        json: the name of the .json file to write out to
        to_exclude: the keywords included in the .cpi file to exclude from the .json file
        calc_param: the koopmans.py keywords to add to the .json file (these are not
                    included in the .cpi file)

    '''
    calc = cp_io.read_espresso_cp_in(open(cpi, 'r')).calc

    to_exclude += keywords_altered_during_workflow
    to_exclude += list(defaults.keys())

    calc_out = copy.deepcopy(calc)

    for blockname, block in calc.parameters['input_data'].items():
        for key in block:
            if key in to_exclude:
                del calc_out.parameters['input_data'][blockname][key]

    write_json(json, calc_out, calc_param=calc_param)

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Converts .cpi files to .json files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
               additional arguments:
                 -<keyword>  <value>  Add '<keyword>: <value>' to the calc_params
                                      block of the JSON file
                 '''))
    parser.add_argument('cpi', metavar='in.cpi', type=str,
                        help='the QE input file to read in')
    parser.add_argument('json', metavar='out.json', type=str,
                        help='the JSON input file to write out')

    parsed, unknown = parser.parse_known_args()

    # Find out the additional keywords provided to the parser
    for arg in unknown:
        if arg.startswith(('-', '--')):
            parser.add_argument(arg)

    args = parser.parse_args()

    settings = {"calc_type": "ki",
                "init_manifold": "pz",
                "n_max_sc_steps": 2,
                "alpha_from_file": "false"}

    # Store any additional arguments in the settings dict
    for arg, value in args.__dict__.items():
        if arg in parsed:
            continue
        settings[arg] = value

    cpi_to_json(args.cpi, args.json, calc_param=settings)
