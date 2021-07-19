#!/usr/bin/env python3

import argparse
from koopmans.io import write_json
from koopmans.defaults import defaults
from koopmans.calculators.cp import CP_calc
import copy
import textwrap


def cpi_to_json(cpi, json, to_exclude=['nat', 'ntyp'], workflow_settings={}, cp_param={}):
    '''

    Converts a cpi file to a json file

    Arguments
    ---------
        cpi: the name of the .cpi file to read in
        json: the name of the .json file to write out to
        to_exclude: the keywords included in the .cpi file to exclude from the .json file
        workflow_settings: the koopmans.py keywords to add to the .json file (these are not
                    included in the .cpi file)
        cp_param: cp flags to alter

    '''
    calc = CP_calc(qe_files=cpi)

    to_exclude += list(defaults[CP_calc].keys())

    calc_out = copy.deepcopy(calc)

    for key in to_exclude:
        setattr(calc_out, key, None)
    for key, value in cp_param.items():
        setattr(calc_out, key, value)

    write_json(json, calc_out, workflow_settings=workflow_settings)

    return calc_out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Converts .cpi files to .json files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
               additional arguments:
                 -<keyword>  <value>  Add '<keyword>: <value>' to the workflow_settingss
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

    cpi_to_json(args.cpi, args.json, workflow_settings=settings)
