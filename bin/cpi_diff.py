#!/usr/bin/env python3

import argparse
from ase.io import espresso_cp as cp_io
from koopmans_cp.io import cpi_diff
from koopmans_cp.calculators.cp import CP_calc

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compare two QE input files')
    parser.add_argument('cpis', metavar='template_1.cpi template_2.cpi', type=str, nargs=2,
                        help='the two QE input files')

    args = parser.parse_args()

    calcs = {}
    for cpi in args.cpis:
        calcs[cpi] = CP_calc(cp_io.read_espresso_cp_in(open(cpi, 'r')).calc)

    cpi_diff(calcs)
