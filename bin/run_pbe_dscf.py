#!/usr/bin/env python3
import argparse
from ase.io import espresso as pw_io
from koopmans.environ import Environ_calc
from koopmans.calculators.calculator import run_qe
import os

'''
Perform delta SCF PBE calculations
'''

if __name__ == '__main__':

    # Construct parser
    parser = argparse.ArgumentParser(
        description='Perform a series of solvated calculations using pw.x')
    parser.add_argument('pwi', metavar='pbe.pwi', type=str,
                        help='an pw.x input file')

    # Parse arguments
    args = parser.parse_args()

    calc = Environ_calc(pw_io.read_espresso_in(args.pwi).calc)

    # Run workflow
    calc_succeeded = True
    calc.name = 'pbe'

    os.system('rm -r neutral charged')

    for charge, label in zip([0, -1], ['neutral', 'charged']):
        # Create working directories
        os.system(f'mkdir {label}')
        os.system(f'mkdir {label}/20')

        # Initialize value of epsilon and other system parameters
        epsilon = 20
        calc.restart_mode = 'from_scratch'
        calc.disk_io = 'medium' # checkpointing files will be required for later restarts
        calc.environ_settings['ENVIRON']['environ_restart'] = False

        # Apply the desired charge
        calc.tot_charge = charge
        calc.tot_magnetization = -charge

        while calc_succeeded:
            calc.directory = f'{label}/{epsilon}'
            calc.environ_settings['ENVIRON']['env_static_permittivity'] = epsilon
            calc.environ_settings['BOUNDARY']['solvent_mode'] = 'ionic'
            # calc.environ_settings['ELECTROSTATIC']['tol'] = 1e-8

            run_qe(calc, silent=False, from_scratch=True)
            # from_scratch = True means that run_qe won't try and skip this calculation
            # if it encounters pre-existing QE output files, and NOT that QE will set
            # restart_mode = 'from_scratch'
            calc_succeeded = calc.is_converged()

            # Preparing for next loop
            if epsilon == 1:
                break
            new_epsilon = max(epsilon - 2, 1)
            os.system(f'cp -r {label}/{epsilon} {label}/{new_epsilon}')
            epsilon = new_epsilon
            calc.restart_mode = 'restart'
            calc.environ_settings['ENVIRON']['environ_restart'] = True
