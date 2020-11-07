"""

Workflow module for performing a single PBE calculation with cp.x

Written by Edward Linscott Oct 2020

"""

import os
import copy
from koopmans import utils, io
from koopmans.calculators.calculator import run_qe

def run(workflow_settings, calcs_dct):
    '''
    This function runs a single PBE calculation with cp.x

    Arguments:
        workflow_settings -- a dictionary containing workflow settings
        calcs_dct         -- a dictionary of calculators, including a CP_calc calculator
                             object containing cp.x settings

    '''

    if 'cp' not in calcs_dct:
        raise ValueError('"functional": "PBE" requires a "cp" block in the input .json file')

    calc = copy.deepcopy(calcs_dct['cp'])

    # Sanitise outdir and define 'outdir'
    calc.outdir = calc.outdir.strip('./')
    if '/' in calc.outdir:
        raise ValueError('"outdir" cannot be a nested directory')
    calc.outdir = os.getcwd() + '/' + calc.outdir

    # Removing old directories
    if workflow_settings['from_scratch']:
        utils.system_call(f'rm -r {calc.outdir} 2>/dev/null', False)

    calc.name = 'pbe'
    calc.ndr = 50
    calc.ndw = 51
    calc.restart_mode = 'from_scratch'
    calc.do_orbdep = False
    calc.fixed_state = False
    calc.do_outerloop = True
    calc.which_compensation = 'tcc'

    if calc.maxiter is None:
        calc.maxiter = 300
    if calc.empty_states_nbnd > 0:
        calc.do_outerloop_empty = True
        if calc.empty_states_maxstep is None:
            calc.empty_states_maxstep = 300

    run_qe(calc, silent=False, from_scratch=workflow_settings['from_scratch'],
               enforce_ss=workflow_settings['enforce_spin_symmetry'])

    if workflow_settings['print_qc']:
        for var in ['energy', 'homo_energy']:
            io.print_qc(var, calc.results[var])

    return calc
