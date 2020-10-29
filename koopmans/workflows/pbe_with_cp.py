"""

Workflow module for performing a single PBE calculation with cp.x

Written by Edward Linscott Oct 2020

"""

import os
import copy
from koopmans import utils, io
from koopmans.calculators.calculator import run_qe

def run(master_calc, workflow_settings):
    '''
    This function runs a single PBE calculation with cp.x

    Arguments:
        master_calc       -- a CP_calc object containing cp.x settings
        workflow_settings -- a dictionary containing workflow settings

    '''

    # Sanitise outdir and define 'outdir'
    master_calc.outdir = master_calc.outdir.strip('./')
    if '/' in master_calc.outdir:
        raise ValueError('"outdir" cannot be a nested directory')
    master_calc.outdir = os.getcwd() + '/' + master_calc.outdir

    # Removing old directories
    if workflow_settings['from_scratch']:
        utils.system_call(f'rm -r {master_calc.outdir} 2>/dev/null', False)

    master_calc = copy.deepcopy(master_calc)
    master_calc.name = 'pbe'
    master_calc.ndr = 50
    master_calc.ndw = 51
    master_calc.restart_mode = 'from_scratch'
    master_calc.do_orbdep = False
    master_calc.fixed_state = False
    master_calc.do_outerloop = True
    master_calc.which_compensation = 'tcc'

    if master_calc.maxiter is None:
        master_calc.maxiter = 300
    if master_calc.empty_states_nbnd > 0:
        master_calc.do_outerloop_empty = True
        if master_calc.empty_states_maxstep is None:
            master_calc.empty_states_maxstep = 300

    run_qe(master_calc, silent=False, from_scratch=workflow_settings['from_scratch'],
               enforce_ss=workflow_settings['enforce_spin_symmetry'])

    if workflow_settings['print_qc']:
        for var in ['energy', 'homo_energy']:
            io.print_qc(var, master_calc.results[var])

    return master_calc
