"""

Workflow module for performing a single DFT calculation with kcp.x

Written by Edward Linscott Oct 2020

"""

import os
import copy
from koopmans import utils
from koopmans.workflows.generic import Workflow


class DFTWorkflow(Workflow):

    def __init__(self, *args):
        super().__init__(*args)
        if 'kcp' not in self.master_calcs:
            raise ValueError(
                '"functional": "PBE" requires a "kcp" block in the input .json file')

    def run(self):

        calc = self.new_calculator('kcp')

        # Removing old directories
        if self.from_scratch:
            utils.system_call(f'rm -r {calc.outdir} 2>/dev/null', False)

        calc.name = 'dft'
        calc.directory = '.'
        calc.ndr = 50
        calc.ndw = 51
        calc.restart_mode = 'from_scratch'
        calc.do_orbdep = False
        calc.fixed_state = False
        calc.do_outerloop = True
        calc.which_compensation = 'tcc'

        if calc.maxiter is None:
            calc.maxiter = 300
        if calc.empty_states_nbnd is not None and calc.empty_states_nbnd > 0:
            calc.do_outerloop_empty = True
            if calc.empty_states_maxstep is None:
                calc.empty_states_maxstep = 300

        self.run_calculator(calc, enforce_ss=self.enforce_spin_symmetry)

        return calc
