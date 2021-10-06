"""

Workflow module for performing a single DFT calculation with either kcp.x or pw.x

Written by Edward Linscott Oct 2020

"""

import os
import copy
from koopmans import utils
from ._generic import Workflow


class DFTCPWorkflow(Workflow):

    def run(self):

        calc = self.new_calculator('kcp')

        # Removing old directories
        if self.parameters.from_scratch:
            utils.system_call(f'rm -r {calc.parameters.outdir} 2>/dev/null', False)

        calc.prefix = 'dft'
        calc.directory = '.'
        calc.parameters.ndr = 50
        calc.parameters.ndw = 51
        calc.parameters.restart_mode = 'from_scratch'
        calc.parameters.do_orbdep = False
        calc.parameters.fixed_state = False
        calc.parameters.do_outerloop = True
        calc.parameters.which_compensation = 'tcc'

        if calc.parameters.maxiter is None:
            calc.parameters.maxiter = 300
        if calc.parameters.empty_states_nbnd is not None and calc.parameters.empty_states_nbnd > 0:
            calc.parameters.do_outerloop_empty = True
            if calc.parameters.empty_states_maxstep is None:
                calc.parameters.empty_states_maxstep = 300

        try:
            self.run_calculator(calc, enforce_ss=self.parameters.enforce_spin_symmetry)
        except:
            raise ValueError()

        return calc


class DFTPWWorkflow(Workflow):

    def run(self):

        # Create the calculator
        calc = self.new_calculator('pw')

        # Update keywords
        calc.prefix = 'dft'
        calc.parameters.ndr = 50
        calc.parameters.ndw = 51
        calc.parameters.restart_mode = 'from_scratch'

        # Remove old directories
        if self.parameters.from_scratch:
            utils.system_call(f'rm -r {calc.parameters.outdir} 2>/dev/null', False)

        # Run the calculator
        self.run_calculator(calc)

        return
