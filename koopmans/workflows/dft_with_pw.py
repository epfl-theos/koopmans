"""

Workflow module for performing a single DFT calculation with pw.x

Written by Edward Linscott Mar 2021

"""

import os
import copy
from koopmans import utils
from koopmans.workflows.generic import Workflow


class DFTPWWorkflow(Workflow):

    def run(self):

        # Create the calculator
        calc = self.new_calculator('pw')

        # Update keywords
        calc.name = 'dft'
        calc.ndr = 50
        calc.ndw = 51
        calc.restart_mode = 'from_scratch'

        # Remove old directories
        if self.from_scratch:
            utils.system_call(f'rm -r {calc.outdir} 2>/dev/null', False)

        # Run the calculator
        self.run_calculator(calc)

        return
