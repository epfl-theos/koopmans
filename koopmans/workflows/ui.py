"""

Workflow module for koopmans, containing the workflow for unfolding and interpolating
results from the supercell to generate a bandstructure

For the moment the code works only with cubic, tetragonal and orthorombic systems.

Originally written by Riccardo De Genarro as the standalone 'unfolding and interpolate' code
Integrated within koopmans by Edward Linscott Jan 2021
"""

from koopmans.workflows.generic import Workflow


class UnfoldAndInterpolateWorkflow(Workflow):

    def run(self):
        '''
        '''
        ui_calc = self.master_calcs['ui']
        ui_calc.name = self.name
        ui_calc.calculate()
        self.all_calcs = [ui_calc]

        # Print quality control
        if self.print_qc and not ui_calc.skip_qc:
            for key in ui_calc.results_for_qc:
                val = ui_calc.results.get(key, None)
                if val is not None:
                    ui_calc.qc_results[key] = val
