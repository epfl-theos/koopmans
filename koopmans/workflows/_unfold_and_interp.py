"""

Workflow module for koopmans, containing the workflow for unfolding and interpolating
results from the supercell to generate a bandstructure

For the moment the code works only with cubic, tetragonal and orthorhombic systems.

Originally written by Riccardo De Genarro as the standalone 'unfolding and interpolate' code
Integrated within koopmans by Edward Linscott Jan 2021
"""

from ._generic import Workflow


class UnfoldAndInterpolateWorkflow(Workflow):

    def run(self):
        '''
        '''
        ui_calc = self.new_calculator('ui')
        ui_calc.prefix = self.name

        ui_calc.calculate()
        self.calculations = [ui_calc]

        # Print quality control
        if self.parameters.print_qc and not ui_calc.skip_qc:
            for key in ui_calc.results_for_qc:
                val = ui_calc.results.get(key, None)
                if val is not None:
                    ui_calc.qc_results[key] = val
