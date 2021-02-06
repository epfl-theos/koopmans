"""

Workflow module for python_KI, containing the workflow for unfolding and interpolating
results from the supercell to generate a bandstructure

For the moment the code works only with cubic, tetragonal and orthorombic systems.

Originally written by Riccardo De Genarro as the standalone 'unfolding and interpolate' code
Integrated within python_KI by Edward Linscott Jan 2021
"""

from koopmans.workflows.generic import Workflow


class UnfoldAndInterpolateWorkflow(Workflow):

    def __init__(self, workflow_settings, calcs_dct):
        super().__init__(workflow_settings, calcs_dct)

    def run(self):
        '''
        '''
        ui_calc = self.master_calcs['ui']
        ui_calc.calculate()
