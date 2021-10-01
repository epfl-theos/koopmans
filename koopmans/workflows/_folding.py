"""

Workflow module for koopmans, containing the workflow for converting W90 files to cp.x friendly-format
using a modified version of pw2wannier.x

Written by Edward Linscott Feb 2021

"""

import numpy as np
from koopmans import utils
from koopmans import calculators
from ._generic import Workflow


class FoldToSupercellWorkflow(Workflow):

    def run(self):
        '''

        Wrapper for the calculation of (maximally localized) Wannier functions
        using PW and Wannier90

        '''

        self.print('Folding to supercell', style='subheading')

        for typ in ['occ', 'emp']:
            # Create the calculator
            kwargs = self.master_calc_params['pw2wannier'].copy()
            kwargs['wan_mode'] = 'wannier2odd'
            calc_w2o = calculators.PW2WannierCalculator(atoms=self.atoms,  **kwargs)
            calc_w2o.directory = f'./{typ}'
            calc_w2o.prefix = 'wan2odd'

            # Checking that gamma_trick is consistent with do_wf_cmplx
            kcp_params = self.master_calc_params['kcp']
            if calc_w2o.parameters.gamma_trick == kcp_params.do_wf_cmplx:
                utils.warn(
                    f'if do_wf_cmplx is {kcp_params.do_wf_cmplx}, gamma_trick cannot be {calc_w2o.parameters.gamma_trick}. '
                    f'Changing gamma_trick to {not kcp_params.do_wf_cmplx}')
                calc_w2o.parameters.gamma_trick = not kcp_params.do_wf_cmplx
            elif calc_w2o.parameters.gamma_trick is None and not kcp_params.do_wf_cmplx:
                calc_w2o.parameters.gamma_trick = True

            if typ == 'emp':
                calc_w2o.parameters.split_evc_file = True

            # Run the calculator
            self.run_calculator(calc_w2o)

        return
