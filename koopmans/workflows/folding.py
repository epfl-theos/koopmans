"""

Workflow module for python_KI, containing the workflow for converting W90 files to cp.x friendly-format
using a modified version of pw2wannier.x

Written by Edward Linscott Feb 2021

"""

from koopmans import utils
from koopmans.workflows.generic import Workflow


class FoldToSupercellWorkflow(Workflow):

    def run(self):
        '''

        Wrapper for the calculation of (maximally localized) Wannier functions
        using PW and Wannier90

        '''

        self.print('Folding to supercell', style='subheading')

        for typ in ['occ', 'emp']:
            calc_w2o = self.new_calculator('pw2wannier', wan_mode='wannier2odd', directory=f'./{typ}', name='wan2odd')
            if typ == 'emp':
                calc_w2o.split_evc_file = True
            self.run_calculator(calc_w2o)

        return

    def new_calculator(self, calc_type, *args, **kwargs):
        calc = super().new_calculator(calc_type, *args, **kwargs)

        # Checking that gamma_trick is consistent with do_wf_cmplx
        kcp_calc = self.master_calcs['kcp']
        if calc.gamma_trick == kcp_calc.do_wf_cmplx:
            utils.warn(
                f'if do_wf_cmplx is {kcp_calc.do_wf_cmplx}, gamma_trick cannot be {calc.gamma_trick}. '
                f'Changing gamma_trick to {not kcp_calc.do_wf_cmplx}')
            calc.gamma_trick = not kcp_calc.do_wf_cmplx
        elif calc.gamma_trick is None and not kcp_calc.do_wf_cmplx:
            calc.gamma_trick = True
        else:
            pass
        return calc
