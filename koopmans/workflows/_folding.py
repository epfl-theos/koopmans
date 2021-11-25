"""

Workflow module for koopmans, containing the workflow for converting W90 or PW files to 
kcp.x friendly-format using a modified version of pw2wannier90.x

Written by Edward Linscott Feb 2021

"""

import itertools
from koopmans import utils
from koopmans import calculators
from ._generic import Workflow


class FoldToSupercellWorkflow(Workflow):

    def run(self):
        '''

        Wrapper for folding Wannier or Kohn-Sham functions from the primitive cell
        to the supercell and convert them to a kcp.x dreindly format.

        '''

        self.print('Folding to supercell', style='subheading')

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs']:
            fillings = ['occ', 'emp']
            if self.parameters.spin_polarised:
                spins = ['up', 'down']
            else:
                spins = [None]

            for filling, spin in itertools.product(fillings, spins):
                if spin is None:
                    typ = filling
                else:
                    typ = f'{filling}_{spin}'

                # Create the calculator
                kwargs = self.master_calc_params['pw2wannier'].copy()
                kwargs['wan_mode'] = 'wannier2odd'
                calc_w2o = calculators.PW2WannierCalculator(atoms=self.atoms, **kwargs)
                calc_w2o.directory = f'./{typ}'
                calc_w2o.prefix = 'wan2odd'
                if spin is not None:
                    calc_w2o.spin_component = spin

                # Checking that gamma_trick is consistent with do_wf_cmplx
                kcp_params = self.master_calc_params['kcp']
                if calc_w2o.parameters.gamma_trick == kcp_params.do_wf_cmplx:
                    utils.warn(
                        f'if do_wf_cmplx is {kcp_params.do_wf_cmplx}, gamma_trick cannot be '
                        f'{calc_w2o.parameters.gamma_trick}. Changing gamma_trick to {not kcp_params.do_wf_cmplx}')
                    calc_w2o.parameters.gamma_trick = not kcp_params.do_wf_cmplx
                elif calc_w2o.parameters.gamma_trick is None and not kcp_params.do_wf_cmplx:
                    calc_w2o.parameters.gamma_trick = True

                # Run the calculator
                self.run_calculator(calc_w2o)
        
        else:
            # Create the calculator
            kwargs = self.master_calc_params['pw2wannier'].copy()
            kwargs['wan_mode'] = 'ks2odd'
            calc_w2o = calculators.PW2WannierCalculator(atoms=self.atoms, **kwargs)
            calc_w2o.directory = 'ks2odd'
            calc_w2o.prefix = 'ks2odd'
            del calc_w2o.parameters['seedname']

            # Run the calculator
            self.run_calculator(calc_w2o)

        return
