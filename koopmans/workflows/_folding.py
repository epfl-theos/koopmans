"""

Workflow module for koopmans, containing the workflow for converting W90 or PW files to
kcp.x friendly-format using a modified version of pw2wannier90.x

Written by Edward Linscott Feb 2021

"""

import numpy as np
from pathlib import Path
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
            # Loop over the various subblocks that we have wannierised separately
            for block in self.projections:
                # Create the calculator
                calc_w2o = self.new_calculator('pw2wannier', spin_component=block.spin, wan_mode='wannier2odd',
                                               directory=block.directory)
                calc_w2o.prefix = 'wan2odd'

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

            if self.parameters.spin_polarised:
                spins = ['up', 'down']
            else:
                spins = [None, None]

            # Merging evcw files
            for occ in [True, False]:
                for spin, evc_index in zip(spins, [1, 2]):
                    subset = self.projections.get_subset(occ=occ, spin=spin)

                    if len(subset) > 1:
                        output_directory = Path(subset[0].merge_directory)
                        output_directory.mkdir(exist_ok=True)
                        if self.parameters.spin_polarised:
                            evc_fname = f'evcw.dat'
                        else:
                            evc_fname = f'evcw{evc_index}.dat'
                        command = ' '.join([f'{calculators.qe_bin_directory}/merge_evc.x -nr {np.prod(self.kgrid)}']
                                           + [f'-i {b.directory}/{evc_fname}' for b in subset]
                                           + [f'-o {output_directory}/{evc_fname}'])
                        if occ:
                            label = 'occupied'
                        else:
                            label = 'empty'
                        label += f' spin {evc_index}'
                        if self.parameters.from_scratch or not (output_directory / evc_fname).exists():
                            self.parameters.from_scratch = True
                            self.print(f'Merging the {label} band blocks... ', end='')
                            utils.system_call(command)
                            self.print('done')
                        else:
                            self.print(f'Not merging the {label} band blocks as this is already complete')

        else:
            # Create the calculator
            calc_w2o = self.new_calculator('pw2wannier', directory='ks2odd', wan_mode='ks2odd', seedname=None)
            calc_w2o.prefix = 'ks2odd'

            # Run the calculator
            self.run_calculator(calc_w2o)

        return
