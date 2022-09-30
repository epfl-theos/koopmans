"""

Workflow module for koopmans, containing the workflow for converting W90 or PW files to
kcp.x friendly-format using wann2kcp.x

Written by Edward Linscott Feb 2021

"""

import os
from pathlib import Path

import numpy as np

from koopmans import calculators, utils

from ._workflow import Workflow


class FoldToSupercellWorkflow(Workflow):

    def _run(self):
        '''

        Wrapper for folding Wannier or Kohn-Sham functions from the primitive cell
        to the supercell and convert them to a kcp.x friendly format.

        '''

        self.print('Folding to supercell', style='subheading')

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs']:
            with utils.chdir('wannier'):
                # Loop over the various subblocks that we have wannierized separately
                for block in self.projections:
                    # Create the calculator
                    calc_w2k = self.new_calculator('wann2kcp', spin_component=block.spin, wan_mode='wannier2kcp',
                                                   directory=block.directory)
                    calc_w2k.prefix = 'w2kcp'

                    # Checking that gamma_trick is consistent with gamma_only
                    if calc_w2k.parameters.gamma_trick and not self.kpoints.gamma_only:
                        calc_w2k.parameters.gamma_trick = False
                    elif not calc_w2k.parameters.gamma_trick and self.kpoints.gamma_only:
                        calc_w2k.parameters.gamma_trick = True
                    else:
                        pass

                    # Run the calculator
                    self.run_calculator(calc_w2k)

                if self.parameters.spin_polarized:
                    spins = ['up', 'down']
                else:
                    spins = [None, None]

                # Merging evcw files
                for output_directory, subset in self.projections.to_merge.items():
                    if len(subset) > 1:
                        output_directory.mkdir(exist_ok=True)
                        if self.parameters.spin_polarized:
                            evc_fname = f'evcw.dat'
                        elif 'up' in output_directory.name:
                            evc_fname = f'evcw1.dat'
                        else:
                            evc_fname = f'evcw2.dat'
                        command = ' '.join([f'merge_evc.x -nr {np.prod(self.kpoints.grid)}']
                                           + [f'-i {b.directory}/{evc_fname}' for b in subset]
                                           + [f'-o {output_directory}/{evc_fname}'])
                        if 'occ' in output_directory.name:
                            label = 'occupied'
                        else:
                            label = 'empty'
                        if subset[0].spin:
                            label += f' spin {subset[0].spin}'
                        if self.parameters.from_scratch or not (output_directory / evc_fname).exists():
                            self.parameters.from_scratch = True
                            self.print(f'Merging the {label} band blocks... ', end='')
                            utils.system_call(command)
                            self.print('done')
                        else:
                            self.print(f'Not merging the {label} band blocks as this is already complete')

        else:
            # Create the calculator
            calc_w2k = self.new_calculator('wann2kcp', directory='ks2kcp', wan_mode='ks2kcp', seedname=None)
            calc_w2k.prefix = 'ks2kcp'

            # Run the calculator
            self.run_calculator(calc_w2k)

        return
