"""

singlepoint workflow object for python_KI

Written by Edward Linscott Oct 2020
Converted to a workflow object Nov 2020

"""

import os
import copy
import numpy as np
import itertools
from koopmans import utils, io
from koopmans.workflows.generic import Workflow
from koopmans.workflows.kc_with_cp import KoopmansWorkflow
from koopmans.workflows.pbe_with_cp import PBEWorkflow


class SinglepointWorkflow(Workflow):

    def run(self):

        if self.functional == 'all':
            # if 'all', create subdirectories and run
            functionals = ['ki', 'pkipz', 'kipz']

            # Make separate directories for KI, pKIPZ, and KIPZ
            for functional in functionals:
                if self.from_scratch and os.path.isdir(functional):
                    utils.system_call(f'rm -r {functional}')
                if not os.path.isdir(functional):
                    utils.system_call(f'mkdir {functional}')

            if self.alpha_from_file:
                utils.system_call('cp file_alpharef*.txt ki/')

            alphas = None
            for functional in functionals:
                print(
                    f'\n{functional.upper().replace("PKIPZ", "pKIPZ")} CALCULATION')

                # Make a copy of the workflow settings to modify
                local_workflow_settings = self.settings

                local_workflow_settings['functional'] = functional

                # For pKIPZ/KIPZ, use KI as a starting point
                if functional == 'pkipz':
                    local_workflow_settings['from_scratch'] = False
                    local_workflow_settings['calculate_alpha'] = False
                    local_workflow_settings['alpha_from_file'] = False
                elif functional == 'kipz':
                    local_workflow_settings['init_density'] = 'ki'
                    local_workflow_settings['init_variational_orbitals'] = 'ki'
                    local_workflow_settings['alpha_from_file'] = True

                # Change to relevant subdirectory
                os.chdir(functional)

                # Create a KC workflow for this particular functional
                kc_workflow = KoopmansWorkflow(local_workflow_settings, self.master_calcs, alphas)

                # Run the workflow
                kc_workflow.run()

                # Save the alpha values and the final calculation of the workflow
                alphas = kc_workflow.alpha_df.iloc[-1].values
                solved_calc = kc_workflow.all_calcs[-1]

                # Return to the base directory
                os.chdir('..')

                # Provide the pKIPZ and KIPZ calculations with a KI starting point
                if functional == 'ki':
                    # pKIPZ
                    utils.system_call('rsync -a ki/ pkipz/')

                    # KIPZ
                    utils.system_call('rsync -a ki/final/ kipz/init/')
                    utils.system_call(f'rsync -a {solved_calc.outdir}/ kipz/')
                    utils.system_call(
                        'mv kipz/init/ki_final.cpi kipz/init/ki_init.cpi')
                    utils.system_call(
                        'mv kipz/init/ki_final.cpo kipz/init/ki_init.cpo')
            return solved_calc

        else:
            if self.functional in ['ki', 'kipz', 'pkipz']:
                workflow = KoopmansWorkflow(self.settings, self.master_calcs)

            elif self.functional == 'pbe':
                workflow = PBEWorkflow(self.settings, self.master_calcs)

            return workflow.run()
