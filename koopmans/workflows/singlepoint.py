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
from koopmans.calculators.kcp import KCP_calc
from koopmans.workflows.generic import Workflow
from koopmans.workflows.kc_with_cp import KoopmansWorkflow
from koopmans.workflows.dft_with_cp import DFTWorkflow


load_results_from_output = True


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
                print(f'\n{functional.upper().replace("PKIPZ", "pKIPZ")} CALCULATION')

                # Make a copy of the workflow settings to modify
                local_workflow_settings = self.settings

                local_workflow_settings['functional'] = functional

                # For pKIPZ/KIPZ, use KI as a starting point
                if functional == 'pkipz':
                    local_workflow_settings['calculate_alpha'] = False
                    local_workflow_settings['alpha_from_file'] = False
                elif functional == 'kipz':
                    local_workflow_settings['init_orbitals'] = 'from old ki'
                    local_workflow_settings['alpha_from_file'] = True

                # Create a KC workflow for this particular functional
                master_calcs_local = copy.deepcopy(self.master_calcs)
                kc_workflow = KoopmansWorkflow(local_workflow_settings, master_calcs_local, alphas)

                # We only need to do the smooth interpolation the first time (i.e. for KI)
                if functional != 'ki':
                    kc_workflow.redo_preexisting_smooth_dft_calcs = False

                # Run the workflow
                if functional == 'pkipz' and self.from_scratch:
                    # We want to run pKIPZ with from_scratch = False, but don't want this to be inherited
                    self.run_subworkflow(kc_workflow, subdirectory=functional, from_scratch=False)
                else:
                    self.run_subworkflow(kc_workflow, subdirectory=functional)

                # Save the alpha values and the calculations of the workflow
                alphas = kc_workflow.alpha_df.iloc[-1].values

                # Provide the pKIPZ and KIPZ calculations with a KI starting point
                if functional == 'ki':
                    # pKIPZ
                    utils.system_call('rsync -a ki/ pkipz/')

                    # KIPZ
                    utils.system_call('rsync -a ki/final/ kipz/init/')
                    utils.system_call('mv kipz/init/ki_final.cpi kipz/init/ki_init.cpi')
                    utils.system_call('mv kipz/init/ki_final.cpo kipz/init/ki_init.cpo')
                    if self.init_orbitals in ['mlwfs', 'projw']:
                        utils.system_call('rsync -a ki/init/wannier kipz/init/')
                    if self.master_calcs['ui'].do_smooth_interpolation:
                        # Copy over the smooth PBE calculation from KI for KIPZ to use
                        utils.system_call('rsync -a ki/postproc kipz/')
                        utils.system_call('find kipz/postproc/ -name "*interpolated.dat" -delete')
            return

        else:
            if self.functional in ['ki', 'kipz', 'pkipz']:
                workflow = KoopmansWorkflow(self.settings, self.master_calcs)

            else:
                workflow = DFTWorkflow(self.settings, self.master_calcs)

            self.run_subworkflow(workflow)
            return
