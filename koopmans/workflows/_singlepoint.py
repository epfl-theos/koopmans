"""

singlepoint workflow object for koopmans

Written by Edward Linscott Oct 2020
Converted to a workflow object Nov 2020

"""

import os
import copy
import numpy as np
from koopmans import utils
from ._generic import Workflow


load_results_from_output = True


class SinglepointWorkflow(Workflow):

    def run(self) -> None:

        # Import it like this so if they have been monkey-patched, we will get the monkey-patched version
        from koopmans.workflows import KoopmansDFPTWorkflow, KoopmansDSCFWorkflow, DFTCPWorkflow

        if self.parameters.method == 'dfpt':
            workflow = KoopmansDFPTWorkflow(**self.wf_kwargs)
            self.run_subworkflow(workflow)

        elif self.parameters.functional == 'all':
            # if 'all', create subdirectories and run
            functionals = ['ki', 'pkipz', 'kipz']

            # Make separate directories for KI, pKIPZ, and KIPZ
            for functional in functionals:
                if self.parameters.from_scratch and os.path.isdir(functional):
                    utils.system_call(f'rm -r {functional}')
                if not os.path.isdir(functional):
                    utils.system_call(f'mkdir {functional}')

            if self.parameters.alpha_from_file:
                utils.system_call('cp file_alpharef*.txt ki/')

            for functional in functionals:
                self.print(f'\n{functional.upper().replace("PKIPZ", "pKIPZ")} calculation', style='heading')

                # Make a copy of the workflow settings to modify
                wf_kwargs = self.wf_kwargs
                local_parameters = wf_kwargs.pop('parameters')

                # Select the functional
                local_parameters.functional = functional

                # For pKIPZ/KIPZ, use KI as a starting point
                if functional == 'pkipz':
                    local_parameters.calculate_alpha = False
                restart_from_old_ki = (functional == 'kipz')

                # We only need to do the smooth interpolation the first time (i.e. for KI)
                redo_smooth_dft = None if functional == 'ki' else False

                # Create a KC workflow for this particular functional
                kc_workflow = KoopmansDSCFWorkflow(parameters=local_parameters,
                                                   restart_from_old_ki=restart_from_old_ki,
                                                   redo_smooth_dft=redo_smooth_dft, **wf_kwargs)

                # Transform to the supercell
                if functional == 'kipz':
                    kc_workflow.primitive_to_supercell()

                # Run the workflow
                if functional == 'pkipz' and self.parameters.from_scratch:
                    # We want to run pKIPZ with from_scratch = False, but don't want this to be inherited
                    self.run_subworkflow(kc_workflow, subdirectory=functional, from_scratch=False)
                else:
                    self.run_subworkflow(kc_workflow, subdirectory=functional)

                # Provide the pKIPZ and KIPZ calculations with a KI starting point
                if functional == 'ki':
                    # files are not overwritten as long as from_scratch is false
                    if self.parameters.from_scratch:
                        rsync_cmd = 'rsync -a'
                    else:
                        rsync_cmd = 'rsync -a --ignore-existing'

                    # pKIPZ
                    utils.system_call(f'{rsync_cmd} ki/ pkipz/')

                    # KIPZ
                    utils.system_call(f'{rsync_cmd} ki/final/ kipz/init/')
                    utils.system_call('mv kipz/init/ki_final.cpi kipz/init/ki_init.cpi')
                    utils.system_call('mv kipz/init/ki_final.cpo kipz/init/ki_init.cpo')
                    if self.parameters.init_orbitals in ['mlwfs', 'projwfs']:
                        utils.system_call(f'{rsync_cmd} ki/init/wannier kipz/init/')
                    if self.parameters.periodic and self.master_calc_params['ui'].do_smooth_interpolation:
                        # Copy over the smooth PBE calculation from KI for KIPZ to use
                        utils.system_call(f'{rsync_cmd} ki/postproc kipz/')
                        if self.parameters.from_scratch:
                            utils.system_call('find kipz/postproc/ -name "*interpolated.dat" -delete')

        else:
            # self.functional != all and self.method != 'dfpt'
            if self.parameters.functional in ['ki', 'pkipz', 'kipz']:
                dscf_workflow = KoopmansDSCFWorkflow(**self.wf_kwargs)
                self.run_subworkflow(dscf_workflow)
            else:
                dft_workflow = DFTCPWorkflow(**self.wf_kwargs)
                self.run_subworkflow(dft_workflow)
