"""

singlepoint workflow object for koopmans

Written by Edward Linscott Oct 2020
Converted to a workflow object Nov 2020

"""

import os
import copy
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
                local_workflow_settings = wf_kwargs.pop('workflow_settings')

                # Select the functional
                local_workflow_settings.functional = functional

                # For pKIPZ/KIPZ, use KI as a starting point
                if functional == 'pkipz':
                    local_workflow_settings.calculate_alpha = False
                elif functional == 'kipz':
                    local_workflow_settings.init_orbitals = 'from old ki'

                # Create a KC workflow for this particular functional
                kc_workflow = KoopmansDSCFWorkflow(workflow_settings=local_workflow_settings, **wf_kwargs)

                # We only need to do the smooth interpolation the first time (i.e. for KI)
                if functional != 'ki':
                    kc_workflow.redo_preexisting_smooth_dft_calcs = False

                # Run the workflow
                if functional == 'pkipz' and self.parameters.from_scratch:
                    # We want to run pKIPZ with from_scratch = False, but don't want this to be inherited
                    self.run_subworkflow(kc_workflow, subdirectory=functional, from_scratch=False)
                else:
                    self.run_subworkflow(kc_workflow, subdirectory=functional)

                # Provide the pKIPZ and KIPZ calculations with a KI starting point
                if functional == 'ki':
                    # pKIPZ
                    utils.system_call('rsync -a ki/ pkipz/')

                    # KIPZ
                    utils.system_call('rsync -a ki/final/ kipz/init/')
                    utils.system_call('mv kipz/init/ki_final.cpi kipz/init/ki_init.cpi')
                    utils.system_call('mv kipz/init/ki_final.cpo kipz/init/ki_init.cpo')
                    if self.parameters.init_orbitals in ['mlwfs', 'projwfs']:
                        utils.system_call('rsync -a ki/init/wannier kipz/init/')
                    if self.parameters.periodic and self.master_calc_params['ui'].do_smooth_interpolation:
                        # Copy over the smooth PBE calculation from KI for KIPZ to use
                        utils.system_call('rsync -a ki/postproc kipz/')
                        utils.system_call('find kipz/postproc/ -name "*interpolated.dat" -delete')

        else:
            # self.functional != all and self.method != 'dfpt'
            if self.parameters.functional in ['ki', 'pkipz', 'kipz']:
                dscf_workflow = KoopmansDSCFWorkflow(**self.wf_kwargs)
                self.run_subworkflow(dscf_workflow)
            else:
                dft_workflow = DFTCPWorkflow(**self.wf_kwargs)
                self.run_subworkflow(dft_workflow)
