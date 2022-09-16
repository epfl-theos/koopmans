"""

singlepoint workflow object for koopmans

Written by Edward Linscott Oct 2020
Converted to a workflow object Nov 2020

"""

import os
import shutil
from pathlib import Path

import numpy as np

from koopmans import utils

from ._workflow import Workflow

load_results_from_output = True


class SinglepointWorkflow(Workflow):

    '''
    Examples
    --------

    Running a Koopmans calculation on ozone

        >>> from ase.build import molecule
        >>> from koopmans.workflows import SinglepointWorkflow
        >>> ozone = molecule('O3', vacuum=5.0, pbc=False)
        >>> wf = SinglepointWorkflow(ozone, ecutwfc = 20.0)
        >>> wf.run()

    Running a Koopmans calculation on GaAs

        >>> from ase.build import bulk
        >>> from koopmans.projections import ProjectionBlocks
        >>> from koopmans.kpoints import Kpoints
        >>> from koopmans.workflows import SinglepointWorkflow
        >>> gaas = bulk('GaAs', crystalstructure='zincblende', a=5.6536)
        >>> projs = ProjectionBlocks.fromprojections([["Ga: d"], ["As: sp3"], ["Ga: sp3"]],
        >>>                                           fillings=[True, True, False],
        >>>                                           spins=[None, None, None],
        >>>                                           atoms=gaas)
        >>> kpoints = Kpoints(grid=[2, 2, 2])
        >>> wf = SinglepointWorkflow(gaas, kpoints=kpoints, projections=projs, init_orbitals='mlwfs',
        >>>                          pseudo_library='sg15_v1.0', ecutwfc=40.0,
        >>>                          calculator_parameters={'pw': {'nbnd': 45},
        >>>                          'w90_emp': {'dis_froz_max': 14.6, 'dis_win_max': 18.6}})
        >>> wf.run()
    '''

    if Workflow.__doc__ and __doc__:
        __doc__ = Workflow.__doc__ + __doc__

    def _run(self) -> None:

        # Import it like this so if they have been monkey-patched, we will get the monkey-patched version
        from koopmans.workflows import (DFTCPWorkflow, DFTPhWorkflow,
                                        KoopmansDFPTWorkflow,
                                        KoopmansDSCFWorkflow)

        if self.parameters.eps_inf == 'auto':
            eps_workflow = DFTPhWorkflow.fromparent(self)
            if self.parameters.from_scratch and Path('calculate_eps').exists():
                shutil.rmtree('calculate_eps')
            eps_workflow.run(subdirectory='calculate_eps')
            self.parameters.eps_inf = np.trace(eps_workflow.calculations[-1].results['dielectric tensor']) / 3

        if self.parameters.method == 'dfpt':
            workflow = KoopmansDFPTWorkflow.fromparent(self)
            workflow.run()

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

                # For pKIPZ/KIPZ, use KI as a starting point
                restart_from_old_ki = (functional == 'kipz')

                # We only need to do the smooth interpolation the first time (i.e. for KI)
                redo_smooth_dft = None if functional == 'ki' else False

                # Create a KC workflow for this particular functional
                kc_workflow = KoopmansDSCFWorkflow.fromparent(self, functional=functional,
                                                              restart_from_old_ki=restart_from_old_ki,
                                                              redo_smooth_dft=redo_smooth_dft,
                                                              calculate_alpha=(functional != 'pkipz'))

                # Transform to the supercell
                if functional == 'kipz':
                    kc_workflow.primitive_to_supercell()

                # Run the workflow
                if functional == 'pkipz' and self.parameters.from_scratch:
                    # We want to run pKIPZ with from_scratch = False, but don't want this to be inherited
                    kc_workflow.run(subdirectory=functional, from_scratch=False)
                else:
                    kc_workflow.run(subdirectory=functional)

                # Provide the pKIPZ and KIPZ calculations with a KI starting point
                if functional == 'ki':
                    if self.parameters.from_scratch:
                        for directory in ['pkipz', 'kipz']:
                            if os.listdir(directory):
                                utils.system_call(f'rm -rf {directory}')

                    # pKIPZ
                    for dir in ['init', 'calc_alpha', 'TMP-CP']:
                        src = Path(f'ki/{dir}/')
                        if src.is_dir():
                            utils.system_call(f'rsync -a {src} pkipz/')
                    for f in ['ki_final.cpi', 'ki_final.cpo', 'ki_final.ham.pkl', 'ki_final.bare_ham.pkl']:
                        file = Path(f'ki/final/{f}')
                        if file.is_file():
                            utils.system_call(f'rsync -a {file} pkipz/final/')
                    if all(self.atoms.pbc) and self.calculator_parameters['ui'].do_smooth_interpolation:
                        if not Path('pkipz/postproc').is_dir():
                            utils.system_call('mkdir pkipz/postproc')
                        for dir in ['wannier', 'TMP', 'pdos']:
                            utils.system_call(f'rsync -a ki/postproc/{dir} pkipz/postproc/')

                    # KIPZ
                    utils.system_call('rsync -a ki/final/ kipz/init/')
                    utils.system_call('mv kipz/init/ki_final.cpi kipz/init/ki_init.cpi')
                    utils.system_call('mv kipz/init/ki_final.cpo kipz/init/ki_init.cpo')
                    if self.parameters.init_orbitals in ['mlwfs', 'projwfs']:
                        utils.system_call('rsync -a ki/init/wannier kipz/init/')
                    if all(self.atoms.pbc) and self.calculator_parameters['ui'].do_smooth_interpolation:
                        # Copy over the smooth PBE calculation from KI for KIPZ to use
                        for dir in ['wannier', 'TMP', 'pdos']:
                            utils.system_call(f'rsync -a ki/postproc/{dir} kipz/postproc/')

        else:
            # self.functional != all and self.method != 'dfpt'
            if self.parameters.functional in ['ki', 'pkipz', 'kipz']:
                dscf_workflow = KoopmansDSCFWorkflow.fromparent(self)
                dscf_workflow.run()
            else:
                dft_workflow = DFTCPWorkflow.fromparent(self)
                dft_workflow.run()
