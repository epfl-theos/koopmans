"""

singlepoint workflow object for koopmans

Written by Edward Linscott Oct 2020
Converted to a workflow object Nov 2020

"""

import os
import shutil
from pathlib import Path
from typing import Generator, List

import numpy as np

from koopmans import utils
from koopmans.calculators import ProjwfcCalculator
from koopmans.files import File
from koopmans.outputs import OutputModel
from koopmans.status import Status

from ._dft import DFTCPWorkflow, DFTPhWorkflow
from ._koopmans_dfpt import KoopmansDFPTWorkflow
from ._koopmans_dscf import KoopmansDSCFWorkflow
from ._workflow import Workflow

load_results_from_output = True


class SinglepointOutputs(OutputModel):
    '''
    Outputs for the SinglepointWorkflow
    '''
    pass


class SinglepointWorkflow(Workflow):

    '''
    Examples
    --------

    Running a Koopmans calculation on ozone

        >>> from ase_koopmans.build import molecule
        >>> from koopmans.workflows import SinglepointWorkflow
        >>> ozone = molecule('O3', vacuum=5.0, pbc=False)
        >>> wf = SinglepointWorkflow(ozone, ecutwfc = 20.0)
        >>> wf.run()

    Running a Koopmans calculation on GaAs

        >>> from ase_koopmans.build import bulk
        >>> from koopmans.projections import ProjectionBlocks
        >>> from koopmans.kpoints import Kpoints
        >>> from koopmans.workflows import SinglepointWorkflow
        >>> gaas = bulk('GaAs', crystalstructure='zincblende', a=5.6536)
        >>> projs = ProjectionBlocks.fromlist([["Ga: d"], ["As: sp3"], ["Ga: sp3"]],
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

    output_model = SinglepointOutputs  # type: ignore

    def _run(self) -> None:

        # Import it like this so if they have been monkey-patched, we will get the monkey-patched version
        if self.parameters.eps_inf == 'auto':
            eps_workflow = DFTPhWorkflow.fromparent(self)
            eps_workflow.run()
            if eps_workflow.status != Status.COMPLETED:
                return
            self.parameters.eps_inf = np.trace(eps_workflow.calculations[-1].results['dielectric tensor']) / 3

        if self.parameters.method == 'dfpt':
            workflow = KoopmansDFPTWorkflow.fromparent(self)
            workflow.run()
            if workflow.status != Status.COMPLETED:
                return

        elif self.parameters.functional == 'all':
            # if 'all', create subdirectories and run
            functionals = ['ki', 'pkipz', 'kipz']

            if self.parameters.alpha_from_file:
                utils.warn('Need to make sure alpharef files are copied over')
                # utils.system_call('cp file_alpharef*.txt ki/')

            ki_workflow = None
            for functional in functionals:
                # For pKIPZ/KIPZ, use KI as a starting point
                restart_from_old_ki = (functional == 'kipz')

                # We only need to do the smooth interpolation the first time (i.e. for KI)
                redo_smooth_dft = None if functional == 'ki' else False

                # For pKIPZ should not be recalculated
                if self.parameters.calculate_alpha:
                    calculate_alpha = (functional != 'pkipz')
                else:
                    calculate_alpha = self.parameters.calculate_alpha

                # Create a KC workflow for this particular functional
                if functional != 'ki':
                    assert ki_workflow is not None
                    variational_orbital_files = ki_workflow.outputs.variational_orbital_files
                    previous_cp_calc = ki_workflow.outputs.final_calc
                    smooth_dft_ham_files = ki_workflow.outputs.smooth_dft_ham_files
                else:
                    variational_orbital_files = None
                    previous_cp_calc = None
                    smooth_dft_ham_files = None

                kc_workflow = KoopmansDSCFWorkflow.fromparent(self, functional=functional,
                                                              initial_variational_orbital_files=variational_orbital_files,
                                                              redo_smooth_dft=redo_smooth_dft,
                                                              smooth_dft_ham_files=smooth_dft_ham_files,
                                                              previous_cp_calc=previous_cp_calc,
                                                              calculate_alpha=calculate_alpha)
                kc_workflow.name += ' ' + functional.upper().replace("PKIPZ", "pKIPZ")

                if functional == 'ki':
                    # Save the KI workflow for later
                    ki_workflow = kc_workflow

                # Transform to the supercell
                if functional == 'kipz':
                    kc_workflow.primitive_to_supercell()

                # Run the workflow
                kc_workflow.run()
                if kc_workflow.status != Status.COMPLETED:
                    return

        else:
            # self.functional != all and self.method != 'dfpt'
            if self.parameters.functional in ['ki', 'pkipz', 'kipz']:
                dscf_workflow = KoopmansDSCFWorkflow.fromparent(self)
                dscf_workflow.run()
                if dscf_workflow.status != Status.COMPLETED:
                    return
            else:
                dft_workflow = DFTCPWorkflow.fromparent(self)
                dft_workflow.run()
                if dft_workflow.status != Status.COMPLETED:
                    return

        self.status = Status.COMPLETED

        return
