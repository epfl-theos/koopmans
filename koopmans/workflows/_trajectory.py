

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from sklearn.metrics import mean_absolute_error, r2_score

from ase import Atoms, io
from koopmans import calculators, utils

from ._workflow import Workflow


class TrajectoryWorkflow(Workflow):

    def __init__(self, indices=None, save_dir=None, get_evs=False, *args, **kwargs):
        snapshots: List[Atoms] = kwargs.pop('snapshots', [])
        if 'atoms' not in kwargs and snapshots != []:
            kwargs['atoms'] = snapshots[0]

        super().__init__(*args, **kwargs)
        self.snapshots = snapshots
        self.number_of_snapshots = len(self.snapshots)

        self.indices: Optional[List[int]] = indices
        self.save_dir: Optional[Path] = save_dir
        self.get_evs: Optional[bool] = get_evs

    @ classmethod
    def _fromjsondct(cls, bigdct: Dict[str, Any]):
        """
        Reads the atomic positions for each snapshot from the xyz-file specified by the user in the snapshots-file.
        """

        try:
            snapshots_file = bigdct['setup'].pop('snapshots')
        except:
            raise ValueError(
                f'To calculate a trajectory, please provide a xyz-file containing the atomic positions of the snapshots in the setup-block of the json-input file.')

        snapshots = io.read(snapshots_file, index=':')
        if isinstance(snapshots, Atoms):
            snapshots = [snapshots]
        bigdct['setup']['atomic_positions'] = utils.construct_atomic_positions_block(snapshots[0])
        wf = super(TrajectoryWorkflow, cls)._fromjsondct(bigdct)
        wf.snapshots = snapshots
        wf.number_of_snapshots = len(snapshots)
        return wf

    def todict(self):

        dct = dict(self.__dict__)

        items_to_pop = ['atoms']
        for item in items_to_pop:
            dct.pop(item)

        return dct

    def _run(self):
        """
        Starts the KoopmansDSCF Workflow for each snapshot indicated in indidses
        """

        # Import it like this so if they have been monkey-patched, we will get the monkey-patched version
        from koopmans.workflows import KoopmansDSCFWorkflow

        if self.indices is None:
            self.indices = list(range(0, self.number_of_snapshots))

        for i in self.indices:
            self.print(
                f'Performing Koopmans calculation on snapshot {i+1} / {self.number_of_snapshots}', style='heading')

            # get the atomic positions for the current snapshot
            subdirectory = f'snapshot_{i+1}'
            snapshot = self.snapshots[i]
            self.parameters.current_snapshot = i
            self.atoms.set_positions(snapshot.positions)

            # if we are interested in the prediction of the eigenvalues,
            # delete the final directory to make sure that the final calcualtion is rerun.
            if self.get_evs:
                from_scratch = self.parameters.from_scratch
                utils.system_call(f'rm -rf {subdirectory}/final')

            # initialize and run the DSCF workflow
            # workflow = KoopmansDSCFWorkflow(**self.wf_kwargs)
            workflow = KoopmansDSCFWorkflow.fromparent(self)
            self.bands = workflow.bands  # reset the bands to the initial guesses
            workflow.run(subdirectory=subdirectory)

            # since we have deleted the final directory and therefore had to rerun, we must now make sure, that from_scratch is again set to its original value
            if self.get_evs:
                self.parameters.from_scratch = from_scratch

            # if necessary, save the results (e.g. for the convergence analysis)
            if self.save_dir is not None:
                alphas = self.bands.alphas
                np.savetxt(self.save_dir / f"alphas_snapshot_{i+1}.txt", alphas)
                if self.get_evs:
                    final_calculator = calculators.KoopmansCPCalculator
                    final_calc = [c for c in workflow.calculations if isinstance(c, final_calculator)][-1]
                    evs = final_calc.results['eigenvalues']
                    np.savetxt(self.save_dir / f"evs_snapshot_{i+1}.txt", evs)
