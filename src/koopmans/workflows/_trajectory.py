"""

A workflow for serially running the Koopmans DSCF workflow on multiple atomic configurations

"""

import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from ase import Atoms, io
from sklearn.metrics import mean_absolute_error, r2_score

from koopmans import calculators, utils

from ._workflow import Workflow


class TrajectoryWorkflow(Workflow):

    def __init__(self, snapshots: List[Atoms], indices: Optional[List[int]] = None, save_dir: Optional[Path] = None,
                 get_evs: bool = False, overwrite_atoms: bool = True, *args, **kwargs):

        if overwrite_atoms:
            if isinstance(snapshots, Atoms):
                kwargs['atoms'] = snapshots
            else:
                kwargs['atoms'] = snapshots[0]
        super().__init__(*args, **kwargs)

        self.snapshots = snapshots
        self.number_of_snapshots = len(self.snapshots)
        self.indices: Optional[List[int]] = indices
        self.save_dir: Optional[Path] = save_dir
        self.get_evs: Optional[bool] = get_evs
        all_alphas: Dict[str, np.ndarray] = {}
        self.all_alphas = all_alphas

    @ classmethod
    def _fromjsondct(cls, bigdct: Dict[str, Any], override: Dict[str, Any] = {}):
        """
        Reads the atomic positions for each snapshot from the xyz file provided by the user
        """

        try:
            snapshots_file = bigdct['atoms']['atomic_positions'].pop('snapshots')
        except:
            raise ValueError('To calculate a trajectory, please provide a "snapshots" entry in the'
                             '"atomic_positions" block, corresponding to the name of an xyz-formatted '
                             'file containing the snapshots')

        snapshots = io.read(snapshots_file, index=':')
        if isinstance(snapshots, Atoms):
            snapshots = [snapshots]
        bigdct['atoms']['atomic_positions'] = utils.construct_atomic_positions_block(snapshots[0])
        wf = super(TrajectoryWorkflow, cls)._fromjsondct(bigdct, override)
        wf.snapshots = snapshots
        wf.number_of_snapshots = len(snapshots)
        return wf

    def toinputjson(self) -> Dict[str, Dict[str, Any]]:
        bigdct = super().toinputjson()
        snapshots_file = "snapshots.json"
        io.write(snapshots_file, self.snapshots)
        bigdct['atoms']['atomic_positions'] = {"snapshots": snapshots_file}
        return bigdct

    def todict(self):

        dct = dict(self.__dict__)

        items_to_pop = ['atoms']
        for item in items_to_pop:
            dct.pop(item)

        return dct

    def _run(self):
        """
        Starts the KoopmansDSCF Workflow for each snapshot indicated in indices
        """

        # Import it like this so if they have been monkey-patched, we will get the monkey-patched version
        from koopmans.workflows import KoopmansDSCFWorkflow

        if self.indices is None:
            self.indices = list(range(0, self.number_of_snapshots))

        for i in self.indices:

            self.print(
                f'Performing Koopmans calculation on snapshot {i+1} / {self.number_of_snapshots}', style='heading')

            # Get the atomic positions for the current snapshot
            subdirectory = f'snapshot_{i+1}'
            snapshot = self.snapshots[i]
            self.ml.current_snapshot = i
            self.atoms.set_positions(snapshot.positions)

            # If we are interested in the prediction of the eigenvalues, delete the final directory to make sure that
            # the final calculation is rerun.
            if self.get_evs:
                from_scratch = self.parameters.from_scratch
                shutil.rmtree(Path(subdirectory) / 'final', ignore_errors=True)

            # Initialize and run the DSCF workflow
            workflow = KoopmansDSCFWorkflow.fromparent(self)
            self.bands = workflow.bands  # reset the bands to the initial guesses
            workflow.run(subdirectory=subdirectory)

            # Since we have deleted the final directory and therefore had to rerun, we must now make sure that
            # from_scratch is set to its original value
            if self.get_evs:
                self.parameters.from_scratch = from_scratch

            # If necessary, save the results (e.g. for the convergence analysis)
            alphas = self.bands.alphas
            self.all_alphas[f'snapshot_{i+1}'] = alphas
            if self.save_dir is not None:
                np.savetxt(self.save_dir / f"alphas_snapshot_{i+1}.txt", alphas)
                if self.get_evs:
                    final_calculator = calculators.KoopmansCPCalculator
                    final_calc = [c for c in workflow.calculations if isinstance(c, final_calculator)][-1]
                    evs = final_calc.results['eigenvalues']
                    np.savetxt(self.save_dir / f"evs_snapshot_{i+1}.txt", evs)
