"""

A workflow for serially running a workflow on multiple atomic configurations

"""

import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Type, TypeVar

from ase import Atoms, io

from koopmans import utils

from ._workflow import Workflow
from ._wrapper import WorkflowWrapper, WorkflowWrapperFactory

W = TypeVar('W', bound='Workflow')


class SnapshotWorkflow(WorkflowWrapper):

    _wrapper_args: List[str] = ['get_eigenvalues']
    snapshots: List[Atoms]

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
        wf = super(SnapshotWorkflow, cls)._fromjsondct(bigdct, override)
        wf.snapshots = snapshots
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
        Starts the subworkflow for each snapshot
        """

        for i, snapshot in enumerate(self.snapshots):

            self.print(
                f'Performing calculation on snapshot {i+1} / {len(self.snapshots)}', style='heading')

            # Get the atomic positions for the current snapshot
            subdirectory = f'snapshot_{i+1}'
            self.ml.current_snapshot = i
            self.atoms.set_positions(snapshot.positions)

            # after each snapshot we want to set the from_scratch_parameter to its original value
            # To do so, we save it here since it might get set from False to True during the calculation
            # of the snapshot
            from_scratch = self.parameters.from_scratch

            # If we are interested in the prediction of the eigenvalues, delete the final directory to make sure that
            # the final calculation is rerun.
            if self.get_eigenvalues:
                shutil.rmtree(Path(subdirectory) / 'final', ignore_errors=True)

            # Initialize and run the subworkflow
            subworkflow = self.subworkflow_class.fromparent(self)
            subworkflow.run(subdirectory=subdirectory)

            # Reset the from_scratch parameter to its original value
            self.parameters.from_scratch = from_scratch

            # Wipe the history of the bands
            delattr(self, '_bands')


def SnapshotWorkflowFactory(subworkflow: Workflow, snapshots: List[Atoms]) -> SnapshotWorkflow:

    wf = WorkflowWrapperFactory(SnapshotWorkflow, subworkflow)

    wf.snapshots = snapshots

    return wf
