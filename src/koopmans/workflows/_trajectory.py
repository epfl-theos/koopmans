"""

A workflow for serially running the Koopmans DSCF workflow on multiple atomic configurations

"""

import shutil
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional

import numpy as np
from ase_koopmans import Atoms, io
from sklearn.metrics import mean_absolute_error, r2_score

from koopmans import calculators, utils
from koopmans.outputs import OutputModel
from koopmans.status import Status
from koopmans.step import Step

from ._koopmans_dscf import KoopmansDSCFOutputs
from ._workflow import Workflow


class TrajectoryOutputs(OutputModel):
    snapshot_outputs: List[KoopmansDSCFOutputs]


class TrajectoryWorkflow(Workflow):

    output_model = TrajectoryOutputs  # type: ignore
    outputs: TrajectoryOutputs

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.parameters.task = 'trajectory'

    def _run(self) -> None:
        """
        Starts the KoopmansDSCF Workflow for each snapshot indicated in indices
        """

        # Import it like this so if they have been monkey-patched, we will get the monkey-patched version
        from koopmans.workflows import KoopmansDSCFWorkflow

        workflows = []
        for i, snapshot in enumerate(self.snapshots):
            # Get the atomic positions for the current snapshot
            self.atoms.set_positions(snapshot.positions)

            # after each snapshot we want to set the from_scratch_parameter to its original value
            # To do so, we save it here since it might get set from False to True during the calculation
            # of the snapshot
            from_scratch = self.parameters.from_scratch

            # Initialize and run the DSCF workflow
            workflow = KoopmansDSCFWorkflow.fromparent(self)
            workflow.name += f' Snapshot {i+1} of {len(self.snapshots)}'
            workflows.append(workflow)

        for w in workflows:
            w.run()

        if not all([w.status == Status.COMPLETED for w in workflows]):
            return

        self.outputs = self.output_model(snapshot_outputs=[w.outputs for w in workflows])
        self.status = Status.COMPLETED

        return
