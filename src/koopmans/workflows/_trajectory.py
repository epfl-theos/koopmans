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
from koopmans.process_io import IOModel
from koopmans.status import Status

from ._koopmans_dscf import KoopmansDSCFOutputs, KoopmansDSCFWorkflow
from ._workflow import Workflow


class TrajectoryOutputs(IOModel):
    snapshot_outputs: List[KoopmansDSCFOutputs]


class TrajectoryWorkflow(Workflow[TrajectoryOutputs]):

    output_model = TrajectoryOutputs

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.parameters.task = 'trajectory'

    def _run(self) -> None:
        """
        Starts the KoopmansDSCF Workflow for each snapshot indicated in indices
        """

        # Import it like this so if they have been monkey-patched, we will get the monkey-patched version

        workflows = []
        for i, snapshot in enumerate(self.snapshots):
            # Get the atomic positions for the current snapshot
            self.atoms.set_positions(snapshot.positions)

            # Initialize and run the DSCF workflow
            workflow = KoopmansDSCFWorkflow.fromparent(self)
            workflow.name += f' Snapshot {i+1} of {len(self.snapshots)}'
            workflows.append(workflow)

        for w in workflows:
            w.proceed(copy_outputs_to_parent=False)

        if not all([w.status == Status.COMPLETED for w in workflows]):
            return

        self.calculations += [c for w in workflows for c in w.calculations]
        self.steps += [s for w in workflows for s in w.steps]

        self.outputs = self.output_model(snapshot_outputs=[w.outputs for w in workflows])
        self.status = Status.COMPLETED

        return
