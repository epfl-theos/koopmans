"""A workflow for running the Koopmans DSCF workflow on multiple atomic configurations."""

from typing import List

from koopmans.pydantic import BaseModel
from koopmans.status import Status

from ._koopmans_dscf import KoopmansDSCFOutputs, KoopmansDSCFWorkflow
from ._workflow import Workflow, WorkflowInput


class TrajectoryOutputs(BaseModel):
    """Pydantic model for the outputs of a `TrajectoryWorkflow`."""

    snapshot_outputs: List[KoopmansDSCFOutputs]


class TrajectoryWorkflow(Workflow[WorkflowInput, TrajectoryOutputs]):
    """A workflow for running the KoopmansDSCF Workflow on multiple snapshots."""

    input_model = WorkflowInput
    output_model = TrajectoryOutputs

    def _run(self) -> None:
        """Run a KoopmansDSCF Workflow for each snapshot indicated in indices."""
        workflows = []
        for i, snapshot in enumerate(self.snapshots):
            # Get the atomic positions for the current snapshot
            self.atoms.set_positions(snapshot.positions)

            # Initialize and run the DSCF workflow
            workflow = KoopmansDSCFWorkflow.fromparent(self)
            workflow.name += f' Snapshot {i + 1} of {len(self.snapshots)}'
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
