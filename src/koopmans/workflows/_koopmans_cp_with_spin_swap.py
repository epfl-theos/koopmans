"""A workflow for managing kcp.x calculations that would want more spin-down electrons than spin-up."""

from koopmans.calculators.kcp import KCPCalculator
from koopmans.files import File
from koopmans.processes.koopmans_cp import SwapSpinFilesProcess
from koopmans.pydantic import BaseModel
from koopmans.status import Status

from ._workflow import Workflow, WorkflowInput


class KoopmansCPWithSpinSwapInputs(WorkflowInput):
    """Input model for the KoopmansCPWithSpinSwapWorkflow."""
    calc: KCPCalculator


class KoopmansCPWithSpinSwapOutput(BaseModel):
    """Output model for the KoopmansCPWithSpinSwapWorkflow."""

    outdir: File


class KoopmansCPWithSpinSwapWorkflow(Workflow[KoopmansCPWithSpinSwapInputs, KoopmansCPWithSpinSwapOutput]):
    """Workflow for running a kcp.x calculation that would want more spin-down electrons than spin-up."""

    input_model = KoopmansCPWithSpinSwapInputs
    output_model = KoopmansCPWithSpinSwapOutput

    def _run(self) -> None:
        """Run the workflow."""
        assert self._calc.parameters.nspin == 2
        assert self._calc.parameters.nelup < self._calc.parameters.neldw

        # Run the calculation (self._calc.swap_spin_files() will take care of swapping the input parameters,
        # the linked files, and the python output data)
        status = self.run_steps(self._calc)
        if status != Status.COMPLETED:
            return

        # Swap the spin-up and spin-down output files using a process
        process = SwapSpinFilesProcess(read_directory=self._calc.write_directory)
        status = self.run_steps(process)
        if status != Status.COMPLETED:
            return

        self.outputs = self.output_model(outdir=process.outputs.write_directory)

        self.status = Status.COMPLETED
