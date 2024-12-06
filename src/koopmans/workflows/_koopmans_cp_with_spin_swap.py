"""A workflow for managing kcp.x calculations that would want more spin-down electrons than spin-up"""

from typing import Generator

from koopmans.calculators import KoopmansCPCalculator
from koopmans.files import FilePointer
from koopmans.outputs import OutputModel
from koopmans.processes.koopmans_cp import SwapSpinFilesProcess
from koopmans.status import Status
from koopmans.step import Step

from ._workflow import Workflow


class KoopmansCPWithSpinSwapOutput(OutputModel):
    outdir: FilePointer

    class Config:
        arbitrary_types_allowed = True


class KoopmansCPWithSpinSwapWorkflow(Workflow):

    output_model = KoopmansCPWithSpinSwapOutput
    outputs: KoopmansCPWithSpinSwapOutput

    def __init__(self, calc: KoopmansCPCalculator, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._calc = calc

    def _run(self) -> None:
        assert self._calc.parameters.nspin == 2
        assert self._calc.parameters.nelup < self._calc.parameters.neldw

        # Run the calculation (self._calc.swap_spin_files() will take care of swapping the input parameters,
        # the linked files, and the python output data)
        status = self.run_steps(self._calc)
        if status != Status.COMPLETED:
            return

        # Swap the spin-up and spin-down output files using a process
        process = SwapSpinFilesProcess(read_directory=FilePointer(self._calc, self._calc.write_directory))
        status = self.run_steps(process)
        if status != Status.COMPLETED:
            return

        self.outputs = self.output_model(outdir=FilePointer(process, process.outputs.write_directory))

        self.status = Status.COMPLETED
