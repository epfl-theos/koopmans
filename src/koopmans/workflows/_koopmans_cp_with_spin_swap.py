"""A workflow for managing kcp.x calculations that would want more spin-down electrons than spin-up"""

from typing import Generator

from pydantic import ConfigDict

from koopmans.calculators import KoopmansCPCalculator
from koopmans.files import File
from koopmans.process_io import IOModel
from koopmans.processes.koopmans_cp import SwapSpinFilesProcess
from koopmans.status import Status

from ._workflow import Workflow


class KoopmansCPWithSpinSwapOutput(IOModel):
    outdir: File
    model_config = ConfigDict(arbitrary_types_allowed=True)


class KoopmansCPWithSpinSwapWorkflow(Workflow[KoopmansCPWithSpinSwapOutput]):

    output_model = KoopmansCPWithSpinSwapOutput

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
        process = SwapSpinFilesProcess(read_directory=self._calc.write_directory)
        status = self.run_steps(process)
        if status != Status.COMPLETED:
            return

        self.outputs = self.output_model(outdir=process.outputs.write_directory)

        self.status = Status.COMPLETED
