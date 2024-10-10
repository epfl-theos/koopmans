"""A workflow for managing kcp.x calculations that would want more spin-down electrons than spin-up"""

from koopmans.calculators import KoopmansCPCalculator
from koopmans.files import FilePointer
from koopmans.outputs import OutputModel
from koopmans.processes.koopmans_cp import SwapSpinFilesProcess

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

    def _run(self):
        assert self._calc.parameters.nspin == 2
        assert self._calc.parameters.nelup < self._calc.parameters.neldw

        # Run the calculation (self._calc.swap_spin_files() will take care of swapping the input parameters,
        # the linked files, and the python output data)
        self.run_calculator(self._calc)

        # Swap the spin-up and spin-down output files using a process
        process = SwapSpinFilesProcess(read_directory=FilePointer(self._calc, self._calc.write_directory))
        self.run_process(process)

        self.outputs = self.output_model(outdir=FilePointer(process, process.outputs.write_directory))
