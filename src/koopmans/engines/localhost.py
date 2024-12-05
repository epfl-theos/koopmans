import os
import sys
from typing import Generator, List

from ase_koopmans.calculators.calculator import CalculationFailed
from upf_tools import UPFDict

from koopmans import utils
from koopmans.calculators import (Calc, ImplementedCalc, PhCalculator,
                                  ProjwfcCalculator, ReturnsBandStructure)
from koopmans.files import FilePointer
from koopmans.processes import Process
from koopmans.pseudopotentials import (fetch_pseudo, pseudos_library_directory,
                                       read_pseudo_file)
from koopmans.status import Status
from koopmans.step import Step

from .engine import Engine


class LocalhostEngine(Engine):
    # Engine for running a workflow locally
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.statuses = {}

    def run(self, step: Step):
        self._step_running_message(step)

        try:
            step.run()
        except CalculationFailed:
            self.set_status(step, Status.FAILED)
            self._step_failed_message(step)
            raise

        self.set_status(step, Status.COMPLETED)
        self._step_completed_message(step)

        # If we reached here, all future steps should be performed from scratch
        self.from_scratch = True

        return

    def load_old_calculator(self, calc: Calc):
        return load_old_calculator(calc)

    def load_results(self, step: Step):
        # For the local calculation, step.run() also loads the results of the calculator
        pass

    def get_status(self, step: Step) -> Status:

        if not self.from_scratch:
            to_run = True
            if isinstance(step, ImplementedCalc):
                calc_file = step.directory / step.prefix

                if calc_file.with_suffix(step.ext_out).is_file():
                    loaded_step = self.load_old_calculator(step)

                    if loaded_step.is_complete():
                        to_run = False
            elif isinstance(step, Process):
                if step.is_complete():
                    to_run = False
                    step.load_outputs()
            else:
                raise ValueError(f'Unknown step type: {type(step)}')

            if to_run:
                # This and subsequent calculations should be performed from scratch
                self.from_scratch = True
            else:
                if step.uid not in self.statuses:
                    self._step_skipped_message(step)
                self.set_status(step, Status.COMPLETED)

        if step.uid not in self.statuses:
            return Status.NOT_STARTED

        return self.statuses[step.uid]

    def set_status(self, step: Step, status: Status):
        self.statuses[step.uid] = status

    def update_statuses(self) -> None:
        pass

    def get_pseudopotential(self, library: str, element: str) -> UPFDict:
        pseudo = fetch_pseudo(library=library, element=element)
        pseudo_path = pseudos_library_directory(pseudo.library) / pseudo.name
        return read_pseudo_file(pseudo_path)

    def read(self, file: FilePointer, binary=False) -> str | bytes:
        return utils.get_content(*file)

    def write(self, content: str | bytes, file: FilePointer) -> None:
        if isinstance(content, bytes):
            utils.write_binary_content(file.aspath(), content)
        else:
            utils.write_content(file.aspath(), content)

    def glob(self, pattern: FilePointer, recursive=False) -> Generator[FilePointer, None, None]:
        raise NotImplementedError()
        # Need to use FilePointer.glob() and FilePointer.rglob() but first splitting pattern into
        # its parent and the pattern itself


def load_old_calculator(calc):
    # This is a separate function so that it can be imported by other engines
    loaded_calc = calc.__class__.fromfile(calc.directory / calc.prefix)

    if loaded_calc.is_complete():
        # If it is complete, load the results
        calc.results = loaded_calc.results

        # Check the convergence of the calculation
        calc.check_convergence()

        # Load k-points if relevant
        if hasattr(loaded_calc, 'kpts'):
            calc.kpts = loaded_calc.kpts

        if isinstance(calc, ReturnsBandStructure):
            calc.generate_band_structure()

        if isinstance(calc, ProjwfcCalculator):
            calc.generate_dos()

        if isinstance(calc, PhCalculator):
            calc.read_dynG()

    return loaded_calc
