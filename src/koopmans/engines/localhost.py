import contextlib
import os
import shutil
import sys
from itertools import chain
from pathlib import Path
from typing import Generator, List, Literal, overload

from ase_koopmans.calculators.calculator import CalculationFailed
from upf_tools import UPFDict

from koopmans import utils
from koopmans.calculators import (Calc, ImplementedCalc, PhCalculator,
                                  ProjwfcCalculator, ReturnsBandStructure)
from koopmans.files import File, LocalFile
from koopmans.processes import Process
from koopmans.pseudopotentials import (element_from_pseudo_filename,
                                       local_libraries)
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

        assert step.directory is not None
        if step.directory.exists():
            assert isinstance(step, utils.HasDirectory)
            self.rmdir(File(step, ''))

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
                assert step.directory is not None
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
        pseudo_dir = LocalFile(Path(__file__).parents[1] / 'pseudopotentials' / library)
        pseudo_name = None
        for pseudo_file in chain(self.glob(pseudo_dir, '*.upf'), self.glob(pseudo_dir, '*.UPF')):
            if element_from_pseudo_filename(pseudo_file.name.name) == element:
                pseudo_name = pseudo_file.name
                break

        if pseudo_name is None:
            raise ValueError(f'Pseudo for {element} not found in library {library}')

        pseudo_path = pseudo_dir / pseudo_name
        return UPFDict.from_upf(pseudo_path.aspath())

    @overload
    def read_file(self, file: File, binary: Literal[True]) -> bytes: ...

    @overload
    def read_file(self, file: File, binary: Literal[False]) -> str: ...

    @overload
    def read_file(self, file: File, binary: bool = False) -> bytes | str: ...

    def read_file(self, file: File, binary: bool = False) -> bytes | str:
        assert file.parent.absolute_directory is not None
        full_path = file.parent.absolute_directory / file.name
        fstring = 'rb' if binary else 'r'
        with open(full_path, fstring) as f:
            flines = f.read()
        return flines

    def write_file(self, content: str | bytes, file: File) -> None:
        fstring = 'wb' if isinstance(content, bytes) else 'w'
        if len(file.name.parents) > 0:
            parent_directory = File(file.parent, file.name.parent)
            self.mkdir(parent_directory, parents=True, exist_ok=True)
        with open(file.aspath(), fstring) as f:
            f.write(content)

    def copy_file(self, source: File, destination: File, exist_ok: bool = False) -> None:
        utils.copy_file(source.aspath(), destination.aspath(), exist_ok=exist_ok)

    def link_file(self, source: File, destination: File, recursive: bool = False, overwrite: bool = False) -> None:
        if recursive:
            utils.symlink_tree(source.aspath(), destination.aspath())
        else:
            utils.symlink(source.aspath(), destination.aspath(), exist_ok=overwrite, force=overwrite)

    def file_exists(self, file: File) -> bool:
        return file.aspath().exists()

    def file_is_dir(self, file: File) -> bool:
        return file.aspath().is_dir()

    @contextlib.contextmanager
    def chdir(self, directory: Path):
        return utils.chdir_logic(directory)

    def rmdir(self, directory: File) -> None:
        path = directory.aspath()
        if isinstance(path, str):
            path = Path(path)
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()

    def mkdir(self, directory: File, parents: bool = False, exist_ok: bool = False) -> None:
        directory.aspath().mkdir(parents=parents, exist_ok=exist_ok)

    def glob(self, directory: File, pattern: str, recursive: bool = False) -> Generator[File, None, None]:
        assert directory.parent is not None
        assert directory.parent.directory is not None
        if recursive:
            generator = directory.aspath().rglob(pattern)
        else:
            generator = directory.aspath().glob(pattern)
        for path in generator:
            yield File(parent=directory.parent, name=path.relative_to(directory.parent.directory))

    def available_pseudo_families(self) -> set[str]:
        return local_libraries


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
