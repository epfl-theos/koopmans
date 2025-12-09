from pathlib import Path
from typing import ClassVar, Generic, TypeVar

from ase_koopmans import Atoms
from pydantic import PrivateAttr

from koopmans.processes.process import Process
from koopmans.pydantic import BaseModel


class CalculatorInputs(BaseModel):
    atoms: Atoms


class CalculatorOutputs(BaseModel):
    is_complete: bool = False
    is_converged: bool = False


GenericCalculatorInputs = TypeVar('GenericCalculatorInput', bound=CalculatorInputs)
GenericCalculatorOutputs = TypeVar('GenericCalculatorOutput', bound=CalculatorOutputs)


class CalculatorProcess(Process[GenericCalculatorInputs, GenericCalculatorOutputs], Generic[GenericCalculatorInputs, GenericCalculatorOutputs]):
    ext_in: ClassVar[str]
    ext_out: ClassVar[str]
    id: ClassVar[str]
    code: ClassVar[str]
    parent_process: Process | None

    def _run():
        raise NotImplementedError()

    def link_pseudopotential(self, src: str | Path, dest: str | Path) -> None:
        """Link a pseudopotential file to the calculator.

        This is engine-dependent because for some engines we don't want to simply link pseudopotential
        files like we do for other files.
        """
        if self.engine is None:
            raise ValueError('Cannot link a pseudopotential to a calculator without first setting its engine')
        if not isinstance(src, Path):
            src = Path(src)
        if not isinstance(dest, Path):
            dest = Path(dest)
        self.engine.link_pseudopotential(self, src, dest)


CalcType = TypeVar('CalcType', bound=CalculatorProcess)
