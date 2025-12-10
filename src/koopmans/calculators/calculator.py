from pathlib import Path
from typing import Any, ClassVar, Generic, TypeVar

import numpy as np
from ase_koopmans import Atoms
from pydantic import field_serializer

from koopmans.processes.process import Process
from koopmans.pydantic import BaseModel


def convert_numpy(obj: Any) -> Any:
    """
    Recursively convert NumPy arrays to Python lists inside nested structures.
    Supports dict, list, tuple, set, and nested combinations.
    """
    if isinstance(obj, dict):
        return {k: convert_numpy(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy(v) for v in obj]
    elif isinstance(obj, tuple):
        return tuple(convert_numpy(v) for v in obj)
    elif isinstance(obj, set):
        return {convert_numpy(v) for v in obj}
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj


class CalculatorInputs(BaseModel):
    atoms: Atoms

    @field_serializer("atoms")
    def serialize_atoms(self, v: Atoms) -> dict[str, Any]:
        return convert_numpy(v.todict())


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
