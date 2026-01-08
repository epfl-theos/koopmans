from abc import ABC, abstractmethod
from typing import Annotated, Literal

import numpy as np
from ase_koopmans.cell import Cell
from ase_koopmans.units import Bohr
from pydantic import BeforeValidator, field_validator, model_validator

from koopmans.cell import parameters_to_cell
from koopmans.pydantic import BaseModel


def tidy_units(value: str) -> str:
    value = value.lower()
    value = value.replace("angstrom", "ang")
    return value


class CellParametersABC(BaseModel, ABC):
    periodic: bool | tuple[bool, bool, bool] = True

    @abstractmethod
    def to_cell(self) -> Cell:
        ...


class CellParametersViaIbrav(CellParametersABC):

    ibrav: int
    celldms: dict[int, float]

    def to_cell(self) -> Cell:
        return parameters_to_cell(celldms=self.celldms, ibrav=self.ibrav)


class CellParametersViaAlat(CellParametersABC):
    celldms: dict[int, float]
    vectors: list[tuple[float, float, float]]
    units: Annotated[Literal['alat'], BeforeValidator(tidy_units)] = 'alat'

    def to_cell(self) -> Cell:
        alat = self.celldms.get(1, None)
        if alat is None:
            raise ValueError("When using `alat` as the units for the cell parameters, `celldms` must be provided.")
        cell = np.array(self.vectors) * alat * Bohr
        return Cell(cell)


class CellParametersViaVectors(CellParametersABC):

    vectors: list[tuple[float, float, float]]
    units: Annotated[Literal['bohr', 'ang'], BeforeValidator(tidy_units)] = 'ang'

    def to_cell(self) -> Cell:
        conversion_factor = 1.0 if self.units == 'ang' else Bohr
        cell = np.array(self.vectors) * conversion_factor
        return Cell(cell)


class CellParametersInput(BaseModel):
    periodic: bool = True
    ibrav: int | None = None
    celldms: dict[int, float] | None = None
    vectors: list[tuple[float, float, float]] | None = None
    units: Annotated[Literal['bohr', 'ang', 'alat'], BeforeValidator(tidy_units)] = 'ang'

    @model_validator(mode="after")
    def check_ibrav_and_celldms_xor_vectors(self):
        if (self.ibrav is None) != (self.celldms is None):
            raise ValueError("'ibrav' and 'celldms' must be provided together.")
        if (self.ibrav is None) == (self.vectors is None):
            raise ValueError("Specify either both 'ibrav' and 'celldms', or 'vectors'.")
        return self

    def to_cell(self) -> Cell:
        if self.vectors is None:
            return parameters_to_cell(celldms=self.celldms, ibrav=self.ibrav)
