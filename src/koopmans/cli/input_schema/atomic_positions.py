from abc import ABC, abstractmethod
from pathlib import Path
from typing import Literal

import numpy as np
from ase_koopmans.io import read as ase_read
from ase_koopmans.io.espresso import label_to_symbol, label_to_tag
from ase_koopmans.units import Bohr
from numpy.typing import NDArray
from pydantic import PrivateAttr

from koopmans.pydantic import BaseModel


class AtomicPositionsABC(BaseModel, ABC):
    positions: list[tuple[str, float, float, float]]
    snapshots: Path | None = None

    @property
    def symbols(self) -> list[str]:
        return [label_to_symbol(p[0]) for p in self.positions]

    @property
    def tags(self) -> list[int]:
        return [label_to_tag(p[0]) for p in self.positions]

    def positions_array(self) -> NDArray[np.float64]:
        return np.array([p[1:] for p in self.positions])

    def _convert_positions(self, conversion_factor: float) -> NDArray[np.float64]:
        return self.positions_array() * conversion_factor

    @abstractmethod
    def in_angstrom(self) -> NDArray[np.float64]:
        ...


class AtomicAlatPositionsInput(AtomicPositionsABC):
    units: Literal['alat'] = 'alat'
    _alat: float | None = PrivateAttr(default=None)

    def in_angstrom(self) -> NDArray[np.float64]:
        if self._alat is None:
            raise ValueError("`alat` must be set before converting to angstroms.")
        return self._convert_positions(self._alat * Bohr)


class AtomicPositionsInput(AtomicPositionsABC):
    units: Literal['crystal', 'ang', 'bohr']
    snapshots: Path | None = None

    def in_angstrom(self) -> NDArray[np.float64]:
        if self.units == 'crystal':
            raise ValueError('Cannot convert fractional coordinates to angstroms without cell information.')
        conversion_factor = 1.0 if self.units == 'ang' else Bohr
        return self._convert_positions(conversion_factor)
