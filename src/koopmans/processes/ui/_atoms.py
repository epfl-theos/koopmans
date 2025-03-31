"""An extended atoms class for use as the Atoms object associated with a UI calculator."""

from typing import Optional

import numpy as np
from ase_koopmans import Atoms
from ase_koopmans.cell import Cell


class UIAtoms(Atoms):
    """An extended atoms class that has information about a supercell."""

    __supercell_matrix: Optional[np.ndarray]

    def __init__(self, supercell_matrix: Optional[np.ndarray] = None, *args, **kwargs):
        if supercell_matrix is not None:
            self._supercell_matrix = supercell_matrix
        super().__init__(*args, **kwargs)

    @property
    def supercell(self) -> Cell:
        """Return the supercell."""
        return Cell(self._supercell_matrix @ self.cell)

    @property
    def acell(self) -> Cell:
        """Return the cell, normalized by alat."""
        return Cell(self.cell / np.linalg.norm(self.cell[0]))

    @property
    def asupercell(self) -> Cell:
        """Return the supercell, normalized by alat."""
        return Cell(self.supercell / np.linalg.norm(self.cell[0]))

    @property
    def _supercell_matrix(self):
        if self.__supercell_matrix is None:
            raise AttributeError(f'`{self}._supercell_matrix` has not yet been set')
        else:
            return self.__supercell_matrix

    @_supercell_matrix.setter
    def _supercell_matrix(self, value: np.ndarray):
        self.__supercell_matrix = value

    @classmethod
    def fromatoms(cls, atoms: Atoms, supercell_matrix: Optional[np.ndarray] = None):
        """Convert an atoms object to a UIAtoms object by taking advantage of Atoms.todict() and Atoms.fromdict."""
        # Convert to a dict
        dct = atoms.todict()

        # Generate a new UIAtoms object from the dict
        ui_atoms = cls.fromdict(dct)

        # Add the supercell matrix information
        if supercell_matrix is not None:
            ui_atoms._supercell_matrix = supercell_matrix

        # Reconnect the calculator
        if atoms.calc is not None:
            ui_atoms.calc = atoms.calc
            ui_atoms.calc.atoms = ui_atoms

        # Return the new UIAtoms object
        return ui_atoms
