"""

An extended atoms class that has the additional properties
  acell: the cell, normalized by alat
  asupercell: a supercell, normalized by alat
for use as the Atoms object associated with a UI calculator

"""

from typing import Optional

import numpy as np
from ase import Atoms
from ase.cell import Cell


class UIAtoms(Atoms):

    __supercell_matrix: Optional[np.ndarray]

    def __init__(self, supercell_matrix: Optional[np.ndarray] = None, *args, **kwargs):
        if supercell_matrix is not None:
            self._supercell_matrix = supercell_matrix
        super().__init__(*args, **kwargs)

    @property
    def supercell(self) -> Cell:
        # The supercell
        return Cell(self._supercell_matrix @ self.cell)

    @property
    def acell(self) -> Cell:
        # The cell, normalized by alat
        return Cell(self.cell / np.linalg.norm(self.cell[0]))

    @property
    def asupercell(self) -> Cell:
        # The supercell, normalized by alat
        return Cell(self.supercell / np.linalg.norm(self.cell[0]))

    @property
    def _supercell_matrix(self):
        if self.__supercell_matrix is None:
            raise AttributeError(f'{self}._supercell_matrix has not yet been set')
        else:
            return self.__supercell_matrix

    @_supercell_matrix.setter
    def _supercell_matrix(self, value: np.ndarray):
        self.__supercell_matrix = value

    # @classmethod
    # def fromdict(cls, dct):
    #     sm = dct.pop('supercell_matrix')
    #     atoms = super(UIAtoms, cls).fromdict(dct)
    #     atoms._supercell_matrix = sm
    #     return atoms

    @classmethod
    def fromatoms(cls, atoms: Atoms, supercell_matrix: Optional[np.ndarray] = None):
        # Crude method for converting an Atoms object to a UIAtoms object by taking advantage
        # of Atoms.todict() and Atoms.fromdict
        # > ui_atoms = UIAtoms.fromatoms(...)

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
