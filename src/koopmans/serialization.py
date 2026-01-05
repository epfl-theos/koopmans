from typing import Any

import numpy as np
from ase import Atoms
from ase.cell import Cell
from ase.dft.kpoints import BandPath
from ase.spectrum.band_structure import BandStructure

from koopmans.utils import remove_numpy_from_dict


def serialize_ase_atoms(self: Atoms) -> dict[str, Any]:
    return remove_numpy_from_dict(self.todict())


def serialize_ase_cell(cell: Cell) -> dict[str, Any]:
    return remove_numpy_from_dict(cell.todict())


def serialize_ase_bandpath(bandpath: BandPath) -> dict[str, Any]:
    dct = bandpath.todict()
    dct.pop('labelseq')
    dct['cell'] = serialize_ase_cell(dct['cell'])
    return remove_numpy_from_dict(dct)


def serialize_ase_bandstructure(bandstructure: BandStructure) -> dict[str, Any]:
    dct = bandstructure.todict()
    dct['path'] = serialize_ase_bandpath(dct['path'])
    return remove_numpy_from_dict(dct)


def serialize_numpy(array: np.ndarray) -> list:
    return array.tolist()


json_encoders = {
    Atoms: serialize_ase_atoms,
    BandPath: serialize_ase_bandpath,
    BandStructure: serialize_ase_bandstructure,
    Cell: serialize_ase_cell,
    np.ndarray: serialize_numpy,
}
