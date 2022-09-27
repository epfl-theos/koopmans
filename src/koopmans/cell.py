import copy
from typing import Dict, Union

import numpy as np
from ase.cell import Cell
from ase.lattice import (BCC, BCT, CUB, FCC, HEX, MCL, MCLC, ORC, ORCC, ORCF,
                         ORCI, RHL, TET, TRI, BravaisLattice)
from ase.units import Bohr


def parameters_to_cell(ibrav: int, celldms: Dict[int, float]) -> Cell:
    '''
    Converts QE cell parameters to the corresponding ASE `Cell` object.

    N.B. `parameters_to_cell` and `cell_to_parameters` are not idempotent for ibravs 12 to 14 because
    ASE's Bravais lattice detection is not reliable

    Parameters
    ----------
    ibrav : int
        the Bravais lattice index
    celldms: Dict[int, float]
        the crystallographic constrains

    Returns
    -------
    Cell
        An ASE `Cell` object corresponding to the provided ibrav and celldms

    '''

    if not isinstance(celldms, dict):
        raise ValueError('Please provide celldms as a dictionary e.g. "celldms": {"1": 10.0, "3": 0.75}')

    # Convert from Bohr to Angstrom
    celldms = copy.deepcopy(celldms)
    celldms[1] *= Bohr

    if ibrav == 1:
        new_array = celldms[1] * np.identity(3)
    elif ibrav == 2:
        new_array = celldms[1] / 2 * np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]])
    elif ibrav == 3:
        new_array = celldms[1] / 2 * np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1]])
    elif ibrav == 4:
        new_array = celldms[1] * np.array(
            [[1, 0, 0], [-0.5, np.sqrt(3) / 2, 0], [0, 0, celldms[3]]]
        )
    elif ibrav == 5:
        c = celldms[4]
        tx = np.sqrt((1 - c) / 2)
        ty = np.sqrt((1 - c) / 6)
        tz = np.sqrt((1 + 2 * c) / 3)
        new_array = celldms[1] * np.array(
            [[tx, -ty, tz], [0, 2 * ty, tz], [-tx, -ty, tz]]
        )
    elif ibrav == 6:
        new_array = celldms[1] * np.array([[1, 0, 0], [0, 1, 0], [0, 0, celldms[3]]])
    elif ibrav == 7:
        # ASE and QE use different conventions
        a = celldms[1]
        c = celldms[3] * celldms[1]
        new_array = 0.5 * np.array([[a, -a, c], [a, a, c], [-a, -a, c]])
    elif ibrav == 8:
        new_array = celldms[1] * np.array(
            [[1, 0, 0], [0, celldms[2], 0], [0, 0, celldms[3]]]
        )
    elif ibrav == 9:
        a = celldms[1]
        b = celldms[1] * celldms[2]
        c = celldms[1] * celldms[3]
        new_array = np.array([[a / 2, b / 2, 0], [-a / 2, b / 2, 0], [0, 0, c]])
    elif ibrav == 10:
        # ASE and QE use different conventions
        a = celldms[1]
        b = celldms[1] * celldms[2]
        c = celldms[1] * celldms[3]
        new_array = 0.5 * np.array([[a, 0, c], [a, b, 0], [0, b, c]])
    elif ibrav == 11:
        a = celldms[1]
        b = celldms[1] * celldms[2]
        c = celldms[1] * celldms[3]
        new_array = 0.5 * np.array([[a, b, c], [-a, b, c], [-a, -b, c]])
    elif ibrav == 12:
        a = celldms[1]
        b = celldms[1] * celldms[2]
        c = celldms[1] * celldms[3]
        new_array = np.array(
            [[a, 0, 0], [b * celldms[4], b * np.sqrt(1 - celldms[4] ** 2), 0], [0, 0, c]]
        )
    elif ibrav == 13:
        a = celldms[1]
        b = celldms[1] * celldms[2]
        c = celldms[1] * celldms[3]
        new_array = np.array(
            [
                [a / 2, 0, -c / 2],
                [b * celldms[4], b * np.sqrt(1 - celldms[4] ** 2), 0],
                [a / 2, 0, c / 2],
            ]
        )
    elif ibrav == 14:
        a = celldms[1]
        b = celldms[1] * celldms[2]
        c = celldms[1] * celldms[3]
        cosalpha = celldms[4]
        cosbeta = celldms[5]
        cosgamma = celldms[6]
        singamma = np.sqrt(1 - cosgamma ** 2)

        a3x = c * cosbeta
        a3y = c / singamma * (cosalpha - cosbeta * cosgamma)
        a3z = (
            c
            / singamma
            * np.sqrt(
                singamma**2
                - cosalpha**2
                - cosbeta**2
                + 2 * cosalpha * cosbeta * cosgamma
            )
        )
        new_array = np.array(
            [[a, 0, 0], [b * cosgamma, b * singamma, 0], [a3x, a3y, a3z]]
        )
    else:
        raise ValueError(f"Unrecognized ibrav {ibrav}")

    return Cell(new_array)


def cell_to_parameters(cell: Cell) -> Dict[str, Union[int, Dict[int, float]]]:
    '''
    Identifies a cell in the language of Quantum ESPRESSO

    Parameters
    ----------
    cell : Cell
        an ASE `Cell`

    Returns
    -------
    Dict
        a dictionary containing the ibrav and celldms

    '''

    # By regenerating the cell we ensure we are always working with the Niggli-reduced lattice
    # Beware: ASE's get_bravais_lattice().tocell() are not idempotent, because tocell() imposes AFlow conventions
    lat: BravaisLattice = cell.get_bravais_lattice()
    new_cell = lat.tocell()

    if abs(cell.volume - new_cell.volume) > 1e-8:
        raise ValueError('You have provided a cell that is not Niggli-reduced.\n'
                         'Try running\n'
                         ' >>> cell.get_bravais_lattice().tocell()\n'
                         'within python to obtain a reduced cell')
    celldms: Dict[int, float] = {}
    [A, B, C, _, _, gamma] = new_cell.cellpar(radians=True)
    a: float
    b: float
    c: float

    if isinstance(lat, CUB):
        ibrav = 1
        celldms[1] = A
    elif isinstance(lat, FCC):
        ibrav = 2
        celldms[1] = np.sqrt(2) * A
    elif isinstance(lat, BCC):
        ibrav = 3
        celldms[1] = 2 * A / np.sqrt(3)
    elif isinstance(lat, HEX):
        ibrav = 4
        celldms[1] = A
        celldms[3] = C / A
    elif isinstance(lat, RHL):
        ibrav = 5
        celldms[1] = A
        celldms[4] = np.cos(gamma)
    elif isinstance(lat, TET):
        ibrav = 6
        celldms[1] = A
        celldms[3] = C / A
    elif isinstance(lat, BCT):
        # ASE and QE use different conventions
        ibrav = 7
        a = lat.a
        c = lat.c
        celldms[1] = a
        celldms[3] = c / a
    elif isinstance(lat, ORC):
        ibrav = 8
        celldms[1] = A
        celldms[2] = B / A
        celldms[3] = C / A
    elif isinstance(lat, ORCC):
        ibrav = 9
        a = lat.a
        b = lat.b
        c = lat.c
        celldms[1] = a
        celldms[2] = b / a
        celldms[3] = c / a
    elif isinstance(lat, ORCF):
        # ASE and QE use different conventions
        ibrav = 10
        a = lat.a
        b = lat.b
        c = lat.c
        celldms[1] = a
        celldms[2] = b / a
        celldms[3] = c / a
    elif isinstance(lat, ORCI):
        ibrav = 11
        a = lat.a
        b = lat.b
        c = lat.c
        celldms[1] = a
        celldms[2] = b / a
        celldms[3] = c / a
    elif isinstance(lat, MCL):
        ibrav = 12
        a = lat.b
        b = lat.c
        c = lat.a
        celldms[1] = a
        celldms[2] = b / a
        celldms[3] = c / a
        celldms[4] = np.cos(lat.alpha * np.pi / 180)
    elif isinstance(lat, MCLC):
        ibrav = 13
        a = lat.b
        b = lat.c
        c = lat.a
        celldms[1] = a
        celldms[2] = b / a
        celldms[3] = c / a
        celldms[4] = np.cos(lat.alpha * np.pi / 180)
    elif isinstance(lat, TRI):
        ibrav = 14
        celldms[1] = lat.a
        celldms[2] = lat.b / lat.a
        celldms[3] = lat.c / lat.a
        cosgamma = np.cos(lat.gamma * np.pi / 180)
        cosbeta = np.cos(lat.beta * np.pi / 180)
        cosalpha = np.cos(lat.alpha * np.pi / 180)
        celldms[4] = cosalpha
        celldms[5] = cosbeta
        celldms[6] = cosgamma
    else:
        raise ValueError(f"Unrecognised Bravais lattice {lat.name}")

    # Convert to Bohr radii
    celldms[1] /= Bohr

    return {'ibrav': ibrav, 'celldms': celldms}


def cell_follows_qe_conventions(cell: Cell) -> bool:
    qe_cell = parameters_to_cell(**cell_to_parameters(cell))  # type: ignore
    return np.allclose(cell, qe_cell)
