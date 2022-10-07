"""

Generic I/O functions that koopmans.calculators and koopmans.workflows can import non-cyclically

Written by Edward Linscott Jan 2020
Moved into utils Sep 2021

"""

import json
from datetime import datetime
from pathlib import Path
from typing import IO, Any, Dict, List, Tuple, Union

import numpy as np
import numpy.typing as npt
from ase.atoms import Atoms
from ase.io.espresso import label_to_symbol, label_to_tag
from ase.units import Bohr

from koopmans.cell import (cell_follows_qe_conventions, cell_to_parameters,
                           parameters_to_cell)


def parse_dict(dct: Dict[str, Any]) -> Dict[str, Any]:
    '''

    Reads in a dict, formatting the values appropriately if they are not already

    '''
    settings: Dict[str, Any] = {}
    for k, v in dct.items():
        # Deal with bools separately since JSON strictly only will interpret
        # 'false' as False, while 'False' will be left as a string and
        # any statement to the effect of 'if param' will evaluate to True if
        # param = 'False'
        if isinstance(v, str) and v.lower() in ['f', 'false']:
            settings[k] = False
        elif isinstance(v, str) and v.lower() in ['t', 'true']:
            settings[k] = True
        else:
            try:
                settings[k] = json.loads(v)
            except (TypeError, json.decoder.JSONDecodeError) as e:
                settings[k] = v
    return settings


def construct_cell_parameters_block(atoms: Atoms) -> Dict[str, Any]:
    params: Dict[str, Any]
    if cell_follows_qe_conventions(atoms.cell):
        params = cell_to_parameters(atoms.cell)
    else:
        params = {'vectors': [list(row) for row in atoms.cell[:]], 'units': 'angstrom'}
    params['periodic'] = all(atoms.pbc)
    return params


def construct_atomic_positions_block(atoms: Atoms, crystal: bool = True) -> dict:
    if len(set(atoms.get_tags())) > 1:
        labels = [s + str(t) if t > 0 else s for s, t in zip(atoms.symbols, atoms.get_tags())]
    else:
        labels = atoms.symbols
    if not crystal:
        dct = {'positions': [
            [label] + [str(x) for x in pos] for label, pos in zip(labels, atoms.get_positions())],
            'units': 'angstrom'}
    else:
        dct = {'positions': [
            [label] + [str(x) for x in pos] for label, pos in zip(labels, atoms.get_scaled_positions())],
            'units': 'crystal'}
    return dct


def construct_atomic_species_block(atoms: Atoms) -> Dict[str, Any]:
    labels: List[str] = atoms.get_chemical_symbols()
    masses: List[float] = atoms.get_masses()
    pseudopotentials = ['Si_ONCV_PBE-1.2.upf', 'Si_ONCV_PBE-1.2.upf']
    species: List[List[Union[str, float]]] = []
    for label, m, pp in zip(labels, masses, pseudopotentials):
        line: List[Union[str, float]] = [label]
        line.append(m)
        line.append(pp)
        species.append(line)

    return {'species': species}


def write_alpha_file(directory: Path, alphas: List[float], filling: List[bool]):
    a_filled = [a for a, f in zip(alphas, filling) if f]
    a_empty = [a for a, f in zip(alphas, filling) if not f]
    for alphas, suffix in zip([a_filled, a_empty], ['', '_empty']):
        with open(directory / f'file_alpharef{suffix}.txt', 'w') as fd:
            fd.write('{}\n'.format(len(alphas)))
            fd.writelines(['{} {} 1.0\n'.format(i + 1, a)
                           for i, a in enumerate(alphas)])


def read_alpha_file(directory: Path) -> List[float]:
    alphas: List[float] = []
    for suffix in ['', '_empty']:
        fname = directory / f'file_alpharef{suffix}.txt'
        if not fname.is_file():
            break
        with open(fname, 'r') as fd:
            flines = fd.readlines()
            n_orbs = int(flines[0])
            alphas += [float(line.split()[1]) for line in flines[1:n_orbs + 1]]
    return alphas


def read_atomic_positions(atoms: Atoms, dct: Dict[str, Any]):

    pos_array = np.array(dct['positions'])
    symbols: List[str] = [label_to_symbol(p) for p in pos_array[:, 0]]
    tags = [label_to_tag(p) for p in pos_array[:, 0]]
    positions = np.array(pos_array[:, 1:], dtype=float)

    scale_positions = False
    units = dct.get('units', 'angstrom').lower()
    if units == 'angstrom':
        pass
    elif units == 'bohr':
        positions /= Bohr
    elif units == 'alat':
        celldms = cell_to_parameters(atoms.cell).get('celldms')
        assert isinstance(celldms, dict)
        positions *= celldms[1] / Bohr
    elif units == 'crystal':
        scale_positions = True
    else:
        raise NotImplementedError(
            f'atomic_positions units = {units} is not yet implemented')

    if not atoms.cell:
        raise ValueError('io.read_atomic_positions() must be called after io.read_cell_parameters()')

    assert len(atoms) == 0, 'Atoms should be length zero at this stage'
    if scale_positions:
        atoms += Atoms(symbols, scaled_positions=positions, cell=atoms.cell)
    else:
        atoms += Atoms(symbols, positions=positions, cell=atoms.cell)
    atoms.set_tags(tags)


def read_cell_parameters(atoms: Atoms, dct: Dict[str, Any]):
    cell = dct.pop('vectors', None)
    units = dct.pop('units', '')
    atoms.pbc = dct.pop('periodic', True)
    if cell is None:
        if 'ibrav' not in dct:
            raise KeyError('Cell has not been defined. Please specify either "ibrav" and related "celldm"s) '
                           ' or a "cell_parameters" block in "setup"')
        celldms = {int(k): v for k, v in dct.pop('celldms', {}).items()}
        cell = parameters_to_cell(celldms=celldms, **dct)
    elif units.lower() == 'angstrom':
        pass
    elif units.lower() == 'bohr':
        cell = np.array(cell) / Bohr
    elif units.lower() == 'alat':
        alat = dct.get('celldms', {}).get(1, None)
        if alat is None:
            raise ValueError('Please provide celldm(1) for a cell specified in units of alat')
        cell = np.array(cell) * alat / Bohr
    else:
        raise ValueError('The combination of vectors, ibrav, & units in the cell_parameter block is not valid')
    atoms.cell = cell
    return


print_call_end = '\n'


def indented_print(text: str = '', indent: int = 0, sep: str = ' ', end: str = '\n',
                   flush: bool = False):
    global print_call_end
    for substring in text.split('\n'):
        if print_call_end == '\n':
            print(' ' * indent + substring, sep=sep, end=end, flush=flush)
        else:
            print(substring, sep=sep, end=end, flush=flush)
    print_call_end = end


def write_wannier_hr_file(fname: Path, ham: np.ndarray, rvect: List[List[int]], weights: List[int]) -> None:

    nrpts = len(rvect)
    num_wann = np.size(ham, -1)
    expected_shape = (nrpts, num_wann, num_wann)
    if ham.shape != expected_shape:
        raise ValueError(f'ham has shape {ham.shape} which does not match the expected shape {expected_shape}')

    flines = [f' Written on {datetime.now().isoformat(timespec="seconds")}']
    flines.append(f'{num_wann:12d}')
    flines.append(f'{nrpts:12d}')

    ints_per_line = 15
    for pos in range(0, len(weights), ints_per_line):
        flines.append(''.join([f'{x:5d}' for x in weights[pos:pos + ints_per_line]]))

    for r, ham_block in zip(rvect, ham):
        flines += [f'{r[0]:5d}{r[1]:5d}{r[2]:5d}{j+1:5d}{i+1:5d}{val.real:12.6f}{val.imag:12.6f}' for i,
                   row in enumerate(ham_block) for j, val in enumerate(row)]

    # Make sure the parent directory exists
    fname.parent.mkdir(exist_ok=True, parents=True)

    # Write the Hamiltonian to file
    with open(fname, 'w') as fd:
        fd.write('\n'.join(flines))


def read_wannier_hr_file(fname: Path) -> Tuple[np.ndarray, np.ndarray, List[int], int]:
    """
    Reads in a hr file, but does not reshape the hamiltonian (because we want to reshape different Hamiltonians
    differently)

    Returns a tuple containing...
        - the hamiltonian
        - the r-vectors
        - the list of weights
        - the number of wannier functions

    """

    with open(fname, 'r') as fd:
        lines = fd.readlines()

    if 'written on' in lines[0].lower():
        pass
    elif 'xml version' in lines[0] or fname == 'hamiltonian_emp.dat':
        raise ValueError(f'The format of {fname} is no longer supported')
    else:
        raise ValueError(f'The format of {fname} is not recognized')

    # Read in the number of r-points and the number of Wannier functions
    nrpts = int(lines[2].split()[0])
    single_R = (nrpts == 1)

    if not single_R:
        num_wann = int(lines[1].split()[0])

    lines_to_skip = 3 + nrpts // 15
    if nrpts % 15 > 0:
        lines_to_skip += 1

    # Read in the weights
    weights = [int(x) for line in lines[3:lines_to_skip] for x in line.split()]

    # Read in the hamiltonian and the unique r-vectors
    hr: List[complex] = []
    rvect: List[List[int]] = []
    for i, line in enumerate(lines[lines_to_skip:]):
        hr.append(float(line.split()[5]) + 1j * float(line.split()[6]))
        if not single_R and i % num_wann**2 == 0:
            rvect.append([int(x) for x in line.split()[0:3]])

    # Convert hr and rvect to numpy arrays
    hr_np = np.array(hr, dtype=complex)
    if single_R:
        rvect_np = np.array([[0, 0, 0]])
    else:
        rvect_np = np.array(rvect, dtype=int)

    return hr_np, rvect_np, weights, nrpts


def read_wannier_u_file(fname: Path) -> Tuple[npt.NDArray[np.complex_], npt.NDArray[np.float_], int]:

    with open(fname, 'r') as fd:
        lines = fd.readlines()

    nk, m, n = [int(x) for x in lines[1].split()]

    kpts = np.empty((nk, 3), dtype=float)
    umat = np.empty((nk, m, n), dtype=complex)

    for i, line in enumerate(lines[3:]):
        ik = i // (2 + m * n)
        if i % (2 + m * n) != 0:
            continue
        kpts[ik, :] = line.split()
        umat[ik, :, :] = np.reshape([complex(*[float(x) for x in line.split()])
                                     for line in lines[i + 4:i + 4 + m * n]], (m, n))

    return umat, kpts, nk


def write_wannier_u_file(fname: Path, umat: npt.NDArray[np.complex_], kpts: npt.NDArray[np.float_]):

    flines = [f' Written on {datetime.now().isoformat(timespec="seconds")}']
    flines.append(''.join([f'{x:12d}' for x in umat.shape]))

    for kpt, umatk in zip(kpts, umat):
        flines.append('')
        flines.append(''.join([f'{k:15.10f}' for k in kpt]))
        flines += [f'{c.real:15.10f}{c.imag:15.10f}' for c in umatk.flatten()]

    with open(fname, 'w') as fd:
        fd.write('\n'.join(flines))


def read_wannier_centers_file(fname: Path):

    with open(fname, 'r') as fd:
        lines = fd.readlines()

    centers = []
    symbols = []
    positions = []
    for line in lines[2:]:
        if line.startswith('X    '):
            centers.append([float(x) for x in line.split()[1:]])
        else:
            symbols.append(line.split()[0])
            positions.append([float(x) for x in line.split()[1:]])
    return centers, Atoms(symbols=symbols, positions=positions, pbc=True)


def write_wannier_centers_file(fname: Path, centers: List[List[float]], atoms: Atoms):
    length = len(centers) + len(atoms)

    # Add the header
    flines = [f'{length:6d}',
              f' Wannier centers, written by koopmans on {datetime.now().isoformat(timespec="seconds")}']

    # Add the centers
    for center in centers:
        flines.append('X    ' + ''.join([f'{x:16.8f}' for x in center]))

    # Add the atoms
    for atom in atoms:
        flines.append(f'{atom.symbol: <5}' + ''.join([f'{x:16.8f}' for x in atom.position]))

    # Write to file
    with open(fname, 'w') as fd:
        fd.write('\n'.join(flines))
