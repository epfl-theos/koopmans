"""

Generic I/O functions that koopmans.calculators and koopmans.workflows can import non-cyclically

Written by Edward Linscott Jan 2020
Moved into utils Sep 2021

"""

from datetime import datetime
import json
import numpy as np
from typing import List, Union, Tuple
from pathlib import Path
from ase.atoms import Atoms
from ase.dft.kpoints import bandpath, BandPath
from ase.io.espresso.utils import label_to_symbol, label_to_tag
from ase.calculators.calculator import Calculator


def parse_dict(dct: dict) -> dict:
    '''

    Reads in a dict, formatting the values appropriately if they are not already

    '''
    settings = {}
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


def construct_cell_parameters_block(atoms: Atoms) -> dict:
    return {'vectors': [list(row) for row in atoms.cell[:]], 'units': 'angstrom'}


def construct_atomic_positions_block(atoms: Atoms) -> dict:
    labels = atoms.get_chemical_symbols()
    positions = atoms.get_scaled_positions()
    return {'positions': list([label] + [x for x in pos] for label, pos in zip(labels, positions)), 'units': 'crystal'}


def construct_atomic_species_block(atoms: Atoms) -> dict:
    labels = atoms.get_chemical_symbols()
    masses = atoms.get_masses()
    pseudopotentials = ['Si_ONCV_PBE-1.2.upf', 'Si_ONCV_PBE-1.2.upf']
    return {'species': list([label] + [m] + [pp] for label, m, pp in zip(labels, masses, pseudopotentials))}


def write_alpha_file(directory: Path, alphas: List[float], filling: List[bool]):
    a_filled = [a for a, f in zip(alphas, filling) if f]
    a_empty = [a for a, f in zip(alphas, filling) if not f]
    for alphas, suffix in zip([a_filled, a_empty], ['', '_empty']):
        with open(directory / f'file_alpharef{suffix}.txt', 'w') as fd:
            fd.write('{}\n'.format(len(alphas)))
            fd.writelines(['{} {} 1.0\n'.format(i + 1, a)
                           for i, a in enumerate(alphas)])


def read_alpha_file(directory: Path) -> List[float]:
    alphas = []
    for suffix in ['', '_empty']:
        fname = directory / f'file_alpharef{suffix}.txt'
        if not fname.is_file():
            break
        with open(fname, 'r') as fd:
            flines = fd.readlines()
            n_orbs = int(flines[0])
            alphas += [float(line.split()[1]) for line in flines[1:n_orbs + 1]]
    return alphas


def read_kpoints_block(calc: Calculator, dct: dict):
    for k, v in dct.items():
        if k in ['gamma_only', 'kgrid', 'koffset']:
            calc.parameters[k] = v
        elif k == 'kpath':
            read_kpath(calc, v)
        else:
            raise KeyError(f'Unrecognised option "{k}" provided in the k_points block')

    return


def read_kpath(calc: Calculator, kpath: Union[str, List[Tuple[float, float, float, int]]]):
    calc.atoms.cell.pbc = True
    if isinstance(kpath, str):
        # Interpret kpath as a string of points in the BZ
        npoints = 10 * len(kpath) - 9 - 29 * kpath.count(',')
        calc.parameters['kpath'] = bandpath(kpath, calc.atoms.cell, npoints=npoints)
    else:
        # Interpret bandpath as using PW syntax (https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm1290)
        kpts = []
        for k1, k2 in zip(kpath[:-1], kpath[1:]):
            # Remove the weights, storing the weight of k1
            npoints = k1[-1]
            # Interpolate the bandpath
            kpts += bandpath([k1[:-1], k2[:-1]], calc.atoms.cell, npoints + 1).kpts[:-1].tolist()
        # Don't forget about the final kpoint
        kpts.append(k2[:-1])
        if len(kpts) != sum([k[-1] for k in kpath[:-1]]) + 1:
            raise AssertionError(
                'Did not get the expected number of kpoints; this suggests there is a bug in the code')
        calc.parameters['kpath'] = BandPath(calc.atoms.cell, kpts)


def read_atomic_species(calc: Calculator, dct: dict):
    calc.parameters['pseudopotentials'] = {l[0]: l[2] for l in dct['species']}


def read_atomic_positions(calc: Calculator, dct: dict):

    pos_array = np.array(dct['positions'])
    symbols = [label_to_symbol(p) for p in pos_array[:, 0]]
    tags = [label_to_tag(p) for p in pos_array[:, 0]]
    positions = np.array(pos_array[:, 1:], dtype=float)

    scale_positions = False
    units = dct.get('units', 'angstrom').lower()
    if units == 'angstrom':
        pass
    elif units == 'crystal':
        scale_positions = True
    else:
        raise NotImplementedError(
            f'atomic_positions units = {units} is not yet implemented')

    if not calc.atoms.cell:
        raise ValueError('io.read_atomic_positions() must be called after io.read_cell_parameters()')

    if scale_positions:
        calc.atoms = Atoms(symbols, scaled_positions=positions, calculator=calc, cell=calc.atoms.cell)
    else:
        calc.atoms = Atoms(symbols, positions=positions, calculator=calc, cell=calc.atoms.cell)
    calc.atoms.set_tags(tags)


def read_cell_parameters(calc: Calculator, dct: dict):
    cell = dct.get('vectors', None)
    units = dct.get('units', None)
    if cell is None and units in [None, 'alat']:
        if 'ibrav' not in calc.parameters:
            raise KeyError('Cell has not been defined. Please specify either "ibrav" and related "celldm"s) '
                           ' or a "cell_parameters" block in "setup"')
    elif cell is not None and units == 'angstrom':
        pass

    else:
        raise NotImplementedError('the combination of vectors, ibrav, & units '
                                  'in the cell_parameter block cannot be read (may not yet be '
                                  'implemented)')
    return cell


print_call_end = '\n'


def indented_print(text: str = '', indent: int = 0, **kwargs):
    global print_call_end
    for substring in text.split('\n'):
        if print_call_end == '\n':
            print(' ' * indent + substring, **kwargs)
        else:
            print(substring, **kwargs)
    print_call_end = kwargs.get('end', '\n')


def write_hr_file(fname: Path, ham: np.ndarray, rvect: List[List[int]], weights: List[int]) -> None:

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


def read_hr_file(fname: Path) -> Tuple[np.ndarray, np.ndarray, List[int], int]:
    """
    Reads in a hr file, but does not reshape the hamiltonian (because we want to reshape different Hamiltonians
    differently)
    """

    with open(fname, 'r') as fd:
        lines = fd.readlines()

    if 'written on' in lines[0].lower():
        pass
    elif 'xml version' in lines[0] or fname == 'hamiltonian_emp.dat':
        raise ValueError(f'The format of {fname} is no longer supported')
    else:
        raise ValueError(f'The format of {fname} is not recognised')

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
