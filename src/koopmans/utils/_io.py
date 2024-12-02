"""

Generic I/O functions that koopmans.calculators and koopmans.workflows can import non-cyclically

Written by Edward Linscott Jan 2020
Moved into utils Sep 2021

"""

import json
import sys
import textwrap
from datetime import datetime
from pathlib import Path
from typing import IO, Any, Dict, List, Optional, Tuple, Union

import numpy as np
import numpy.typing as npt
from ase_koopmans.atoms import Atoms
from ase_koopmans.io.espresso import label_to_symbol, label_to_tag
from ase_koopmans.units import Bohr

from koopmans.cell import (cell_follows_qe_conventions, cell_to_parameters,
                           parameters_to_cell)

from ._os import HasDirectory, get_content


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
        params = dict(**cell_to_parameters(atoms.cell))
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
        positions *= Bohr
    elif units == 'alat':
        celldms = cell_to_parameters(atoms.cell).get('celldms')
        assert isinstance(celldms, dict)
        positions *= celldms[1] * Bohr
    elif units == 'crystal':
        scale_positions = True
    else:
        raise NotImplementedError(
            f'atomic_positions `units = {units}` is not yet implemented')

    if not atoms.cell:
        raise ValueError('`io.read_atomic_positions()` must be called after `io.read_cell_parameters()`')

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
            raise KeyError('Cell has not been defined. Please specify either `ibrav` and related `celldms`) '
                           ' or a `cell_parameters` block')
        celldms = {int(k): v for k, v in dct.pop('celldms', {}).items()}
        cell = parameters_to_cell(celldms=celldms, **dct)
    elif units.lower() == 'angstrom':
        pass
    elif units.lower() == 'bohr':
        cell = np.array(cell) * Bohr
    elif units.lower() == 'alat':
        alat = dct.get('celldms', {}).get(1, None)
        if alat is None:
            raise ValueError('Please provide `celldm(1)` for a cell specified in units of `alat`')
        cell = np.array(cell) * alat * Bohr
    else:
        raise ValueError('The combination of vectors, ibrav, & units in the cell_parameter block is not valid')
    atoms.cell = cell
    return


print_indent = 0
print_call_end = '\n'
previous_indent = 0


def indented_print(text: str = '', indent: Optional[int] = None, style='body', parse_asterisks=True, flush=True, wrap=True, **kwargs: Any):
    if sys.stdout.isatty():
        if style == 'heading' and '**' not in text:
            text = f'**{text}**'
        if parse_asterisks:
            while '**' in text:
                text = text.replace('**', '\033[1m', 1).replace('**', '\033[0m', 1)
            while '*' in text:
                text = text.replace('*', '\033[3m', 1).replace('*', '\033[0m', 1)
    global print_indent
    indent = print_indent if indent is None else indent
    if style == 'body':
        _indented_print(text, indent=indent, flush=flush, wrap=wrap, **kwargs)
    elif style == 'heading':
        assert kwargs.get('end', '\n') == '\n'
        _indented_print(text, indent=indent, flush=flush, wrap=wrap, **kwargs)
    else:
        raise ValueError(f'Invalid choice `{style}` for style; must be `heading`/`body`')


def _indented_print(text: str = '', indent: int = 0, sep: str = ' ', end: str = '\n',
                    flush: bool = False, initial_indent: str | None = None, subsequent_indent: str | None = None, wrap=True):
    global print_call_end
    global previous_indent
    if indent < 0:
        indent = previous_indent
    for substring in text.split('\n'):
        if print_call_end in ['\n', '\r']:
            initial_indent = ' ' * indent if initial_indent is None else initial_indent
            subsequent_indent = ' ' * indent if subsequent_indent is None else subsequent_indent
            width = 120 if wrap else len(text) + indent
            message = textwrap.fill(substring, width=width, initial_indent=initial_indent,
                                    subsequent_indent=subsequent_indent, break_long_words=False, break_on_hyphens=False,
                                    drop_whitespace=False)
            print(message, sep=sep, end=end, flush=flush)
        else:
            print(substring, sep=sep, end=end, flush=flush)
    print_call_end = end
    previous_indent = indent


def print_alert(kind, message, header=None, indent=-1, **kwargs):
    global previous_indent
    allowed_kinds = {'note': "ℹ️ ", 'tip': "💡", 'important': "❕", 'warning': "🚨", 'caution': "❗"}
    if kind not in allowed_kinds:
        raise ValueError('`kind` must be one of ' + '/'.join(allowed_kinds.keys()))
    if sys.stdout.isatty():
        if indent < 0:
            indent = previous_indent
        header = "" if header is None else header + ': ' if message else header
        indented_print('\n' + allowed_kinds[kind] + ' ' + header + message + '\n', indent, **kwargs)
    else:
        if indent >= 0:
            width = 120 - indent
        else:
            width = 120 - previous_indent
        header = "" if header is None else header
        message = "\n".join(["",
                             f'> [!{kind.upper()}] {header} ',
                             textwrap.fill(str(message), width=width, initial_indent='> ', subsequent_indent='> '),
                             ""
                             ])
        indented_print(message, indent=indent)


def generate_wannier_hr_file_contents(ham: np.ndarray, rvect: List[List[int]], weights: List[int]) -> List[str]:

    nrpts = len(rvect)
    num_wann = np.size(ham, -1)
    expected_shape = (nrpts, num_wann, num_wann)
    if ham.shape != expected_shape:
        raise ValueError(f'`ham` has shape {ham.shape} which does not match the expected shape {expected_shape}')

    flines = [f' Written on {datetime.now().isoformat(timespec="seconds")}']
    flines.append(f'{num_wann:12d}')
    flines.append(f'{nrpts:12d}')

    ints_per_line = 15
    for pos in range(0, len(weights), ints_per_line):
        flines.append(''.join([f'{x:5d}' for x in weights[pos:pos + ints_per_line]]))

    for r, ham_block in zip(rvect, ham):
        flines += [f'{r[0]:5d}{r[1]:5d}{r[2]:5d}{j+1:5d}{i+1:5d}{val.real:12.6f}{val.imag:12.6f}' for i,
                   row in enumerate(ham_block) for j, val in enumerate(row)]

    return flines


def parse_wannier_hr_file_contents(lines: List[str]) -> Tuple[np.ndarray, np.ndarray, List[int], int]:
    """ Parse the contents of a Hamiltonian file

    Returns a tuple containing...
        - the hamiltonian (not reshaped, because we want to reshape different Hamiltonians differently)
        - the r-vectors
        - the list of weights
        - the number of wannier functions
    """

    if 'written on' in lines[0].lower():
        pass
    elif 'xml version' in lines[0]:
        raise ValueError(f'The format of Hamiltonian file contents no longer supported')
    else:
        raise ValueError(f'The format of the Hamiltonian file contents are not recognized')

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


def read_wannier_hr_file(src_calc_path: Tuple[HasDirectory, Path]) -> Tuple[np.ndarray, np.ndarray, List[int], int]:
    lines = get_content(*src_calc_path)
    return parse_wannier_hr_file_contents(lines)


def parse_wannier_u_file_contents(lines: List[str]) -> Tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64], int]:

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


def generate_wannier_u_file_contents(umat: npt.NDArray[np.complex128], kpts: npt.NDArray[np.float64]) -> List[str]:

    flines = [f' Written on {datetime.now().isoformat(timespec="seconds")}']
    flines.append(''.join([f'{x:12d}' for x in umat.shape]))

    for kpt, umatk in zip(kpts, umat):
        flines.append('')
        flines.append(''.join([f'{k:15.10f}' for k in kpt]))
        flines += [f'{c.real:15.10f}{c.imag:15.10f}' for c in umatk.flatten()]

    return flines


def parse_wannier_centers_file_contents(lines: List[str]) -> Tuple[List[List[float]], Atoms]:

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


def generate_wannier_centers_file_contents(centers: List[List[float]], atoms: Atoms) -> List[str]:
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

    return flines
