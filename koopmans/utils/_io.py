"""

Generic I/O functions that koopmans.calculators and koopmans.workflows can import non-cyclically

Written by Edward Linscott Jan 2020
Moved into utils Sep 2021

"""

import os
import json
import numpy as np
from typing import List, Union, Tuple
from ase.atoms import Atoms
from ase.dft.kpoints import bandpath, BandPath
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


def write_alpha_file(directory: str, alphas: List[float], filling: List[bool]):
    a_filled = [a for a, f in zip(alphas, filling) if f]
    a_empty = [a for a, f in zip(alphas, filling) if not f]
    for alphas, suffix in zip([a_filled, a_empty], ['', '_empty']):
        with open(f'{directory}/file_alpharef{suffix}.txt', 'w') as fd:
            fd.write('{}\n'.format(len(alphas)))
            fd.writelines(['{} {} 1.0\n'.format(i + 1, a)
                           for i, a in enumerate(alphas)])


def read_alpha_file(directory: str) -> List[float]:
    alphas = []
    for suffix in ['', '_empty']:
        fname = f'{directory}/file_alpharef{suffix}.txt'
        if not os.path.isfile(fname):
            break
        with open(fname, 'r') as fd:
            flines = fd.readlines()
            n_orbs = int(flines[0])
            alphas += [float(line.split()[1]) for line in flines[1:n_orbs + 1]]
    return alphas


def read_kpoints_block(calc: Calculator, dct: dict):
    for k, v in dct.items():
        if k in ['kgrid', 'koffset']:
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
        calc.parameters['kpath'] = bandpath(kpath, calc.atoms.cell, npoints=len(kpath) * 10 - 9)
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
    labels = pos_array[:, 0]
    symbols = [''.join([c for c in label if c.isalpha()])
               for label in labels]
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
    calc.atoms.set_array('labels', labels)


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
