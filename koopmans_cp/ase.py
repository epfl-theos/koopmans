"""

Module for python_KI containing ASE-specific functions
These functions will need tp be replaced if/when we migrate away from ASE

Written by Edward Linscott Jan 2020


"""

import json
import numpy as np
from ase.io.espresso_cp import Espresso_cp, KEYS, ibrav_to_cell
from ase.atoms import Atoms
from koopmans_cp.utils import warn
from koopmans_cp.io import input_dft_from_pseudos, nelec_from_pseudos
import os


def read_json(fd):
    '''

    Reads in settings listed in JSON file

    '''

    if isinstance(fd, str):
        fd = open(fd, 'r')

    bigdct = json.loads(fd.read())

    calc = Espresso_cp()
    calc.parameters['input_data'] = {k: {} for k in KEYS.keys()}

    symbols = []
    positions = []
    cell = None
    labels = None
    settings = {}
    scale_positions = False

    for block, dct in bigdct.items():
        block = block.lower()
        if block == 'calc_param':
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
                    except:
                        settings[k] = v
        elif block == 'atomic_species':
            calc.parameters['pseudopotentials'] = {
                l[0]: l[2] for l in dct['species']}
        elif block == 'atomic_positions':
            pos_array = np.array(dct['positions'])
            labels = pos_array[:, 0]
            symbols = [''.join([c for c in label if c.isalpha()])
                       for label in labels]
            positions = np.array(pos_array[:, 1:], dtype=float)
            if dct.get('units', 'alat') == 'crystal':
                scale_positions = True
        elif block == 'cell_parameters':
            cell = dct['vectors']
        elif block in KEYS:
            for key, value in dct.items():
                if value == "":
                    continue
                try:
                    value = json.loads(value)
                except:
                    pass
                calc.parameters['input_data'][block][key] = value
        else:
            warn(f'The {block} block is not yet implemented and will be ignored')

    # Generating cell if it is missing
    if cell is None:
        _, cell = ibrav_to_cell(calc.parameters['input_data']['system'])

    # Defining our atoms object and tethering it to the calculator
    if scale_positions:
        atoms = Atoms(symbols, scaled_positions=positions,
                      cell=cell, calculator=calc)
    else:
        atoms = Atoms(symbols, positions, cell=cell, calculator=calc)
    atoms.set_array('labels', labels)
    calc.atoms = atoms

    # If they are missing, fill in nelec and input_dft fields using information
    # contained in the pseudopotential files

    if 'pseudopotentials' in calc.parameters:
        if 'input_dft' not in calc.parameters['input_data']['control']:
            calc.parameters['input_data']['system']['input_dft'] = input_dft_from_pseudos(
                calc)
        if 'nelec' not in calc.parameters['input_data']['system']:
            nelec = nelec_from_pseudos(calc)
            calc.parameters['input_data']['system']['nelec'] = nelec
        else:
            nelec = calc.parameters['input_data']['system']['nelec']
        if 'tot_charge' in calc.parameters['input_data']['system']:
            tot_charge = calc.parameters['input_data']['system']['tot_charge']
            calc.parameters['input_data']['system']['nelec'] -= tot_charge
            nelec -= tot_charge
        if 'tot_magnetization' in calc.parameters['input_data']['system']:
            tot_mag = calc.parameters['input_data']['system']['tot_magnetization']
        else:
            tot_mag = nelec%2
            calc.parameters['input_data']['system']['tot_magnetization'] = tot_mag
        if 'nelup' not in calc.parameters['input_data']['system']:
            calc.parameters['input_data']['system']['nelup'] = int(nelec/2 + tot_mag/2)
        if 'neldw' not in calc.parameters['input_data']['system']:
            calc.parameters['input_data']['system']['neldw'] = int(nelec/2 - tot_mag/2)

    return calc, settings


def write_json(fd, calc, calc_param={}):
    '''

    Reads in settings listed in JSON file

    '''

    if isinstance(fd, str):
        fd = open(fd, 'w')

    bigdct = {}
    bigdct['calc_param'] = calc_param

    for key, block in calc.parameters['input_data'].items():
        if len(block) > 0:
            bigdct[key] = dict(block)

    # cell parameters
    if calc.parameters['input_data']['system']['ibrav'] == 0:
        bigdct['cell_parameters'] = {'vectors': [
            list(row) for row in calc.atoms.cell[:]]}

    # atomic positions
    try:
        labels = calc.atoms.get_array('labels')
    except:
        labels = calc.atoms.get_chemical_symbols()
    if calc.parameters['input_data']['system']['ibrav'] == 0:
        bigdct['atomic_positions'] = {'positions': [
            [label] + [str(x) for x in pos] for label, pos in zip(labels, calc.atoms.get_positions())],
            'units': 'alat'}
    else:
        bigdct['atomic_positions'] = {'positions': [
            [label] + [str(x) for x in pos] for label, pos in zip(labels, calc.atoms.get_scaled_positions())],
            'units': 'crystal'}

    # atomic species
    bigdct['atomic_species'] = {'species': [[key, 1.0, val]
                                            for key, val in calc.parameters.pseudopotentials.items()]}

    json.dump(bigdct, fd, indent=2)

    return
