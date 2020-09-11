"""

Module for python_KI containing ASE-specific functions
These functions will need tp be replaced if/when we migrate away from ASE

Written by Edward Linscott Jan 2020


"""

import json
import numpy as np
from ase import Atoms
from ase.io.espresso_cp import Espresso_cp, KEYS, ibrav_to_cell
from ase.io.wannier90 import Wannier90
from ase.io.pw2wannier import PW2Wannier
from ase.atoms import Atoms
from koopmans_cp.utils import warn
from koopmans_cp.io import input_dft_from_pseudos, nelec_from_pseudos
from koopmans_cp.calculators import cp, wannier90, pw2wannier
import os

def read_w90_dict(dct):
    # Setting up ASE atoms and calc objects
    calc = Wannier90()
    calc.atoms = Atoms()
    calc.atoms.calc = calc

    for k, v in dct.items():
        calc.parameters[k] = v

    # Return python_KI-type calculator object rather than ASE calculator
    return wannier90.W90_calc(calc)

def read_pw2wannier_dict(dct):
    # Setting up ASE atoms and calc objects
    calc = PW2Wannier()
    calc.atoms = Atoms()
    calc.atoms.calc = calc

    calc.parameters['inputpp'] = {}
    for k, v in dct.items():
        try:
            v = json.loads(v)
        except:
            pass
        calc.parameters['inputpp'][k] = v
    
    # Return python_KI-type calculator object rather than ASE calculator
    return pw2wannier.PW2Wannier_calc(calc)

def read_cp_dict(dct):
    '''

    Reads in dict of cp input file and returns a CP_calc object

    '''
    calc = Espresso_cp()
    calc.parameters['input_data'] = {k: {} for k in KEYS.keys()}

    symbols = []
    positions = []
    cell = None
    labels = None
    scale_positions = False

    for block, subdct in dct.items():
        if block == 'atomic_species':
            calc.parameters['pseudopotentials'] = {
                l[0]: l[2] for l in subdct['species']}
        elif block == 'atomic_positions':
            pos_array = np.array(subdct['positions'])
            labels = pos_array[:, 0]
            symbols = [''.join([c for c in label if c.isalpha()])
                       for label in labels]
            positions = np.array(pos_array[:, 1:], dtype=float)
            if subdct.get('units', 'alat') == 'crystal':
                scale_positions = True
        elif block == 'cell_parameters':
            cell = subdct['vectors']
        elif block in KEYS:
            for key, value in subdct.items():
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

    return cp.CP_calc(calc)

def read_workflow_dict(dct):
    '''

    Reads in dict of workflow settings and returns a tidied dict

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
            except:
                settings[k] = v
    return settings

def read_json(fd):
    '''

    Reads in settings listed in JSON file

    '''

    if isinstance(fd, str):
        fd = open(fd, 'r')

    bigdct = json.loads(fd.read())

    readers = {'cp': read_cp_dict, 'w90': read_w90_dict, 'pw2wannier': read_pw2wannier_dict}
    calcs = {}
    for block, dct in bigdct.items():
        block = block.lower()
        if block == 'workflow':
            workflow_settings = read_workflow_dict(dct)
        elif block in readers:
            calcs[block] = readers[block](dct)
        else:
            raise ValueError(f'Unrecognised block "{block}" in json input file; '
                            'valid options are workflow/' + '/'.join(readers.keys()))

    return workflow_settings, calcs


def write_json(fd, calcs=[], workflow_settings={}):
    '''

    Writes out settings to a JSON file

    '''

    if isinstance(fd, str):
        fd = open(fd, 'w')

    if not isinstance(calcs, list):
        calcs = [calcs]

    bigdct = {}
    bigdct['workflow_settings'] = workflow_settings

    for calc in calcs:
        if isinstance(calc, CP_calc):
            bigdct['cp'] = {}
            calc = calc._ase_calc
            for key, block in calc.parameters['input_data'].items():
                if len(block) > 0:
                    bigdct['cp'][key] = dict(block)

            # cell parameters
            if calc.parameters['input_data']['system']['ibrav'] == 0:
                bigdct['cp']['cell_parameters'] = {'vectors': [
                    list(row) for row in calc.atoms.cell[:]]}

            # atomic positions
            try:
                labels = calc.atoms.get_array('labels')
            except:
                labels = calc.atoms.get_chemical_symbols()
            if calc.parameters['input_data']['system']['ibrav'] == 0:
                bigdct['cp']['atomic_positions'] = {'positions': [
                    [label] + [str(x) for x in pos] for label, pos in zip(labels, calc.atoms.get_positions())],
                    'units': 'alat'}
            else:
                bigdct['cp']['atomic_positions'] = {'positions': [
                    [label] + [str(x) for x in pos] for label, pos in zip(labels, calc.atoms.get_scaled_positions())],
                    'units': 'crystal'}

            # atomic species
            bigdct['cp']['atomic_species'] = {'species': [[key, 1.0, val]
                                                for key, val in calc.parameters.pseudopotentials.items()]}
        else:
            raise ValueError(f'Writing of {calc.__class__} with write_json is not yet implemented')

    json.dump(bigdct, fd, indent=2)

    return
