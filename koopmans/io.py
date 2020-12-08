"""

I/O module for python_KI

Written by Edward Linscott Jan 2020

Contains some ASE-specific functions. These functions will need to be
replaced if/when we migrate away from ASE

"""

import os
import copy
import numpy as np
import json as json_ext
import xml.etree.ElementTree as ET
from ase.atoms import Atoms
from ase.units import create_units
from ase.io.espresso import Espresso
from ase.io.espresso import KEYS as pw_keys
from ase.io.espresso_cp import Espresso_cp, ibrav_to_cell
from ase.io.espresso_cp import KEYS as cp_keys
from ase.io.wannier90 import Wannier90
from ase.io.pw2wannier import PW2Wannier
from koopmans import defaults
from koopmans.utils import warn


# Quantum Espresso -- and python_KI -- uses CODATA 2006 internally
ase_units = create_units('2006')


def read_w90_dict(dct, generic_atoms=Atoms()):

    from koopmans.calculators import wannier90

    # Setting up ASE atoms and calc objects
    calc = Wannier90()

    # Loading defaults
    if generic_atoms is None:
        generic_atoms = Atoms()
    symbols = generic_atoms.get_chemical_symbols()
    positions = generic_atoms.get_positions()
    if np.any(generic_atoms.cell.lengths() != 0):
        cell = generic_atoms.cell
    else:
        cell = None
    if generic_atoms.has('labels'):
        labels = generic_atoms.get_array('labels')
    else:
        labels = None
    scale_positions = False

    # Add any generic parameters (calc.atoms.calc points to the generic calculator
    # at the moment, but will be overwritten)
    if generic_atoms.calc is not None:
        for key, val in generic_atoms.calc.parameters.items():
            if key not in calc.parameters and key not in ['pseudopotentials']:
                calc.parameters[key] = val

    dct = {k.lower(): v for k, v in dct.items()}

    for k, v in dct.items():
        if k == "atoms_frac":
            scaled_positions = True
            labels = [line.split()[0] for line in v['positions']]
            symbols = [''.join([c for c in label if c.isalpha()])
                       for label in labels]
            positions = np.array(
                [[float(x) for x in line.split()[1:]] for line in v['positions']])
        elif k == "unit_cell_cart":
            if v['units'].lower() != 'ang':
                raise NotImplementedError(
                    'unit_cell_cart with units != "ang" is not yet implemented')
            cell = np.array(v['vectors'])
        elif k == 'mp_grid':
            calc.parameters['koffset'] = [0, 0, 0]
            calc.parameters['kpts'] = v
        else:
            try:
                v = json_ext.loads(v)
            except (TypeError, json_ext.decoder.JSONDecodeError) as e:
                pass
            calc.parameters[k] = v

    # Defining our atoms object and tethering it to the calculator
    if scale_positions:
        calc.atoms = Atoms(symbols, scaled_positions=positions,
                           cell=cell, calculator=calc)
    else:
        calc.atoms = Atoms(symbols, positions, cell=cell, calculator=calc)
    calc.atoms.set_array('labels', labels)

    # Return python_KI-type calculator object rather than ASE calculator
    return wannier90.W90_calc(calc)


def read_pw2wannier_dict(dct, generic_atoms=Atoms()):

    from koopmans.calculators import pw2wannier

    # Setting up ASE atoms and calc objects
    calc = PW2Wannier()

    calc.parameters['inputpp'] = {}
    for k, v in dct.items():
        try:
            v = json_ext.loads(v)
        except (TypeError, json_ext.decoder.JSONDecodeError) as e:
            pass
        calc.parameters['inputpp'][k] = v

    # Attaching an atoms object (though pw2wannier doesn't need this information,
    # ASE will complain if it's missing)
    calc.atoms = copy.deepcopy(generic_atoms)
    calc.atoms.calc = calc

    # Return python_KI-type calculator object rather than ASE calculator
    return pw2wannier.PW2Wannier_calc(calc)


def read_qe_generic_dict(dct, calc, generic_atoms=Atoms()):
    '''

    Reads in dict of cp/pw input file and adds the settings to the provided calc object

    '''

    if isinstance(calc, Espresso_cp):
        qe_keys = cp_keys
    elif isinstance(calc, Espresso):
        qe_keys = pw_keys
    else:
        raise TypeError(f'io.read_qe_generic_dict() does not accept "calc" with class {calc.__class__}')

    calc.parameters['input_data'] = {k: {} for k in qe_keys.keys()}

    if generic_atoms is None:
        generic_atoms = Atoms()

    symbols = generic_atoms.get_chemical_symbols()
    positions = generic_atoms.get_positions()
    if np.any(generic_atoms.cell.lengths() != 0):
        cell = generic_atoms.cell
    else:
        cell = None
    if generic_atoms.has('labels'):
        labels = generic_atoms.get_array('labels')
    else:
        labels = None
    scale_positions = False
    skipped_blocks = []

    # Add any generic parameters (calc.atoms.calc points to the generic calculator
    # at the moment, but will be overwritten)
    if generic_atoms.calc is not None:
        for key, val in generic_atoms.calc.parameters.items():
            if key not in calc.parameters:
                calc.parameters[key] = val

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
            units = subdct.get('units', 'angstrom').lower()
            if units == 'angstrom':
                pass
            elif units == 'crystal':
                scale_positions = True
            else:
                raise NotImplementedError(
                    f'atomic_positions units = {units} is not yet implemented')

        elif block == 'cell_parameters':
            cell = subdct.get('vectors', None)
            units = subdct.get('units', None)
            if cell is None and units in [None, 'alat']:
                if 'ibrav' not in subdct:
                    raise KeyError('"ibrav" (and any other cell-related '
                                   'parameters) should be listed in the '
                                   '"cell_parameters" block when not '
                                   'explicitly providing the cell')
                # Copy over all the rest of the information to the system block
                for k, v in subdct.items():
                    if k not in ['vectors', 'units']:
                        calc.parameters['input_data']['system'][k] = v
            elif cell is not None and units == 'angstrom':
                pass
            else:
                raise NotImplementedError('the combination of vectors, ibrav, & units '
                                          'in the cell_parameter block cannot be read (may not yet be '
                                          'implemented)')
        elif block in qe_keys:
            for key, value in subdct.items():
                if value == "":
                    continue

                # Force pseudo_dir to be an absolute path
                if key == 'pseudo_dir' and value[0] != '/':
                    descend_depth = value.count('../')
                    splitdir = os.getcwd().rsplit('/', descend_depth)
                    value = splitdir[0] + '/' + value.lstrip('../')
                    if not os.path.isdir(value):
                        raise ValueError(
                            f'Could not parse pseudo_dir; {value} is not a directory')

                try:
                    value = json_ext.loads(value)
                except (TypeError, json_ext.decoder.JSONDecodeError) as e:
                    pass
                calc.parameters['input_data'][block][key] = value
        else:
            skipped_blocks.append(block)

    # Generating cell if it is missing
    if cell is None:
        _, cell = ibrav_to_cell(calc.parameters['input_data']['system'])

    # Defining our atoms object and tethering it to the calculator
    if scale_positions:
        calc.atoms = Atoms(symbols, scaled_positions=positions,
                           cell=cell, calculator=calc)
    else:
        calc.atoms = Atoms(symbols, positions, cell=cell, calculator=calc)
    calc.atoms.set_array('labels', labels)

    return skipped_blocks


def read_cp_dict(dct, generic_atoms=None):
    '''

    Reads in dict of cp input file and returns a CP_calc object

    '''

    from koopmans.calculators import cp

    calc = Espresso_cp()
    skipped_blocks = read_qe_generic_dict(dct, calc, generic_atoms)

    set_up_pseudos(calc)

    for block in skipped_blocks:
        warn(f'The {block} block is not yet implemented and will be ignored')

    # If they are missing, fill the nelec/nelup/neldw field using information contained
    # in the pseudopotential files
    if 'pseudopotentials' in calc.parameters:
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
            tot_mag = nelec % 2
            calc.parameters['input_data']['system']['tot_magnetization'] = tot_mag
        if 'nelup' not in calc.parameters['input_data']['system']:
            calc.parameters['input_data']['system']['nelup'] = int(
                nelec / 2 + tot_mag / 2)
        if 'neldw' not in calc.parameters['input_data']['system']:
            calc.parameters['input_data']['system']['neldw'] = int(
                nelec / 2 - tot_mag / 2)

    return cp.CP_calc(calc)


def read_pw_dict(dct, generic_atoms=None):
    '''

    Reads in dict of pw input file and returns a PW_calc object

    '''

    from koopmans.calculators import pw

    calc = Espresso()
    skipped_blocks = read_qe_generic_dict(dct, calc, generic_atoms)

    set_up_pseudos(calc)

    for block in skipped_blocks:
        if block == 'k_points':
            read_kpoints_block(calc, dct[block])
        else:
            warn(f'The {block} block is not yet implemented and will be ignored')

    return pw.PW_calc(calc)


def read_atoms_dct(dct):
    '''

    Reads the atoms block

    '''

    # Take advantage of the PW dict reader in order to read the block settings
    calc = Espresso()
    skipped_blocks = read_qe_generic_dict(dct, calc)
    if 'k_points' in skipped_blocks:
        read_kpoints_block(calc, dct['k_points'])

    atoms = calc.atoms

    # Remove the code-specific settings
    del atoms.calc.parameters['input_data']

    return atoms


def read_generic_dict(dct):
    '''

    Reads in dict of generic settings and returns a tidied dict

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
                settings[k] = json_ext.loads(v)
            except (TypeError, json_ext.decoder.JSONDecodeError) as e:
                settings[k] = v
    return settings


def read_json(fd):
    '''

    Reads in settings listed in JSON file

    '''

    from koopmans.workflows.singlepoint import SinglepointWorkflow
    from koopmans.workflows.convergence import ConvergenceWorkflow
    from koopmans.workflows.pbe_dscf_with_pw import DeltaSCFWorkflow

    if isinstance(fd, str):
        fd = open(fd, 'r')

    bigdct = json_ext.loads(fd.read())

    readers = {'cp': read_cp_dict, 'w90_occ': read_w90_dict, 'w90_emp': read_w90_dict,
               'pw2wannier': read_pw2wannier_dict, 'pw': read_pw_dict}
    calcs_dct = {}

    # Load default values
    generic_atoms = None
    if 'atoms' in bigdct:
        generic_atoms = read_atoms_dct(bigdct['atoms'])
        del bigdct['atoms']

    if 'w90' in bigdct:
        bigdct['w90_occ'] = bigdct['w90']['occ']
        bigdct['w90_emp'] = bigdct['w90']['emp']
        del bigdct['w90']

    for block, dct in bigdct.items():
        block = block.lower()
        if block == 'workflow':
            workflow_settings = read_generic_dict(dct)
        elif block in readers:
            # Read in the block
            calcs_dct[block] = readers[block](dct, generic_atoms)
        else:
            raise ValueError(f'Unrecognised block "{block}" in json input file; '
                             'valid options are workflow/' + '/'.join(readers.keys()))

    task_name = workflow_settings.pop('task', 'singlepoint')
    if task_name == 'singlepoint':
        workflow = SinglepointWorkflow(workflow_settings, calcs_dct)
    elif task_name == 'convergence':
        workflow = ConvergenceWorkflow(workflow_settings, calcs_dct)
    elif task_name == 'environ_dscf':
        workflow = DeltaSCFWorkflow(workflow_settings, calcs_dct)
    else:
        raise ValueError('Invalid task name "{task_name}"')

    return workflow


def write_json(fd, calcs=[], workflow_settings={}):
    '''

    Writes out settings to a JSON file

    '''

    from koopmans.calculators import cp, pw

    if isinstance(fd, str):
        fd = open(fd, 'w')

    if not isinstance(calcs, list):
        calcs = [calcs]

    bigdct = {}
    bigdct['workflow'] = {k: v for k,
                          v in workflow_settings.items() if v is not None}

    for calc in calcs:
        if isinstance(calc, (cp.CP_calc, pw.PW_calc)):
            if isinstance(calc, cp.CP_calc):
                code = 'cp'
            else:
                code = 'pw'
            bigdct[code] = {}
            # Update the ase calculator
            calc = calc.calc
            for key, block in calc.parameters['input_data'].items():
                if len(block) > 0:
                    bigdct[code][key] = {k: v for k, v in dict(
                        block).items() if v is not None}

            # cell parameters
            if calc.parameters['input_data']['system'].get('ibrav', None) == 0:
                bigdct[code]['cell_parameters'] = {'vectors': [
                    list(row) for row in calc.atoms.cell[:]], 'units': 'angstrom'}

            # atomic positions
            if calc.atoms.has('labels'):
                labels = calc.atoms.get_array('labels')
            else:
                labels = calc.atoms.get_chemical_symbols()
            if calc.parameters['input_data']['system'].get('ibrav', None) == 0:
                bigdct[code]['atomic_positions'] = {'positions': [
                    [label] + [str(x) for x in pos] for label, pos in zip(labels, calc.atoms.get_positions())],
                    'units': 'angstrom'}
            else:
                bigdct[code]['atomic_positions'] = {'positions': [
                    [label] + [str(x) for x in pos] for label, pos in zip(labels, calc.atoms.get_scaled_positions())],
                    'units': 'alat'}

            # atomic species
            bigdct[code]['atomic_species'] = {'species': [[key, 1.0, val]
                                                          for key, val in calc.parameters.pseudopotentials.items()]}

            # PW-specific blocks
            if code == 'pw':
                if 'kpts' in calc.parameters:
                    kpts = {'kind': 'automatic',
                            'kpts': calc.parameters['kpts'],
                            'koffset': calc.parameters['koffset']}
                else:
                    kpts = {'kind': 'gamma'}
                bigdct[code]['k_points'] = kpts
        else:
            raise ValueError(
                f'Writing of {calc.__class__} with write_json is not yet implemented')

    json_ext.dump(bigdct, fd, indent=2)

    return


def parse_physical(value):
    '''
    Takes in a value that potentially has a unit following a float,
    converts the value to ASE's default units (Ang, eV), and returns
    that value
    '''

    if isinstance(value, float):
        return value
    elif isinstance(value, int):
        return float(value)
    elif isinstance(value, str):
        splitline = value.strip().split()
        if len(splitline) == 1:
            return float(splitline[0])
        elif len(splitline) == 2:
            [value, units] = splitline
            value = float(value)
            matching_units = [
                u for u in ase_units if u.lower() == units.lower()]
            if len(matching_units) == 1:
                return value * ase_units[matching_units[0]]
            elif len(matching_units) > 1:
                raise ValueError(
                    f'Multiple matches for {units} found; this should not happen')
            else:
                raise NotImplementedError(
                    f'{units} not implemented in koopmans.io.parse_physical')


def read_pseudo_file(fd):
    '''

    Reads in settings from a .upf file using XML parser

    '''

    upf = ET.parse(fd).getroot()

    return upf


def get_pseudo_dir(calc):
    '''
    Works out the pseudo directory (pseudo_dir is given precedence over $ESPRESSO_PSEUDO)
    '''

    pseudo_dir = None
    if 'control' in calc.parameters['input_data']:
        if 'pseudo_dir' in calc.parameters['input_data']['control']:
            pseudo_dir = calc.parameters['input_data']['control']['pseudo_dir']
    if pseudo_dir is not None:
        return pseudo_dir
    elif 'ESPRESSO_PSEUDO' in os.environ:
        return os.environ['ESPRESSO_PSEUDO']
    else:
        return './'


def set_up_pseudos(calc):

    # Set up pseudopotentials, by...
    #  1. trying to locating the directory as currently specified by the calculator
    #  2. if that fails, checking if $ESPRESSO_PSEUDO is set
    #  3. if that fails, raising an error
    pseudo_dir = calc.parameters['input_data']['control'].get(
        'pseudo_dir', None)
    if pseudo_dir is None or not os.path.isdir(pseudo_dir):
        try:
            calc.parameters['input_data']['control']['pseudo_dir'] = os.environ.get(
                'ESPRESSO_PSEUDO')
        except KeyError:
            raise NotADirectoryError('Directory for pseudopotentials not found. Please define '
                                     'the environment variable ESPRESSO_PSEUDO or provide a pseudo_dir in '
                                     'the cp block of your json input file.')

    # If it is missing, fill the input_dft field using information contained
    # in the pseudopotential files
    if 'pseudopotentials' in calc.parameters:
        if 'input_dft' not in calc.parameters['input_data']['control']:
            calc.parameters['input_data']['system']['input_dft'] = input_dft_from_pseudos(
                calc)


def nelec_from_pseudos(calc):
    '''
    Determines the number of electrons in the system using information from pseudopotential files
    '''

    directory = get_pseudo_dir(calc)
    valences_dct = {key: read_pseudo_file(directory + value).find('PP_HEADER').get(
        'z_valence') for key, value in calc.parameters['pseudopotentials'].items()}
    valences = [int(float(valences_dct[l]))
                for l in calc.atoms.get_array('labels')]
    return sum(valences)


def input_dft_from_pseudos(calc):
    '''
    Determines input_dft using information from pseudopotential files
    '''

    directory = get_pseudo_dir(calc)
    input_dft = list(set([read_pseudo_file(directory + '/' + fname).find('PP_HEADER').get('functional')
                          for fname in calc.parameters['pseudopotentials'].values()]))

    if len(input_dft) != 1:
        warn('The listed pseudopotentials do not use the same functional; they are using '
             ', '.join(input_dft))

    return input_dft[0]


def write_alpha_file(directory, alphas, filling):
    a_filled = [a for a, f in zip(alphas, filling) if f]
    a_empty = [a for a, f in zip(alphas, filling) if not f]
    for alphas, suffix in zip([a_filled, a_empty], ['', '_empty']):
        with open(f'{directory}/file_alpharef{suffix}.txt', 'w') as fd:
            fd.write('{}\n'.format(len(alphas)))
            fd.writelines(['{} {} 1.0\n'.format(i + 1, a)
                           for i, a in enumerate(alphas)])


def read_alpha_file(directory):
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


def read_kpoints_block(calc, dct):
    if dct['kind'] == 'gamma':
        kpts = [(1, 1, 1)]
        koffset = 0
    elif dct['kind'] == 'automatic':
        kpts = dct['kpts']
        koffset = dct['koffset']
    calc.parameters['kpts'] = kpts
    calc.parameters['koffset'] = koffset
    return
