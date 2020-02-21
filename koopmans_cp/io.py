import json
import numpy as np
import warnings
from ase.io.espresso_cp import Espresso_cp, KEYS, ibrav_to_cell
from ase.atoms import Atoms
import xml.etree.ElementTree as ET

"""

I/O module for python_KI

Written by Edward Linscott Jan 2020

"""

import os


def cpi_diff(calcs):

    # If calcs is a dict, convert it to a list (we only need the values)
    if isinstance(calcs, dict):
        calcs = calcs.values()

    diffs = []

    blocks = set([b for c in calcs for b in c.parameters['input_data'].keys()])
    for block in sorted(blocks):
        keys = set(
            [k for c in calcs for k in c.parameters['input_data'].get(block, {}).keys()])
        for key in sorted(keys):
            vals = [c.parameters['input_data']
                    [block].get(key, None) for c in calcs]
            if len(set(vals)) > 1:
                print(f'{block}.{key}: ' + ', '.join(map(str, vals)))
                diffs.append(key)

    return diffs


def print_summary(alpha_df, error_df):
    # Printing out a progress summary
    print('\nalpha')
    print(alpha_df)
    print('\ndelta E - lambda^alpha_ii (eV)')
    print(error_df)
    print('')


def write_alpharef(alphas, filling, directory='.', duplicate=True):
    '''
    Generates file_alpharef.txt and file_alpharef_empty.txt from a list of alpha values

    Arguments:
       alphas    -- a list of alpha values (floats)
       filling   -- a list of booleans; true if the corresponding orbital is filled
       directory -- the directory within which to write the files
       duplicate -- if True, use the same result for spin-up and spin-down
    '''

    a_filled = [a for a, f in zip(alphas, filling) if f]
    a_empty = [a for a, f in zip(alphas, filling) if not f]
    for alphas, suffix in zip([a_filled, a_empty], ['', '_empty']):
        if duplicate:
            alphas += alphas
        with open(f'{directory}/file_alpharef{suffix}.txt', 'w') as fd:
            fd.write('{}\n'.format(len(alphas)))
            fd.writelines(['{} {} 1.0\n'.format(i+1, a)
                           for i, a in enumerate(alphas)])


def read_alpharef(calc=None, directory=None):
    '''
    Reads in file_alpharef.txt and file_alpharef_empty.txt from a calculation's directory

    Arguments:
       calc      -- an ASE calculator
       directory -- a directory

    Output:
       alphas -- a list of alpha values (1 per orbital)
    '''

    if calc is not None:
        directory = calc.label.rsplit('/', 1)[0]
    elif directory is None:
        raise ValueError(
            'read_alpharef called without a calculator or a directory. Please provide at least one.')

    alphas = []
    for suffix in ['', '_empty']:
        fname = f'{directory}/file_alpharef{suffix}.txt'
        if not os.path.isfile(fname):
            break
        with open(fname, 'r') as fd:
            flines = fd.readlines()
            n_orbs = int(flines[0]) // 2
            alphas += [float(line.split()[1]) for line in flines[1:n_orbs + 1]]
    return alphas


def parse_algebraic_expression(expr, calc):
    if not isinstance(expr, str):
        return expr
    if all([c.isalpha() for c in expr]):
        return expr

    expr = expr.replace('/', ' / ').replace('*', ' * ').split()
    for i, term in enumerate(expr):
        if term in ['*', '/']:
            continue
        elif all([c.isalpha() for c in term]):
            if getattr(calc, term, None) is None:
                raise ValueError(f'Failed to parse ' + ''.join(expr))
            else:
                expr[i] = getattr(calc, term)
        else:
            expr[i] = float(term)

    value = float(expr[0])
    for op, term in zip(expr[1::2], expr[2::2]):
        if op == '*':
            value *= float(term)
        elif op == '/':
            value /= float(term)
        else:
            raise ValueError('Failed to parse ' +
                             ''.join([str(e) for e in expr]))

    return value


def parse_algebraic_expressions(calc):

    cp_settings = {key: value for block in calc.parameters['input_data'].values()
                   for key, value in block.items()}

    for key, value in cp_settings.items():
        if key in ['pseudo_dir', 'outdir']:
            continue
        setattr(calc, key, parse_algebraic_expression(value, calc))

    return calc


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
            calc.parameters['input_data']['system']['nelec'] = nelec_from_pseudos(
                calc)
        if 'nelup' not in calc.parameters['input_data']['system']:
            calc.parameters['input_data']['system']['nelup'] = calc.parameters['input_data']['system']['nelec']//2
        if 'neldw' not in calc.parameters['input_data']['system']:
            calc.parameters['input_data']['system']['neldw'] = calc.parameters['input_data']['system']['nelec'] - \
                calc.parameters['input_data']['system']['nelup']

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


def read_pseudopotential(fd):
    '''

    Reads in settings from a .upf file using XML parser

    '''

    upf = ET.parse(fd).getroot()

    return upf


def get_pseudo_dir(calc):
    '''
    Works out the pseudo directory (pseudo_dir is given precedence over $ESPRESSO_PSEUDO)
    '''

    if getattr(calc, 'pseudo_dir', None) is not None:
        return calc.pseudo_dir
    elif 'ESPRESSO_PSEUDO' in os.environ:
        return os.environ['ESPRESSO_PSEUDO']
    else:
        return './'


def nelec_from_pseudos(calc):
    '''
    Determines the number of electrons in the system using information from pseudopotential files
    '''

    directory = get_pseudo_dir(calc)
    valences_dct = {key: read_pseudopotential(directory + value).find('PP_HEADER').get(
        'z_valence') for key, value in calc.parameters['pseudopotentials'].items()}
    valences = [int(float(valences_dct[l]))
                for l in calc.atoms.get_array('labels')]
    return sum(valences)


def input_dft_from_pseudos(calc):
    '''
    Determines input_dft using information from pseudopotential files
    '''

    directory = get_pseudo_dir(calc)
    input_dft = list(set([read_pseudopotential(directory + fname).find('PP_HEADER').get('functional')
                          for fname in calc.parameters['pseudopotentials'].values()]))

    if len(input_dft) != 1:
        warn('The listed pseudopotentials do not use the same functional; they are using ' +
             ', '.join(input_dft))
    return input_dft[0]


def _warning(message, category=UserWarning, filename='', lineno=-1, file=None, line=None):
    '''
    Monkey-patching warnings.warn
    '''
    print(f'{category.__name__}: {message}')


warnings.showwarning = _warning


def warn(message):
    '''
    Allowing the monkey-patched warnings.warn to be imported as io.warn
    '''
    warnings.warn(message)


def print_qc(key, value):
    '''
    Prints out a quality control message for testcode to evaluate
    '''
    print(f'<QC> {key} {value}')
