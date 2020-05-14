"""

I/O module for python_KI

Written by Edward Linscott Jan 2020

Major modifications
May 2020: converted the module to use CP_calc rather than an ASE calculator

"""

import json
import numpy as np
from ase.io.espresso_cp import Espresso_cp, KEYS, ibrav_to_cell
from ase.atoms import Atoms
from koopmans_cp.utils import warn
import xml.etree.ElementTree as ET
import os


def cpi_diff(calcs):

    # If calcs is a dict, convert it to a list (we only need the values)
    if isinstance(calcs, dict):
        calcs = calcs.values()

    diffs = []

    settings = [c.construct_namelist() for c in calcs]

    blocks = set([b for s in settings for b in s.keys()])
    for block in sorted(blocks):
        keys = set(
            [k for s in settings for k in s.get(block, {}).keys()])
        for key in sorted(keys):
            vals = [s[block].get(key, None) for s in settings]
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
        directory = calc.directory
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


def print_qc(key, value):
    '''
    Prints out a quality control message for testcode to evaluate
    '''
    print(f'<QC> {key} {value}')
