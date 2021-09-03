"""

Generic I/O functions for koopmans

Written by Edward Linscott Jan 2020
Moved into _utils Aug 2021

"""

import os
import glob
import numpy as np
import xml.etree.ElementTree as ET
from ase.atoms import Atoms
from ase.dft.kpoints import bandpath, BandPath


def construct_cell_parameters_block(calc):
    return {'vectors': [list(row) for row in calc.atoms.cell[:]], 'units': 'angstrom'}


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
        return '.'


def set_up_pseudos(calc):

    # Set up pseudopotentials, by...
    #  1. trying to locating the directory as currently specified by the calculator
    #  2. if that fails, checking if $ESPRESSO_PSEUDO is set
    #  3. if that fails, raising an error
    pseudo_dir = calc.parameters['input_data']['control'].get('pseudo_dir', None)
    if pseudo_dir is None or not os.path.isdir(pseudo_dir):
        try:
            calc.parameters['input_data']['control']['pseudo_dir'] = os.environ.get('ESPRESSO_PSEUDO')
        except KeyError:
            raise NotADirectoryError('Directory for pseudopotentials not found. Please define '
                                     'the environment variable ESPRESSO_PSEUDO or provide a pseudo_dir in '
                                     'the kcp block of your json input file.')


def nelec_from_pseudos(calc):
    '''
    Determines the number of electrons in the system using information from pseudopotential files
    '''

    directory = get_pseudo_dir(calc)
    valences_dct = {key: read_pseudo_file(directory + value).find('PP_HEADER').get(
        'z_valence') for key, value in calc.parameters['pseudopotentials'].items()}
    if calc.atoms.has('labels'):
        labels = calc.atoms.get_array('labels')
    else:
        labels = calc.atoms.symbols
    valences = [int(float(valences_dct[l])) for l in labels]
    return sum(valences)


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
        kpts = None
        koffset = None
    elif dct['kind'] == 'automatic':
        kpts = dct['kpts']
        koffset = dct['koffset']
    calc.parameters['kpts'] = kpts
    calc.parameters['koffset'] = koffset

    if 'kpath' in dct:
        read_kpath(calc, dct['kpath'])

    return


def read_kpath(calc, kpath):
    calc.atoms.cell.pbc = True
    if isinstance(kpath, str):
        # Interpret kpath as a string of points in the BZ
        calc.parameters['kpath'] = bandpath(kpath, calc.atoms.cell, npoints=len(kpath) * 10 - 9)
    else:
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


def read_atomic_species(calc, dct):
    calc.parameters['pseudopotentials'] = {l[0]: l[2] for l in dct['species']}


def read_atomic_positions(calc, dct):

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


def read_cell_parameters(calc, dct):
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


def indented_print(text='', indent=0, **kwargs):
    global print_call_end
    for substring in text.split('\n'):
        if print_call_end == '\n':
            print(' ' * indent + substring, **kwargs)
        else:
            print(substring, **kwargs)
    print_call_end = kwargs.get('end', '\n')
