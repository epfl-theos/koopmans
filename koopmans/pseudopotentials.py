"""

Module containing functions relating to pseudopotentials

Written by Edward Linscott Jan 2020
Moved into _utils Aug 2021
Split into a separate module Sep 2021

"""

import os
import xml.etree.ElementTree as ET


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
