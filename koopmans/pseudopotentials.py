"""

Module containing functions relating to pseudopotentials

Written by Edward Linscott Jan 2020
Moved into _utils Aug 2021
Split into a separate module Sep 2021

"""

import os
from pathlib import Path
import xml.etree.ElementTree as ET
from ase import Atoms
from typing import Dict, Optional


def read_pseudo_file(fd):
    '''

    Reads in settings from a .upf file using XML parser

    '''

    upf = ET.parse(fd).getroot()

    return upf


def set_up_pseudos(calc):

    # Set up pseudopotentials, by...
    #  1. trying to locating the directory as currently specified by the calculator
    #  2. if that fails, checking if $ESPRESSO_PSEUDO is set
    #  3. if that fails, raising an error
    pseudo_dir = calc.parameters.get('pseudo_dir', None)
    if pseudo_dir is None:
        try:
            calc.parameters.pseudo_dir = os.environ.get('ESPRESSO_PSEUDO')
        except KeyError:
            raise NotADirectoryError('Directory for pseudopotentials not found. Please define '
                                     'the environment variable ESPRESSO_PSEUDO or provide a pseudo_dir in '
                                     'the kcp block of your json input file.')
    else:
        if not os.path.isdir(pseudo_dir):
            raise NotADirectoryError(f'The pseudo_dir you provided ({pseudo_dir}) does not exist')


def nelec_from_pseudos(atoms: Atoms, pseudopotentials: Dict[str, str], pseudo_dir: Optional[Path] = None) -> int:
    '''
    Determines the number of electrons in the system using information from pseudopotential files
    '''

    # Works out the pseudo directory (pseudo_dir is given precedence over $ESPRESSO_PSEUDO)
    if pseudo_dir is None:
        if 'ESPRESSO_PSEUDO' in os.environ:
            pseudo_dir = Path(os.environ['ESPRESSO_PSEUDO'])
        else:
            pseudo_dir = Path.cwd()

    valences_dct = {key: read_pseudo_file(pseudo_dir / value).find('PP_HEADER').get(
        'z_valence') for key, value in pseudopotentials.items()}

    if atoms.has('labels'):
        labels = atoms.get_array('labels')
    else:
        labels = atoms.symbols

    valences = [int(float(valences_dct[l])) for l in labels]
    return sum(valences)
