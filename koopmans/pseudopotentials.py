"""

Module containing functions relating to pseudopotentials

Written by Edward Linscott Jan 2020
Moved into _utils Aug 2021
Split into a separate module Sep 2021

"""

import os
import re
from pathlib import Path
from typing import Dict, Optional, Union
import xml.etree.ElementTree as ET
from ase import Atoms


pseudos_directory = Path(__file__).parents[1] / 'pseudos'
available_pseudo_libraries = {func.lower(): [directory.name for directory in pseudos_directory.rglob(
    '*') if (directory / func).exists()] for func in ['LDA', 'PBE', 'PBEsol']}


def fetch_pseudo_from_library(element: str, pseudo_library: str, base_functional: str):
    '''
    Fetches the appropriate pseudopotential for a particular element from a particular pseudopotential library
    corresponding to a particular base functional
    '''

    directory = pseudos_library_directory(pseudo_library, base_functional)
    pseudo_matches = [psp.name for psp in directory.glob(
        '*') if re.split(r'\.|_|-', psp.name)[0].lower() == element.lower()]
    if len(pseudo_matches) == 0:
        raise ValueError(f'Failed to find a pseudopotential corresponding to {element} in {directory}')
    else:
        return sorted(pseudo_matches)[-1]


def pseudos_library_directory(pseudo_library: str, base_functional: str) -> Path:
    directory = pseudos_directory / pseudo_library
    if directory.is_symlink():
        directory = Path(os.path.realpath(directory))
    return directory / base_functional


def read_pseudo_file(fd):
    '''

    Reads in settings from a .upf file using XML parser

    '''

    upf = ET.parse(fd).getroot()

    return upf


def nelec_from_pseudos(atoms: Atoms, pseudopotentials: Dict[str, str],
                       pseudo_dir: Path) -> int:
    '''
    Determines the number of electrons in the system using information from pseudopotential files
    '''

    valences_dct = {key: read_pseudo_file(pseudo_dir / value).find('PP_HEADER').get(
        'z_valence') for key, value in pseudopotentials.items()}

    if atoms.has('labels'):
        labels = atoms.get_array('labels')
    else:
        labels = atoms.symbols

    valences = [int(float(valences_dct[l])) for l in labels]
    return sum(valences)
