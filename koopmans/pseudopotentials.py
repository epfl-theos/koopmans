"""

Module containing functions relating to pseudopotentials

Written by Edward Linscott Jan 2020
Moved into _utils Aug 2021
Split into a separate module Sep 2021

"""

import os
import re
import json
from itertools import chain
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Optional, List
import xml.etree.ElementTree as ET
from ase import Atoms


@dataclass
class Pseudopotential:
    name: str
    element: str
    path: Path
    functional: str
    library: str
    kind: str
    citations: List[str]
    cutoff_wfc: Optional[float] = None
    cutoff_rho: Optional[float] = None


pseudos_directory = Path(__file__).parents[1] / 'pseudos'

# A database containing all the available pseudopotentials
pseudo_database = []
for pseudo_file in chain(pseudos_directory.rglob('*.UPF'), pseudos_directory.rglob('*.upf')):
    name = pseudo_file.name
    splitname = re.split(r'\.|_|-', name)[0]
    element = splitname[0].upper() + splitname[1:].lower()
    library = pseudo_file.parents[1].name
    functional = pseudo_file.parent.name
    citations = []

    kwargs = {}
    if library.startswith('sssp'):
        [json_name] = list(pseudo_file.parent.glob('*.json'))
        metadata = json.load(open(json_name, 'r'))[element]
        original_library = metadata['pseudopotential'].replace('SG15', 'sg15').replace('Dojo', 'pseudo_dojo')
        if original_library.startswith('sg15') or original_library.startswith('pseudo_dojo'):
            kind = 'norm-conserving'
        elif original_library.startswith('GBRV') or original_library in ['031US', '100US', 'THEOS']:
            kind = 'ultrasoft'
        elif 'PAW' in original_library or original_library == 'Wentzcovitch':
            kind = 'projector-augmented wave'
        else:
            raise ValueError(f'Unrecognised library {original_library}')
        citations += ['Lehaeghere2016', 'Prandini2018']
        for key in ['cutoff_wfc', 'cutoff_rho']:
            kwargs[key] = metadata[key]
    else:
        original_library = library
        if original_library.startswith('sg15') or original_library.startswith('pseudo_dojo'):
            kind = 'norm-conserving'
        else:
            kind = 'unknown'

    if original_library.startswith('sg15'):
        citations.append('Hamann2013')
        citations.append('Schlipf2015')
        if 'relativistic' in original_library:
            citations.append('Scherpelz2016')
    elif original_library.startswith('pseudo_dojo'):
        citations.append('Hamann2013')

    pseudo_database.append(Pseudopotential(name, element, pseudo_file.parent, functional, library, kind, citations))


def pseudos_library_directory(pseudo_library: str, base_functional: str) -> Path:
    return pseudos_directory / pseudo_library / base_functional


def fetch_pseudo(**kwargs):
    matches = [psp for psp in pseudo_database if all([getattr(psp, k) == v for k, v in kwargs.items()])]
    request_str = ', '.join([f'{k} = {v}' for k, v in kwargs.items()])
    if len(matches) == 0:
        raise ValueError('Could not find a pseudopotential in the database matching ' + request_str)
    elif len(matches) > 1:
        raise ValueError('Found multiple pseudopotentials in the database matching ' + request_str)
    else:
        return matches[0]


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
