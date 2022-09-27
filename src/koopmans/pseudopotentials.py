"""

Module containing functions relating to pseudopotentials

Written by Edward Linscott Jan 2020
Moved into _utils Aug 2021
Split into a separate module Sep 2021

"""

import json
import os
import re
from dataclasses import dataclass
from itertools import chain
from pathlib import Path
from typing import Any, Dict, List, Optional

from ase import Atoms
from upf_to_json import upf_to_json


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


pseudos_directory = Path(__file__).parent / 'pseudopotentials'

# A database containing all the available pseudopotentials
pseudo_database: List[Pseudopotential] = []
for pseudo_file in chain(pseudos_directory.rglob('*.UPF'), pseudos_directory.rglob('*.upf')):
    name = pseudo_file.name
    splitname = re.split(r'\.|_|-', name)[0]
    element = splitname[0].upper() + splitname[1:].lower()
    library = pseudo_file.parents[1].name
    functional = pseudo_file.parent.name
    citations: List[str] = []

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
            raise ValueError(f'Unrecognized library {original_library}')
        citations += ['Lejaeghere2016', 'Prandini2018']
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
        citations.append('vanSetten2018')

    pseudo_database.append(Pseudopotential(name, element, pseudo_file.parent,
                                           functional, library, kind, citations, **kwargs))


def pseudos_library_directory(pseudo_library: str, base_functional: str) -> Path:
    return pseudos_directory / pseudo_library / base_functional


def fetch_pseudo(**kwargs: Any) -> Pseudopotential:
    matches = [psp for psp in pseudo_database if all([getattr(psp, k) == v for k, v in kwargs.items()])]
    request_str = ', '.join([f'{k} = {v}' for k, v in kwargs.items()])
    if len(matches) == 0:
        raise ValueError('Could not find a pseudopotential in the database matching ' + request_str)
    elif len(matches) > 1:
        raise ValueError('Found multiple pseudopotentials in the database matching ' + request_str)
    else:
        return matches[0]


def read_pseudo_file(filename: Path) -> Dict[str, Any]:
    '''

    Reads in settings from a .upf file

    '''

    upf: Dict[str, Any] = upf_to_json(open(filename, 'r').read(), filename.name)

    return upf['pseudo_potential']


def valence_from_pseudo(filename: str, pseudo_dir: Optional[Path] = None) -> int:
    '''
    Determines the valence of a pseudopotential
    '''

    # Works out the pseudo directory (pseudo_dir is given precedence over $ESPRESSO_PSEUDO)
    if pseudo_dir is None:
        if 'ESPRESSO_PSEUDO' in os.environ:
            pseudo_dir = Path(os.environ['ESPRESSO_PSEUDO'])
        else:
            pseudo_dir = Path.cwd()
    elif isinstance(pseudo_dir, str):
        pseudo_dir = Path(pseudo_dir)

    return int(read_pseudo_file(pseudo_dir / filename)['header']['z_valence'])


def nelec_from_pseudos(atoms: Atoms, pseudopotentials: Dict[str, str],
                       pseudo_dir: Optional[Path] = None) -> int:
    '''
    Determines the number of electrons in the system using information from pseudopotential files
    '''

    valences_dct = {key: valence_from_pseudo(value, pseudo_dir) for key, value in pseudopotentials.items()}

    if len(set(atoms.get_tags())) > 1:
        labels = [s + str(t) if t > 0 else s for s, t in zip(atoms.symbols, atoms.get_tags())]
    else:
        labels = atoms.symbols

    valences = [valences_dct[l] for l in labels]
    return sum(valences)


def expected_subshells(atoms: Atoms, pseudopotentials: Dict[str, str],
                       pseudo_dir: Optional[Path] = None) -> Dict[str, List[str]]:
    """
    Determine which subshells will make up the valences of a set of pseudopotentials.

    Returns
    -------
    Dict[str, List[str]]
        a dict mapping element names to a corresponding list of suborbitals that *might* be in the pseudopotential
        valence (depending on how many bands are included)

    """

    z_core_to_first_orbital = {0: '1s', 2: '2s', 4: '2p', 10: '3s', 12: '3p', 18: '3d', 28: '4s', 30: '4p',
                               36: '4d', 46: '4f', 60: '5s', 62: '5p', 68: '6s'}

    expected_orbitals = {}
    for atom in atoms:
        if atom.symbol in expected_orbitals:
            continue
        pseudo_file = pseudopotentials[atom.symbol]
        z_core = atom.number - valence_from_pseudo(pseudo_file, pseudo_dir)
        if z_core in z_core_to_first_orbital:
            first_orbital = z_core_to_first_orbital[z_core]
        else:
            raise ValueError(f'Failed to identify the subshells of the valence of {pseudo_file}')
        all_orbitals = list(z_core_to_first_orbital.values()) + ['5d', '6p', '6d']
        expected_orbitals[atom.symbol] = sorted(all_orbitals[all_orbitals.index(first_orbital):])
    return expected_orbitals
