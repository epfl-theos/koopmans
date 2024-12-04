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

from ase_koopmans import Atoms
from upf_tools import UPFDict


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
    pseudo_file = pseudo_file.relative_to(pseudos_directory)
    library = str(pseudo_file.parent)
    version, functional, protocol = pseudo_file.parts[1:4]
    splitname = re.split(r'\.|_|-', pseudo_file.name)[0]
    element = splitname[0].upper() + splitname[1:].lower()
    citations: List[str] = []

    if library.startswith('SG15') or library.startswith('PseudoDojo'):
        kind = 'norm-conserving'
    else:
        kind = 'unknown'

    if library.startswith('SG15'):
        citations.append('Hamann2013')
        citations.append('Schlipf2015')
        if 'FR' in library:
            citations.append('Scherpelz2016')
    elif library.startswith('PseudoDojo'):
        citations.append('Hamann2013')
        citations.append('vanSetten2018')

    pseudo_database.append(Pseudopotential(pseudo_file.name, element, pseudo_file.parent,
                                           functional.lower(), library, kind, citations))


def pseudos_library_directory(pseudo_library: str) -> Path:
    return pseudos_directory / pseudo_library


def fetch_pseudo(**kwargs: Any) -> Pseudopotential:
    matches = [psp for psp in pseudo_database if all([getattr(psp, k) == v for k, v in kwargs.items()])]
    request_str = ', '.join([f'{k} = {v}' for k, v in kwargs.items()])
    if len(matches) == 0:
        raise ValueError(f'Could not find a pseudopotential in the database matching `{request_str}`')
    elif len(matches) > 1:
        raise ValueError(f'Found multiple pseudopotentials in the database matching `{request_str}`')
    else:
        return matches[0]


def read_pseudo_file(filename: Path) -> UPFDict:
    '''

    Reads in settings from a .upf file

    '''

    if not filename.exists():
        raise FileNotFoundError(f'Could not find the pseudopotential file `{filename}`')

    upf = UPFDict.from_upf(filename)

    return upf


def nelec_from_pseudos(atoms: Atoms, pseudopotentials: Dict[str, UPFDict]) -> int:
    '''
    Determines the number of electrons in the system using information from pseudopotential files
    '''

    valences_dct = {key: int(value['header']['z_valence']) for key, value in pseudopotentials.items()}

    if len(set(atoms.get_tags())) > 1:
        labels = [s + str(t) if t > 0 else s for s, t in zip(atoms.symbols, atoms.get_tags())]
    else:
        labels = atoms.symbols

    valences = [valences_dct[l] for l in labels]
    return sum(valences)


def expected_subshells(atoms: Atoms, pseudopotentials: Dict[str, UPFDict]) -> Dict[str, List[str]]:
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
        label = atom.symbol + str(atom.tag) if atom.tag > 0 else atom.symbol
        if label in expected_orbitals:
            continue
        z_core = atom.number - int(pseudopotentials[label]['header']['z_valence'])
        if z_core in z_core_to_first_orbital:
            first_orbital = z_core_to_first_orbital[z_core]
        else:
            raise ValueError(f'Failed to identify the subshells of the valence of `{pseudopotentials[label].filename}`')
        all_orbitals = list(z_core_to_first_orbital.values()) + ['5d', '6p', '6d']
        expected_orbitals[label] = sorted(all_orbitals[all_orbitals.index(first_orbital):])
    return expected_orbitals
