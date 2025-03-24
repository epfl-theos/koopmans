"""

Module containing functions relating to pseudopotentials

Written by Edward Linscott Jan 2020
Moved into _utils Aug 2021
Split into a separate module Sep 2021

"""

import re
from itertools import chain
from pathlib import Path
from typing import Dict, List, Optional, OrderedDict

from ase_koopmans import Atoms
from upf_tools import UPFDict
from upf_tools.projectors import Projectors

local_base_directory = (Path(__file__).parent / 'pseudopotentials').resolve()
local_libraries = set([str(f.parent.relative_to(local_base_directory))
                       for f in chain(local_base_directory.rglob('*.upf'),
                                      local_base_directory.rglob('*.UPF'))])


def pseudopotential_library_citations(library: str) -> List[str]:
    citations = []
    if library.startswith('SG15'):
        citations.append('Hamann2013')
        citations.append('Schlipf2015')
        if 'FR' in library:
            citations.append('Scherpelz2016')
    elif library.startswith('PseudoDojo'):
        citations.append('Hamann2013')
        citations.append('vanSetten2018')
    return citations


def element_from_pseudo_filename(filename: str) -> str:
    splitname = re.split(r'\.|_|-', filename)[0]
    element = splitname[0].upper() + splitname[1:].lower()
    return element


def read_pseudo_file(filename: Path) -> UPFDict:
    '''

    Reads in settings from a .upf file

    '''

    if not filename.exists():
        raise FileNotFoundError(f'Could not find the pseudopotential file `{filename}`')

    upf = UPFDict.from_upf(filename)

    return upf


def nelec_from_pseudos(atoms: Atoms, pseudopotentials: OrderedDict[str, UPFDict]) -> int:
    '''
    Determines the number of electrons in the system using information from pseudopotential files
    '''

    raise NotImplementedError('Remove this')
    # valences_dct = {key: int(value['header']['z_valence']) for key, value in pseudopotentials.items()}


def nwfcs_from_pseudos(atoms: Atoms, pseudopotentials: Dict[str, str],
                       pseudo_dir: Optional[Path] = None) -> int:
    '''
    Determines the number of wfcs in the system using information from pseudopotential files
    '''

    raise NotImplementedError('This function is not yet implemented')
    # for psp in pseudopotentials.values():
    #     if pseudo_contents(psp, pseudo_dir)['header']['number_of_wfc'] == 0:
    #         raise ValueError(f'The pseudopotential {psp} does not contain wavefunctions')

    # def extract_nwfc(dct: Dict[str, List[Dict[str, Any]]]) -> int:
    #     return sum([2 * int(wfc['angular_momentum']) + 1 for wfc in dct['atomic_wave_functions']])

    # return _quantity_from_pseudos(extract_nwfc, atoms, pseudopotentials, pseudo_dir)


def nwfcs_from_projectors(atoms: Atoms, pseudopotentials: Dict[str, str], projector_dir: Path) -> int:
    '''
    Determines the number of wfcs in the system using information from projector files
    '''

    # Construct a dict mapping the element name to the projector file
    proj_files = {key: projector_dir / (value.rsplit('.', 1)[0] + '.dat') for key, value in pseudopotentials.items()}

    # Construct a dict mapping the element name to the number of projectors
    nproj = {key: sum([2 * p.l + 1 for p in Projectors.from_file(value)]) for key, value in proj_files.items()}

    # Calculate the total number of projectors
    if len(set(atoms.get_tags())) > 1:
        labels = [s + str(t) if t > 0 else s for s, t in zip(atoms.symbols, atoms.get_tags())]
    else:
        labels = atoms.symbols

    return sum([nproj[l] for l in labels])


def cutoffs_from_pseudos(atoms: Atoms, pseudo_dir_in: Optional[Path] = None) -> Dict[str, float]:
    """
    Works out a recommended ecutwfc and ecutrho based on the contents of a cutoffs.json file
    located in the pseudopotential directory
    """

    raise NotImplementedError('This function is not yet implemented')
    # pseudo_dir = get_pseudo_dir(pseudo_dir_in)

    # cutoff_database_file = pseudo_dir / 'cutoffs.json'
    # if not cutoff_database_file.exists():
    #     return {}

    # with open(cutoff_database_file, 'r') as fd:
    #     cutoff_database = json.load(fd)

    # cutoffs = {}
    # for symbol in atoms.symbols:
    #     if symbol not in cutoff_database:
    #         return {}
    #     for k, v in cutoff_database[symbol].items():
    #         if k not in cutoffs:
    #             cutoffs[k] = v
    #         elif cutoffs[k] < v:
    #             cutoffs[k] = v

    # return cutoffs


def expected_subshells(atoms: Atoms, pseudopotentials: OrderedDict[str, UPFDict]) -> Dict[str, List[str]]:
    """
    Determine which subshells will make up the valences of a set of pseudopotentials.

    Returns
    -------
    OrderedDict[str, List[str]]
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
