from itertools import chain

import pytest
from ase_koopmans import Atoms

from koopmans import pseudopotentials

all_local_pseudos = [x for x in chain(pseudopotentials.local_base_directory.rglob(
    '*.upf'), pseudopotentials.local_base_directory.rglob('*.UPF'))]


@pytest.mark.parametrize('pseudo', all_local_pseudos)
def test_reading_pseudos(pseudo):
    '''
    Check that the function read_pseudo_file doesn't raise any exceptions
    '''
    pseudopotentials.read_pseudo_file(pseudo)


@pytest.mark.parametrize('pseudo', all_local_pseudos)
def test_element_from_pseudo_filename(pseudo):
    '''
    Check that the function element_from_pseudo_filename doesn't raise any exceptions
    '''
    element = pseudopotentials.element_from_pseudo_filename(pseudo.name)
    assert len(element) in [1, 2]


@pytest.mark.parametrize('library', pseudopotentials.local_libraries)
def test_pseudopotential_library_citations(library):
    '''
    Check that the function pseudopotential_library_citations doesn't raise any exceptions
    '''
    citations = pseudopotentials.pseudopotential_library_citations(library)
    assert len(citations) > 0


@pytest.mark.parametrize('pseudo', all_local_pseudos)
def test_expected_subshells(pseudo):
    '''
    Check that the function expected_subshells doesn't encounter KeyErrors for any element
    '''

    element = pseudopotentials.element_from_pseudo_filename(pseudo.name)
    pseudo_dict = pseudopotentials.read_pseudo_file(pseudo)
    subshells = pseudopotentials.expected_subshells(Atoms(element), {element: pseudo_dict})
    assert len(subshells) > 0
