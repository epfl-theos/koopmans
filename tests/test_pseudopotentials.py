from ase import Atoms

from koopmans import pseudopotentials


def test_expected_subshells():
    '''
    Check that the function expected_subshells doesn't encounter KeyErrors for any element
    '''
    for pseudo in pseudopotentials.pseudo_database:
        atoms = Atoms(pseudo.element)
        psps = {pseudo.element: pseudo.name}
        subshells = pseudopotentials.expected_subshells(atoms, psps, pseudo.path)
        assert len(subshells) > 0
