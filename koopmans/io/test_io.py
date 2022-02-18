'''

Tests for koopmans.io

'''

from ase.build import molecule
from koopmans.workflows import SinglepointWorkflow
from koopmans.io import write


def test_write_json():

    atoms = molecule('H2O', vacuum=10, tags=['0', '1', '1'])

    workflow = SinglepointWorkflow(atoms, parameters={'pseudo_library': 'sg15'})

    write(workflow, 'test.json')

    return
