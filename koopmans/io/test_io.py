'''

Tests for koopmans.io

'''

import os
from ase.build import molecule
from koopmans.workflows import SinglepointWorkflow
from koopmans.io import write, read


def test_write_then_read_json():

    atoms = molecule('H2O', vacuum=10, tags=['0', '1', '1'])

    workflow_out = SinglepointWorkflow(atoms, parameters={'pseudo_library': 'sg15'}, name='test')

    write(workflow_out, 'test.json')

    workflow_in = read('test.json')

    assert workflow_out == workflow_in

    os.remove('test.json')

    return
