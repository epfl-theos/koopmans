'''

Tests for koopmans.io._json

'''

from ase.build import molecule
from koopmans import utils
from koopmans.io import read, write
from koopmans.workflows import (ConvergenceMLWorkflow, SinglepointWorkflow,
                                TrajectoryWorkflow)


def test_write_then_read_json(tmpdir):
    with utils.chdir(tmpdir):
        atoms = molecule('H2O', vacuum=10, tags=['0', '1', '1'])
        workflow_out = SinglepointWorkflow(atoms, parameters={'pseudo_library': 'sg15'}, name='test')
        write(workflow_out, 'test.json')
        workflow_in = read('test.json')
        assert workflow_out == workflow_in


def test_write_then_read_trajectory_json(tmpdir):
    with utils.chdir(tmpdir):
        atoms = molecule('H2O', vacuum=10, tags=['0', '1', '1'])
        snapshots = 2*[atoms]
        workflow_out = TrajectoryWorkflow(
            snapshots, parameters={'task': 'trajectory', 'pseudo_library': 'sg15'}, name='test')

        write(workflow_out, 'test.json')
        workflow_in = read('test.json')

        assert workflow_out == workflow_in


def test_write_then_read_convergence_ml_json(tmpdir):
    with utils.chdir(tmpdir):
        atoms = molecule('H2O', vacuum=10, tags=['0', '1', '1'])
        snapshots = 2*[atoms]
        workflow_out = ConvergenceMLWorkflow(
            snapshots, parameters={'task': 'convergence_ml', 'pseudo_library': 'sg15'}, name='test')

        write(workflow_out, 'test.json')
        workflow_in = read('test.json')

        assert workflow_out == workflow_in
