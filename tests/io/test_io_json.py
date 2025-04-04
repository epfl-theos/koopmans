"""Tests for koopmans.io._json."""

from ase_koopmans.build import molecule

from koopmans import utils
from koopmans.engines.localhost import LocalhostEngine
from koopmans.io import read, write
from koopmans.workflows import SinglepointWorkflow, TrajectoryWorkflow


def test_write_then_read_json(tmpdir):
    """Test writing and reading a SinglepointWorkflow to and from a json file."""
    with utils.chdir(tmpdir):
        atoms = molecule('H2O', vacuum=10, tags=['0', '1', '1'])
        workflow_out = SinglepointWorkflow(atoms, engine=LocalhostEngine(), parameters={
                                           'pseudo_library': 'SG15/1.2/PBE/SR'}, name='test')
        write(workflow_out, 'test.json')
        workflow_in = read('test.json')
        assert workflow_out == workflow_in


def test_write_then_read_trajectory_json(tmpdir):
    """Test writing and reading a TrajectoryWorkflow to and from a json file."""
    with utils.chdir(tmpdir):
        atoms = molecule('H2O', vacuum=10, tags=['0', '1', '1'])
        snapshots = 2 * [atoms]
        workflow_out = TrajectoryWorkflow(
            atoms, engine=LocalhostEngine(), snapshots=snapshots,
            parameters={'task': 'trajectory', 'pseudo_library': 'SG15/1.2/PBE/SR'}, name='test'
        )

        write(workflow_out, 'test.json')
        workflow_in = read('test.json')

        # It is known that ASE's writing/reading from xyz files (which happens when we write/read a trajectory) is not
        # idempotent (because the atomic positions are written with fixed width), so for the purposes of this test we
        # will override the atoms and snapshots with the same object
        workflow_out.atoms = workflow_in.atoms
        workflow_out.snapshots = workflow_in.snapshots
        workflow_in.projections._atoms = workflow_out.projections._atoms

        assert workflow_out == workflow_in
