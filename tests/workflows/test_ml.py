"""
"""

import pytest

from koopmans.io import write
from koopmans.projections import ProjectionBlocks
from koopmans.utils import chdir
from koopmans.workflows import TrajectoryWorkflow


@pytest.mark.parametrize('descriptor', ['orbital_density', 'self_hartree'])
@pytest.mark.parametrize('estimator', ['mean', 'ridge_regression'])
@pytest.mark.parametrize('occ_and_emp_together', [True, False])
def test_ml_train_water(tmp_path, workflow_patch, water_snapshots, descriptor, estimator, occ_and_emp_together):
    with chdir(tmp_path):
        projs = ProjectionBlocks.fromlist([[{'site': 'O', 'ang_mtm': 'sp3'}], [{'site': 'H', 'ang_mtm': 's'}]], spins=[
                                          None, None], atoms=water_snapshots['atoms'])
        water_snapshots.pop('nbnd')
        water_snapshots['ecutwfc'] = 50.0
        wf = TrajectoryWorkflow(**water_snapshots, descriptor=descriptor, estimator=estimator, occ_and_emp_together=occ_and_emp_together, train=True,
                                init_orbitals='mlwfs', init_empty_orbitals='mlwfs', projections=projs, dis_froz_max=0.0, calculator_parameters={'pw': {'nbnd': 12}}, calculate_bands=False)
        write(wf, 'input.json')
        wf.run()
