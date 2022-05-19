import pytest
from ase.spacegroup import crystal
from koopmans.utils import chdir
from koopmans.io import write
from koopmans import workflows
from koopmans.projections import ProjectionBlocks


a = 4.6
c = 2.95

# Rutile TiO2:
tio2 = crystal(['Ti', 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
               spacegroup=136, cellpar=[a, a, c, 90, 90, 90])

tio2_projs = ProjectionBlocks.fromprojections([["Ti:l=0"], ["Ti:l=1"], ["O:l=0"], ["O:l=1"], ["Ti:l=0", "Ti:l=2"]],
                                              fillings=[True, True, True, True, False],
                                              spins=[None, None, None, None, None],
                                              atoms=tio2)

tio2_master_calc_params = {'pw': {'ecutwfc': 40, 'nbnd': 36}}


def test_wannierize_tio2(tmp_path, sys2file):
    with chdir(tmp_path):
        parameters = {
            "init_orbitals": "mlwfs",
            "init_empty_orbitals": "projwfs",
            "pseudo_library": "sg15_v1.0"}
        wf = workflows.WannierizeWorkflow(tio2, parameters, tio2_master_calc_params, kgrid=[
            2, 2, 2], kpath='GXG', projections=tio2_projs)
        wf.run()
