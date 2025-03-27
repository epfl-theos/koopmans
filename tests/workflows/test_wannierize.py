"""Tests for the Wannierize workflow."""

import pytest  # noqa

from koopmans import workflows
from koopmans.utils import chdir


def test_wannierize_tio2(tio2, tmp_path, sys2file, workflow_patch):
    """Test the Wannierize workflow on TiO2."""
    with chdir(tmp_path):
        parameters = {
            "init_orbitals": "mlwfs",
            "init_empty_orbitals": "projwfs",
        }
        tio2["pseudo_library"] = "PseudoDojo/0.4/PBE/SR/standard/upf"
        wf = workflows.WannierizeWorkflow(parameters=parameters, **tio2)
        wf.run()
