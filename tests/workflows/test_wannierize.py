from pathlib import Path

import numpy as np
import pytest

from koopmans import workflows
from koopmans.kpoints import Kpoints
from koopmans.utils import chdir


def test_wannierize_tio2(tio2, tmp_path, sys2file, workflow_patch):
    with chdir(tmp_path):
        parameters = {
            "init_orbitals": "mlwfs",
            "init_empty_orbitals": "projwfs",
        }
        tio2["pseudo_library"] = "PseudoDojo/0.4/PBE/SR/standard/upf"
        wf = workflows.WannierizeWorkflow(parameters=parameters, **tio2)
        wf.run()
