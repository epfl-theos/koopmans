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
            "keep_tmpdirs": False,
            "pseudo_library": "pseudo_dojo_standard"}
        wf = workflows.WannierizeWorkflow(parameters=parameters, **tio2)
        wf.run()
