"""Testing `koopmans.calculators.ph`."""

import shutil
from pathlib import Path

import numpy as np
import pytest  # noqa: F401

from koopmans import utils, workflows
from koopmans.io import read_pkl, write_pkl
from tests.helpers.patches import benchmark_filename


def test_read_dynG(tio2, tmp_path, datadir, pytestconfig):
    """Test reading dynG files."""
    with utils.chdir(tmp_path):
        # Create a ph calculator to match the one that was used to generate the dynG file
        wf = workflows.DFTPhWorkflow(
            parameters={'pseudo_library': 'PseudoDojo/0.4/PBEsol/SR/standard/upf'},
            name='tio2', **tio2)
        wf.directory = Path()

        calc = wf.new_calculator('ph', epsil=True, fildyn=f'{wf.name}.dynG')
        calc.directory = Path()

        # Copy over ph files
        for f in (datadir / 'ph').glob('*.dynG'):
            shutil.copy(f, f.name)

        # Attempt to read dynG file
        calc.read_dynG()
        eps = calc.results['dielectric tensor']

        if pytestconfig.getoption('generate_benchmark'):
            # Write the dielectric tensor to file
            write_pkl(eps, benchmark_filename(calc))
        else:
            # Compare with the dielectric tensor on file
            eps_ref = read_pkl(benchmark_filename(calc))
            assert np.allclose(eps, eps_ref)
