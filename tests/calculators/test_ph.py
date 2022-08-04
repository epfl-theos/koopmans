import shutil

import numpy as np
import pytest

from koopmans import base_directory, utils, workflows
from koopmans.io import read_kwf as read_encoded_json
from koopmans.io import write_kwf as write_encoded_json
from koopmans.testing import benchmark_filename


def test_read_dynG(tio2, tmp_path, datadir, pytestconfig):
    with utils.chdir(tmp_path):
        # Create a ph calculator to match the one that was used to generate the dynG file
        wf = workflows.DFTPhWorkflow(
            parameters={'pseudo_library': 'pseudo_dojo_standard', 'base_functional': 'pbesol', 'from_scratch': True},
            name='tio2', **tio2)

        calc = wf.new_calculator('ph', epsil=True, fildyn=f'{wf.name}.dynG')

        # Copy over ph files
        for f in (datadir / 'ph').glob('*.dynG'):
            shutil.copy(f, f.name)

        # Attempt to read dynG file
        calc.read_dynG()
        eps = calc.results['dielectric tensor']

        if pytestconfig.getoption('generate_benchmark'):
            # Write the dielectric tensor to file
            with open(benchmark_filename(calc), 'w') as fd:
                write_encoded_json(eps, fd)
        else:
            # Compare with the dielectric tensor on file
            with open(benchmark_filename(calc), 'r') as fd:
                eps_ref = read_encoded_json(fd)
            assert np.allclose(eps, eps_ref)
