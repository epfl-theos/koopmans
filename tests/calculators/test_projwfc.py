import shutil
from pathlib import Path

import pytest

from koopmans import __path__ as koopmans_src
from koopmans import utils, workflows
from koopmans.io import read_kwf as read_encoded_json
from koopmans.io import write_kwf as write_encoded_json


from tests.helpers.patches import benchmark_filename


def test_generate_dos(silicon, tmp_path, datadir, pytestconfig):
    with utils.chdir(tmp_path):
        # Create a projwfc calculator to match the one that was used to generate the pdos files
        wf = workflows.DFTBandsWorkflow(
            parameters={'pseudo_library': 'pseudo_dojo_standard', 'base_functional': 'pbesol', 'from_scratch': True},
            name='si', **silicon)

        calc = wf.new_calculator('projwfc')
        calc.pseudo_dir = Path(koopmans_src[0]) / 'pseudopotentials/pseudo_dojo_standard_v0.4.1/pbesol'

        # Copy over pdos files
        for f in (datadir / 'projwfc').glob('*.pdos*'):
            shutil.copy(f, f.name)

        # Attempt to read pdos files
        calc.generate_dos()
        dos = calc.results['dos']

        if pytestconfig.getoption('generate_benchmark'):
            # Write the DOS to file
            with open(benchmark_filename(calc), 'w') as fd:
                write_encoded_json(dos, fd)
        else:
            # Compare with the DOS on file
            with open(benchmark_filename(calc), 'r') as fd:
                dos_ref = read_encoded_json(fd)
            assert dos == dos_ref
