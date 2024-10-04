import shutil
from pathlib import Path

import pytest

from koopmans import __path__ as koopmans_src
from koopmans import utils, workflows
from koopmans.io import read_pkl, write_pkl
from tests.helpers.patches import benchmark_filename


def test_generate_dos(silicon, tmp_path, datadir, pytestconfig):
    with utils.chdir(tmp_path):
        # Create a projwfc calculator to match the one that was used to generate the pdos files
        wf = workflows.DFTBandsWorkflow(
            parameters={'pseudo_library': 'pseudo_dojo_standard', 'base_functional': 'pbesol', 'from_scratch': True},
            name='si', **silicon)
        wf.directory = Path()
        calc = wf.new_calculator('projwfc')
        calc.directory = Path()

        # Make sure the pseudopotential files exist where the calculator will expect them to be
        pseudo_dir = calc.directory / calc.parameters.outdir / (calc.parameters.prefix + '.save')
        pseudo_dir.mkdir(parents=True)
        assert wf.parameters.pseudo_directory is not None
        for psp in wf.pseudopotentials.values():
            utils.copy(wf.parameters.pseudo_directory / psp, pseudo_dir)

        # Copy over pdos files
        for f in (datadir / 'projwfc').glob('*.pdos*'):
            utils.copy(f, f.name)

        # Attempt to read pdos files
        calc.generate_dos()
        dos = calc.results['dos']

        if pytestconfig.getoption('generate_benchmark'):
            # Write the DOS to file
            write_pkl(dos, benchmark_filename(calc))
        else:
            # Compare with the DOS on file
            dos_ref = read_pkl(benchmark_filename(calc))
            assert dos == dos_ref
