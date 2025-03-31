"""Tests for `koopmans.calculators.projwfc`."""

from pathlib import Path

import pytest  # noqa: F401

from koopmans import utils, workflows
from koopmans.files import File, LocalFile
from koopmans.io import read_pkl, write_pkl
from tests.helpers.patches import benchmark_filename


def test_generate_dos(silicon, tmp_path, datadir, pytestconfig):
    """Test the generation of DOS files."""
    with utils.chdir(tmp_path):
        # Create a projwfc calculator to match the one that was used to generate the pdos files
        wf = workflows.DFTBandsWorkflow(
            parameters={'pseudo_library': 'PseudoDojo/0.4/PBEsol/SR/standard/upf'},
            name='si', **silicon)
        wf.directory = Path()
        calc = wf.new_calculator('projwfc')
        calc.directory = Path()

        # Make sure the pseudopotential files exist where the calculator will expect them to be
        pseudo_dir = LocalFile(calc.directory / calc.parameters.outdir / (calc.parameters.prefix + '.save'))
        wf.engine.mkdir(pseudo_dir, parents=True)
        for psp in wf.pseudopotentials.values():
            wf.engine.copy_file(LocalFile(psp.filename), pseudo_dir)

        # Copy over pdos files
        for f in (datadir / 'projwfc').glob('*.pdos*'):
            wf.engine.copy_file(LocalFile(f), File(wf, f.name))

        # Attempt to read pdos files
        calc.generate_dos()
        dos = calc.results['dos']

        if pytestconfig.getoption('generate_benchmark'):
            # Write the DOS to file
            write_pkl(dos, benchmark_filename(calc))
        else:
            # Compare with the DOS on file
            dos_ref = read_pkl(benchmark_filename(calc))
            assert dos._almost_equals(dos_ref)
