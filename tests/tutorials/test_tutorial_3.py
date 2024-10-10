import shutil
from pathlib import Path

import pytest

from koopmans.io import read
from koopmans.utils import chdir

tutorial_dir = Path(__file__).parents[2] / 'tutorials' / 'tutorial_3'


@pytest.mark.filterwarnings("ignore::koopmans.utils.CalculatorNotConvergedWarning")
@pytest.mark.tutorials
def test_zno_ki(tutorial_patch):
    with chdir(tutorial_dir / '01-ki'):
        wf = read('zno.json')
        wf.run()


@pytest.mark.tutorials
def test_zno_plot_bandstructures():
    with chdir(tutorial_dir / '01-ki'):
        exec(open('plot_bandstructures.py', 'r').read())


@pytest.mark.tutorials
def test_zno_pdos(tutorial_patch):
    with chdir(tutorial_dir / '02-dft_bands'):
        wf = read('zno.json')
        wf.run()


@pytest.mark.tutorials
def test_replot_dft_bandstructure():
    with chdir(tutorial_dir / '02-dft_bands'):
        exec(open('replot_dft_bandstructure.py', 'r').read())


@pytest.mark.filterwarnings("ignore::koopmans.utils.CalculatorNotConvergedWarning")
@pytest.mark.tutorials
def test_zno_wannierize(tutorial_patch):
    with chdir(tutorial_dir / '03-wannierize'):
        wf = read('zno.json')
        wf.run()
