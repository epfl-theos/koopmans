import shutil

import pytest

from koopmans import base_directory
from koopmans.io import read, write_json
from koopmans.utils import chdir

tutorial_dir = base_directory / 'tutorials' / 'tutorial_3'


@pytest.mark.tutorials
def test_zno_pdos(tutorial_patch):
    with chdir(tutorial_dir):
        wf = read('zno.json', override={'workflow': {'task': 'dft_bands'}})
        wf.run()


@pytest.mark.tutorials
def test_replot_dft_bandstructure():
    with chdir(tutorial_dir):
        exec(open('replot_dft_bandstructure.py', 'r').read())


@pytest.mark.tutorials
def test_zno_wannierize(tutorial_patch):
    with chdir(tutorial_dir):
        wf = read('zno.json', override={'workflow': {'task': 'wannierize'}})
        wf.run()


@pytest.mark.tutorials
def test_zno_ki(tutorial_patch):
    with chdir(tutorial_dir / 'dfpt'):
        shutil.copy('../zno.json', '.')
        wf = read('zno.json')
        wf.run()


@pytest.mark.tutorials
def test_zno_plot_bandstructures():
    with chdir(tutorial_dir / 'dfpt'):
        exec(open('../plot_bandstructures.py', 'r').read())
