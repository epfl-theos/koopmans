import shutil
from pathlib import Path

import pytest

from koopmans.io import read
from koopmans.utils import chdir

tutorial_dir = Path(__file__).parents[2] / 'tutorials' / 'tutorial_6'


@pytest.mark.tutorials
def test_cri3(tutorial_patch):
    with chdir(tutorial_dir):
        wf = read('cri3.json')
        wf.run()


@pytest.mark.tutorials
def test_cri3_plot_bandstructure():
    with chdir(tutorial_dir):
        exec(open('plot_bands.py', 'r').read())
