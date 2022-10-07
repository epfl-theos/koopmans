from pathlib import Path

import pytest

from koopmans.io import read
from koopmans.utils import chdir

tutorial_dir = Path(__file__).parents[2] / 'tutorials' / 'tutorial_2'


@pytest.mark.tutorials
def test_wannierize(tutorial_patch):
    with chdir(tutorial_dir):
        exec(open('wannierize.py', 'r').read())


@pytest.mark.tutorials
def test_plot_pickled_fig():
    with chdir(tutorial_dir):
        exec(open('load_pickled_figure.py', 'r').read())


@pytest.mark.tutorials
def test_run_ki(tutorial_patch):
    with chdir(tutorial_dir):
        wf = read('si.json', override={'workflow': {'task': 'singlepoint'}})
        wf.parameters.task = 'singlepoint'
        wf.run()


@pytest.mark.tutorials
def test_plot_bs():
    with chdir(tutorial_dir):
        exec(open('plot_bandstructure.py', 'r').read())
