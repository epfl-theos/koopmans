"""Test tutorial 2."""

from pathlib import Path

import pytest

from koopmans.io import read
from koopmans.utils import chdir

tutorial_dir = Path(__file__).parents[2] / 'tutorials' / 'tutorial_2'


@pytest.mark.tutorials
def test_wannierize(tutorial_patch):
    """Test the Wannierization calculation."""
    with chdir(tutorial_dir):
        exec(open('wannierize.py', 'r').read())


@pytest.mark.tutorials
def test_plot_pickled_fig():
    """Test the plotting of the pickled figure."""
    with chdir(tutorial_dir):
        exec(open('load_pickled_figure.py', 'r').read())


@pytest.mark.tutorials
def test_run_ki(tutorial_patch):
    """Test the KI calculation."""
    with chdir(tutorial_dir):
        wf = read('si.json', override={'workflow': {'task': 'singlepoint'}})
        wf.parameters.task = 'singlepoint'
        wf.run()


@pytest.mark.tutorials
def test_plot_bs():
    """Test the script for plotting the band structure."""
    with chdir(tutorial_dir):
        exec(open('plot_bandstructure.py', 'r').read())
