"""Test tutorial 4."""

from pathlib import Path

import pytest

from koopmans.io import read
from koopmans.utils import chdir

tutorial_dir = Path(__file__).parents[2] / 'tutorials' / 'tutorial_4'


@pytest.mark.tutorials
def test_run(tutorial_patch):
    """Test the calculation."""
    with chdir(tutorial_dir):
        wf = read('h2o_conv.json')
        wf.run()


@pytest.mark.tutorials
def test_plot():
    """Test the plottings script."""
    with chdir(tutorial_dir):
        exec(open('plot.py', 'r').read(), globals())


@pytest.mark.tutorials
def test_custom_convergence(tutorial_patch):
    """Test the custom convergence script."""
    with chdir(tutorial_dir):
        exec(open('custom_convergence.py', 'r').read(), globals())
