from pathlib import Path

import pytest

from koopmans import base_directory
from koopmans.io import read
from koopmans.utils import chdir

tutorial_dir = base_directory / 'tutorials' / 'tutorial_3'


@pytest.mark.tutorials
def test_run(tutorial_patch):
    with chdir(tutorial_dir):
        wf = read('h2o_conv.json')
        wf.run()


@pytest.mark.tutorials
def test_plot():
    with chdir(tutorial_dir):
        exec(open('plot.py', 'r').read())
