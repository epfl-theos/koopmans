from pathlib import Path

import pytest

from koopmans import base_directory
from koopmans.io import read
from koopmans.utils import chdir

tutorial_dir = base_directory / 'tutorials' / 'tutorial_5'

# TODO: Yannick


@pytest.mark.tutorials
def test_run_trajectory(tutorial_patch):
    with chdir(tutorial_dir):
        wf = read('h2o_trajectory_ml.json')
        wf.run()
