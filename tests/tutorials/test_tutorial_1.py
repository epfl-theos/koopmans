import pytest

from koopmans import base_directory
from koopmans.io import read
from koopmans.utils import chdir

tutorial_dir = base_directory / 'tutorials' / 'tutorial_1'


@pytest.mark.tutorials
def test_run(tutorial_patch):
    with chdir(tutorial_dir):
        exec(open('run.py', 'r').read())


@pytest.mark.tutorials
def test_read():
    with chdir(tutorial_dir):
        exec(open('read.py', 'r').read())
