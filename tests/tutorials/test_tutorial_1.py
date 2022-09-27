import pytest

from pathlib import Path
from koopmans.utils import chdir

tutorial_dir = Path(__path__[0]).parents[2] / 'tutorials' / 'tutorial_1'


@pytest.mark.tutorials
def test_run(tutorial_patch):
    with chdir(tutorial_dir):
        exec(open('run.py', 'r').read())


@pytest.mark.tutorials
def test_read():
    with chdir(tutorial_dir):
        exec(open('read.py', 'r').read())
