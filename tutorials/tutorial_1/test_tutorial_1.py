import pytest
from koopmans.io import read
from koopmans.utils import chdir
from pathlib import Path


def test_run():
    with chdir(Path(__file__).parent):
        exec(open('run.py', 'r').read())


def test_read():
    with chdir(Path(__file__).parent):
        exec(open('run.py', 'r').read())
