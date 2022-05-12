import pytest
from pathlib import Path
from koopmans.io import read
from koopmans.utils import chdir


def test_run():
    with chdir(Path(__file__).parent):
        wf = read('h2o_conv.json')
        wf.run()


def test_plot():
    with chdir(Path(__file__).parent):
        exec(open('plot.py', 'r').read())
