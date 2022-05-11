import pytest
from koopmans.io import read
from koopmans.utils import chdir
from pathlib import Path


def test_wannierise():
    with chdir(Path(__file__).parent):
        exec(open('wannierise.py', 'r').read())


def test_plot_pickled_fig():
    with chdir(Path(__file__).parent):
        exec(open('load_pickled_figure.py', 'r').read())


def test_run_ki():
    with chdir(Path(__file__).parent):
        wf = read('si.json', override={'workflow': {'task': 'singlepoint'}})
        wf.parameters.task = 'singlepoint'
        wf.run()


def test_plot_bs():
    with chdir(Path(__file__).parent):
        exec(open('plot_bandstructure.py', 'r').read())
