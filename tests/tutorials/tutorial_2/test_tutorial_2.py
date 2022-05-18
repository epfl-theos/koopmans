import pytest
from koopmans.io import read
from koopmans.utils import chdir
from koopmans import base_directory


tutorial_dir = base_directory / 'tutorials' / 'tutorial_2'


def test_wannierise():
    with chdir(tutorial_dir):
        exec(open('wannierise.py', 'r').read())


def test_plot_pickled_fig():
    with chdir(tutorial_dir):
        exec(open('load_pickled_figure.py', 'r').read())


def test_run_ki():
    with chdir(tutorial_dir):
        wf = read('si.json', override={'workflow': {'task': 'singlepoint'}})
        wf.parameters.task = 'singlepoint'
        wf.run()


def test_plot_bs():
    with chdir(tutorial_dir):
        exec(open('plot_bandstructure.py', 'r').read())
