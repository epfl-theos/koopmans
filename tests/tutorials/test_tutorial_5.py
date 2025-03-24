from pathlib import Path

import numpy as np
import pytest

from koopmans.io import read
from koopmans.utils import chdir

tutorial_dir = Path(__file__).parents[2] / 'tutorials' / 'tutorial_5'
benchmark_dir = Path(__file__).parents[1] / 'benchmarks'


@pytest.mark.tutorials
def test_train(tutorial_patch, tmpdir, pytestconfig, sys2file):
    with chdir(tutorial_dir / '01-train'):
        wf = read('h2o_train.json')
        wf.run()


@pytest.mark.tutorials
def test_predict(tutorial_patch, tmpdir, pytestconfig, sys2file):
    with chdir(tutorial_dir / '02-predict'):
        wf = read('h2o_predict.json')
        wf.run()


@pytest.mark.tutorials
def test_predict_plot(tutorial_patch, tmpdir, pytestconfig, sys2file):
    with chdir(tutorial_dir / '02-predict'):
        exec(open('plot.py').read())


@pytest.mark.tutorials
def test_testing(tutorial_patch, tmpdir, pytestconfig):
    with chdir(tutorial_dir / '03-test'):
        wf = read('h2o_test.json')
        wf.run()


@pytest.mark.tutorials
def test_testing_plot(tutorial_patch, tmpdir, pytestconfig):
    with chdir(tutorial_dir / '03-test'):
        exec(open('plot.py').read())


@pytest.mark.tutorials
def test_advanced_testing(tutorial_patch, tmpdir, pytestconfig):
    with chdir(tutorial_dir / '04-advanced-testing'):
        exec(open('train_and_verify.py').read())
