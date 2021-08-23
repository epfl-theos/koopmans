'''
Script for running test with pytest

Written by Edward Linscott, Jan 2021
'''

import pytest
from conftest import WorkflowTest
import matplotlib
matplotlib.use('Agg')


json = 'tests/test_04/test_h2o_pkipz.json'


@pytest.mark.mock
def test_mock_h2o_pkipz(capsys, mock_quantum_espresso):
    test = WorkflowTest(json, capsys, mock=True)
    test.run()


def test_h2o_pkipz(capsys):
    test = WorkflowTest(json, capsys)
    test.run()
