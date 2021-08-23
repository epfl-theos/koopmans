'''
Script for running test with pytest

Written by Edward Linscott, Jan 2021
'''

import pytest
from conftest import WorkflowTest
import matplotlib
matplotlib.use('Agg')


json = 'tests/test_07/test_si_all.json'


@pytest.mark.mock
def test_mock_si_all(capsys, mock_quantum_espresso):
    test = WorkflowTest(json, capsys, mock=True)
    test.run()


def test_si_all(capsys):
    test = WorkflowTest(json, capsys)
    test.run()
