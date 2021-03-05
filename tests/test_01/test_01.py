'''
Script for running test with pytest

Written by Edward Linscott, Jan 2021
'''

import pytest
from tests.conftest import WorkflowTest


json = 'tests/test_01/test_h2o_kipz.json'


@pytest.mark.mock
def test_mock_h2o_kipz(capsys, mock_quantum_espresso):
    test = WorkflowTest(json, capsys, mock=True)
    test.run()


def test_h2o_kipz(capsys):
    test = WorkflowTest(json, capsys)
    test.run()

