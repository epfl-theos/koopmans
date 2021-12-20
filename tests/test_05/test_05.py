'''
Script for running test with pytest

Written by Edward Linscott, Jan 2021
'''

import pytest
from conftest import WorkflowTest


json = 'tests/test_05/test_h2o_all.json'


@pytest.mark.mock
def test_mock_h2o_all(capsys, mock_quantum_espresso):
    test = WorkflowTest(json, capsys, mock=True)
    test.run()


@pytest.mark.stumble
def test_h2o_all_stumble(capsys, stumble):
    test = WorkflowTest(json, capsys)
    test.run()


@pytest.mark.standard
def test_h2o_all(capsys):
    test = WorkflowTest(json, capsys)
    test.run()
