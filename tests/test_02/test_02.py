'''
Script for running test with pytest

Written by Edward Linscott, Jan 2021
'''

import pytest
from conftest import WorkflowTest


json = 'tests/test_02/test_co2_ki_with_empty.json'


@pytest.mark.mock
def test_mock_co2_ki_with_empty(capsys, mock_quantum_espresso):
    test = WorkflowTest(json, capsys, mock=True)
    test.run()


@pytest.mark.stumble
def test_co2_ki_with_empty_stumble(capsys, stumble):
    test = WorkflowTest(json, capsys)
    test.run()


@pytest.mark.standard
def test_co2_ki_with_empty(capsys):
    test = WorkflowTest(json, capsys)
    test.run()
