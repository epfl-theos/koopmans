'''
Script for running test with pytest

Written by Edward Linscott, Jan 2021
'''

import pytest
from conftest import WorkflowTest


json = 'tests/test_12/test_si_ks2odd.json'


@pytest.mark.mock
def test_mock_si_ks2odd(capsys, mock_quantum_espresso):
    test = WorkflowTest(json, capsys, mock=True)
    test.run()


@pytest.mark.stumble
def test_si_ks2odd_stumble(capsys, stumble):
    test = WorkflowTest(json, capsys)
    test.run()


@pytest.mark.standard
def test_si_ks2odd(capsys):
    test = WorkflowTest(json, capsys)
    test.run()
