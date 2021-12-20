'''
Script for running test with pytest

Written by Edward Linscott, Jan 2021
'''

import pytest
from conftest import WorkflowTest


json = 'tests/test_11/test_tio2_wannierise.json'


@pytest.mark.mock
def test_mock_tio2_wannierise(capsys, mock_quantum_espresso):
    test = WorkflowTest(json, capsys, mock=True)
    test.run()


@pytest.mark.stumble
def test_tio2_wannierise_stumble(capsys, stumble):
    test = WorkflowTest(json, capsys)
    test.run()


@pytest.mark.standard
def test_tio2_wannierise(capsys):
    test = WorkflowTest(json, capsys)
    test.run()
