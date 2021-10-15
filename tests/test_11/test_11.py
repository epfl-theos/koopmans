'''
Script for running test with pytest

Written by Edward Linscott, Jan 2021
'''

import pytest
from conftest import WorkflowTest


json = 'tests/test_11/test_si_wannierise.json'


@pytest.mark.mock
def test_mock_si_wannierise(capsys, mock_quantum_espresso):
    test = WorkflowTest(json, capsys, mock=True)
    test.run()


def test_si_wannierise(capsys):
    test = WorkflowTest(json, capsys)
    test.run()
