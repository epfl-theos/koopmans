'''
Script for running test with pytest

Written by Edward Linscott, Jan 2021
'''

import pytest
from tests.conftest import WorkflowTest


json = 'tests/test_03/test_si_ki.json'


@pytest.mark.mock
def test_mock_si_ki(capsys, mock_quantum_espresso):
    test = WorkflowTest(json, capsys, mock=True)
    test.run()


def test_si_ki(capsys):
    test = WorkflowTest(json, capsys)
    test.run()
