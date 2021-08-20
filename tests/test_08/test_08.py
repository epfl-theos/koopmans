'''
Script for running test with pytest

Written by Edward Linscott, Jan 2021
'''

import pytest
from conftest import WorkflowTest


json = 'tests/test_08/test_si_ui.json'


# Note that this is identical to the non-mock test because this calculation relies on python alone, and not QE
@pytest.mark.mock
def test_mock_si_ui(capsys):
    test = WorkflowTest(json, capsys)
    test.run()


def test_si_ui(capsys):
    test = WorkflowTest(json, capsys)
    test.run()
