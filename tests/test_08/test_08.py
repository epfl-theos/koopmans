'''
Script for running test with pytest

Written by Edward Linscott, Jan 2021
'''

import pytest
from conftest import WorkflowTest


json = 'tests/test_08/test_si_ui.json'


@pytest.mark.mock
def test_mock_si_ui(capsys):
    test = WorkflowTest(json, capsys)
    test.run()


@pytest.mark.stumble
def test_si_ui_stumble(capsys, stumble):
    test = WorkflowTest(json, capsys)
    test.run()


@pytest.mark.standard
def test_si_ui(capsys):
    test = WorkflowTest(json, capsys)
    test.run()
