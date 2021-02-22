'''
Script for running test with pytest

Written by Edward Linscott, Jan 2021
'''

import os
import pytest
from conftest import WorkflowTest
from koopmans.utils import find_executable
from koopmans.calculators.environ import environ_addon_is_installed


# Work out if environ is installed
pw_command = os.environ.get('ASE_ESPRESSO_COMMAND', 'pw.x')
if pw_command.startswith('mpirun'):
    isplit = 3
elif pw_command.startswith('srun'):
    isplit = 1
else:
    isplit = 0
pw_executable = find_executable(pw_command.split()[isplit])
if pw_executable is None:
    environ_installed = False
    reason = 'pw not installed'
else:
    qe_directory, _, _ = pw_executable.rsplit('/', 2)
    environ_installed = environ_addon_is_installed(qe_directory)
    reason = f'Environ addon is not installed for {pw_executable}'


json = 'tests/test_07/test_so2_pbe_delta_scf.json'


@pytest.mark.mock
def test_mock_so2_pbe_delta_scf(capsys, mock_quantum_espresso):
    test = WorkflowTest(json, capsys, mock=True)
    test.run()


@pytest.mark.skipif(not environ_installed, reason=reason)
def test_so2_pbe_delta_scf(capsys):
    test = WorkflowTest(json, capsys)
    test.run()

