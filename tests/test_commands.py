import os

import pytest

from koopmans import utils
from koopmans.commands import (Command, ParallelCommand,
                               ParallelCommandWithPostfix)


def test_command():
    # Creation from a string
    c_str = 'pw.x -in in.pwi > out.pwi'
    c = Command(c_str)
    assert c.executable == 'pw.x'
    assert str(c) == c_str

    # Creation from another command object
    c2 = Command(c)
    assert c.executable == 'pw.x'
    assert str(c) == c_str


parallel_commands = [('mpirun -n 16', 'pw.x', '-in in.pwi > out.pwi'),
                     ('srun', 'pw.x', '-in in.pwi > out.pwo 2>&1'),
                     ('srun --mpi=pmi2 -n 48', 'pw.x', '-in in.pwi > out.pwo 2>&1')]


@pytest.mark.parametrize('mpi_command,executable,suffix', parallel_commands)
def test_parallel_command(mpi_command, executable, suffix):
    c_str = ' '.join((mpi_command, executable, suffix))
    c = ParallelCommand(c_str)

    assert c.executable == executable
    assert c.mpi_command == mpi_command
    assert str(c) == c_str


def test_parallel_command_with_postfix():
    with utils.set_env(PARA_POSTFIX='-npool 4'):
        c_str = 'mpirun -n 16 pw.x -in in.pwi > out.pwi'
        postfix = '-npool 4'
        os.environ['PARA_POSTFIX'] = postfix
        c = ParallelCommandWithPostfix(c_str)

        assert c.postfix == postfix
        assert c.flags == postfix
