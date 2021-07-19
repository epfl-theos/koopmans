"""
Defines a smart "command" class in place of storing these simply as strings

Written by Edward Linscott, Feb 2021
"""

import os


class Command(object):
    """
    A more flexible class for storing commands
    It has the attributes...
        Command.path
        Command.executable
        Command.flags
        Command.suffix
    """

    def __init__(self, value, **kwargs):
        if isinstance(value, Command):
            value = str(value)
        self.__set__(value)

        # Accepts setting of public attributes as kwargs
        for k, v in kwargs.items():
            assert hasattr(self, k), f'Unrecognised argument {k} provided to {self.__class__.__name__}'
            assert not k.startswith(
                '_'), f'Do not attempt to set private variables via {self.__class__.__name__}.__init__()'
            setattr(self, k, v)

    def __set__(self, value):
        self.flags = ''
        if isinstance(value, str):
            if value.startswith(('srun', 'mpirun')):
                raise ValueError(
                    'You tried to set the command for the serial calculator {self.__class__.__name__} with an MPI call')
            [path_plus_executable, self.suffix] = value.split(' ', 1)
            if '/' in path_plus_executable:
                self.path, self.executable = path_plus_executable.rsplit('/', 1)
                self.path += '/'
            else:
                self.executable = path_plus_executable
                self.path = ''
        else:
            raise NotImplementedError(f'{self.__class__.__name__} must be set via a string')

    def __repr__(self):
        return ' '.join([self.path + self.executable, self.flags, self.suffix]).replace('  ', ' ')


class ParallelCommand(Command):
    """
    An extension to the Command class for mpi-parallelized executables
    """

    def __get__(self):
        return self.mpi_command + ' ' + super().__get__()

    def __set__(self, value):
        if isinstance(value, str):
            default_mpi_command = os.environ.get('PARA_PREFIX', None)
            if value.startswith('srun'):
                [self.mpi_command, rest_of_command] = value.split(' ', 1)
            elif value.startswith('mpirun -np '):
                [mpirun, np, np_num, rest_of_command] = value.split(' ', 3)
                self.mpi_command = f'mpirun -np {np_num}'
            elif default_mpi_command is not None:
                self.mpi_command = default_mpi_command
                rest_of_command = value
            else:
                self.mpi_command = ''
                rest_of_command = value
            super().__set__(rest_of_command)
        else:
            raise NotImplementedError(f'{self.__class__.__name__} must be set via a string')

    def __repr__(self):
        return self.mpi_command + ' ' + super().__repr__()


class ParallelCommandWithPostfix(ParallelCommand):
    """
    An extension to the Parallel Command class that supports Quantum ESPRESSO's $PARA_POSTFIX
    """

    @property
    def flags(self):
        return self.postfix + ' ' + self._flags

    @flags.setter
    def flags(self, value):
        self._flags = value

    def __init__(self, *args, **kwargs):
        self.postfix = os.environ.get('PARA_POSTFIX', '')
        super().__init__(*args, **kwargs)
