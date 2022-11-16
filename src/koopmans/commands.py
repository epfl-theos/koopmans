"""
Defines a smart "command" class in place of storing these simply as strings

Written by Edward Linscott, Feb 2021
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Union


class Command(object):
    """
    A more flexible class for storing commands
    It has the attributes...
        Command.path
        Command.executable
        Command.flags
        Command.suffix
    """

    def __init__(self, value: Optional[Union[str, Command]] = None, **kwargs) -> None:
        self._path = Path()
        self.executable: str = ''
        self._flags: str = ''
        self.suffix: str = ''

        if isinstance(value, Command):
            value = str(value)
        if value is not None:
            self.__set__(value)

        # Accepts setting of public attributes as kwargs
        for k, v in kwargs.items():
            assert hasattr(self, k), f'Unrecognized argument {k} provided to {self.__class__.__name__}'
            assert not k.startswith(
                '_'), f'Do not attempt to set private variables via {self.__class__.__name__}.__init__()'
            setattr(self, k, v)

    def __set__(self, value: str):
        self.flags = ''
        if isinstance(value, str):
            if value.startswith(('srun', 'mpirun')):
                raise ValueError(
                    'You tried to set the command for the serial calculator {self.__class__.__name__} with an MPI call')
            [path_plus_executable, self.suffix] = value.split(' ', 1)
            if '/' in path_plus_executable:
                path, self.executable = path_plus_executable.rsplit('/', 1)
                self.path = Path(path)
            else:
                self.executable = path_plus_executable
                self.path = Path()
        else:
            raise NotImplementedError(f'{self.__class__.__name__} must be set via a string')

    def __repr__(self):
        return ' '.join([str(self.path / self.executable), self.flags, self.suffix]).replace('  ', ' ')

    @property
    def path(self) -> Path:
        return self._path

    @path.setter
    def path(self, value: Union[Path, str]):
        self._path = Path(value)

    @property
    def flags(self) -> str:
        return self._flags

    @flags.setter
    def flags(self, value: str):
        self._flags = value

    def todict(self):
        dct = self.__dict__.copy()
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @classmethod
    def fromdict(cls, dct):
        command = cls(**{k.lstrip('_'): v for k, v in dct.items()})
        return command

    def __eq__(self, other):
        return self.__dict__ == other.__dict__


class ParallelCommand(Command):
    """
    An extension to the Command class for mpi-parallelized executables
    """

    def __init__(self, *args, **kwargs) -> None:
        self.mpi_command: str = ''
        super().__init__(*args, **kwargs)

    def __get__(self):
        return self.mpi_command + ' ' + super().__get__()

    def __set__(self, value: str):
        if isinstance(value, str):
            default_mpi_command = os.environ.get('PARA_PREFIX', None)
            if value.startswith('srun') or value.startswith('mpirun'):
                splitval = value.split()

                i_command = None
                for i, val in enumerate(splitval[1:]):
                    prev_val = splitval[i]
                    if val.startswith('-'):
                        continue
                    elif prev_val.startswith('-') and '=' not in prev_val:
                        continue
                    else:
                        i_command = i + 1
                        break
                if i_command is None:
                    raise ValueError(f'Failed to parse {value}')

                self.mpi_command = ' '.join(splitval[:i_command])
                rest_of_command = ' '.join(splitval[i_command:])
            elif default_mpi_command is not None:
                self.mpi_command = default_mpi_command
                rest_of_command = value
            else:
                self.mpi_command = ''
                rest_of_command = value
            super().__set__(rest_of_command)
        else:
            raise NotImplementedError(f'{self.__class__.__name__} must be set via a string')

    def __repr__(self) -> str:
        return self.mpi_command + ' ' + super().__repr__()


class ParallelCommandWithPostfix(ParallelCommand):
    """
    An extension to the Parallel Command class that supports Quantum ESPRESSO's $PARA_POSTFIX
    """

    @property
    def flags(self) -> str:
        return (self.postfix + ' ' + self._flags).strip()

    @flags.setter
    def flags(self, value: str):
        self._flags = value

    def __init__(self, *args, **kwargs):
        self.postfix = os.environ.get('PARA_POSTFIX', '')
        super().__init__(*args, **kwargs)
