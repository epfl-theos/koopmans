import subprocess
from abc import abstractmethod
from typing import Generic

from ._process import InputModel, OutputModel, Process


class CommandLineTool(Process[InputModel, OutputModel], Generic[InputModel, OutputModel]):
    """A Process that involves running a command on the command line."""

    @property
    @abstractmethod
    def command(self) -> str:
        """Return the command the executes the process."""
        ...

    @abstractmethod
    def _set_outputs(self):
        """Set the outputs of the process."""
        ...

    def _run(self):
        with self.engine.chdir(self.directory):
            # Run the command within self.directory
            ierr = subprocess.call(self.command, shell=True)
            if ierr > 0:
                raise OSError(f'`{self.command}` exited with exit code {ierr}')
        self._set_outputs()
