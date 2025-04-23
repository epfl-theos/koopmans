import subprocess
from abc import abstractmethod

from koopmans.commands import Command

from ._process import Process


class CommandLineTool(Process):
    """A Process that involves running a command on the command line."""

    @property
    @abstractmethod
    def command(self) -> Command:
        """Return the command the executes the process."""
        ...

    @abstractmethod
    def _set_outputs(self):
        """Set the outputs of the process."""
        ...

    def _run(self):
        with self.engine.chdir(self.directory):
            # Run the command within self.directory
            ierr = subprocess.call(str(self.command), shell=True)
            if ierr > 0:
                raise OSError(f'`{self.command}` exited with exit code {ierr}')
        self._set_outputs()
