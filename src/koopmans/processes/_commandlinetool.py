import subprocess
from abc import abstractmethod

from koopmans import utils
from koopmans.commands import Command

from ._process import Process


class CommandLineTool(Process):

    @property
    @abstractmethod
    def command(self) -> Command:
        ...

    @abstractmethod
    def _set_outputs(self):
        # This method should be used to set self.outputs
        ...

    def _run(self):
        with self.engine.chdir(self.directory):
            # Run the command within self.directory
            ierr = subprocess.call(str(self.command), shell=True)
            if ierr > 0:
                raise OSError(f'`{self.command}` exited with exit code {ierr}')
        self._set_outputs()
