import subprocess
from abc import abstractmethod

from koopmans.commands import Command
from koopmans.utils import system_call

from ._process import Process


class CommandLineTool(Process):

    @property
    @abstractmethod
    def command(self) -> Command:
        ...

    @abstractmethod
    def _pre_run(self):
        ...

    @abstractmethod
    def _post_run(self):
        ...

    def _run(self):
        self._pre_run()
        ierr = subprocess.call(str(self.command), shell=True)
        if ierr > 0:
            raise OSError(f'`{self.command}` exited with exit code {ierr}')
        self._post_run()
