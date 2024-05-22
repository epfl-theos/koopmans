"""Defines Process, a class for abstractly representing a Workflow/Calculator/other operation.

Inspired by CWL."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Tuple, Type

from pydantic import BaseModel

from koopmans import calculators, utils


class Process(ABC):

    def __init__(self, **kwargs):
        self.inputs = self._input_model(**kwargs)
        self.outputs = None
        self.directory: Path = Path()
        self._linked_files: Dict[str, Tuple[Process | calculators.Calc | None, Path]] = {}

    @property
    @abstractmethod
    def _input_model(self) -> Type[BaseModel]:
        ...

    @property
    @abstractmethod
    def _output_model(self) -> Type[BaseModel]:
        ...

    def run(self):
        with utils.chdir(self.directory):
            # self._fetch_linked_files()
            self._run()

    @abstractmethod
    def _run(self):
        ...

    # def _fetch_linked_files(self):
    #    """Link all files provided in self._linked files

    #    This function is called in run() immediately before the process is executed
    #    """

    #    for dest_filename, (src_process, src_filename) in self._linked_files.items():
    #       if src_process is None:
    #          src_filename = src_filename.resolve()
    #       else:
    #          src_filename = (src_process.directory / src_filename).resolve()
    #
    #       if not src_filename.exists():
    #          raise FileNotFoundError(f'Tried to link {src_filename} with a process, but it does not exist')
    #
    #       dest_filename = (self.directory / dest_filename).resolve()

    #       if not dest_filename.exists():
    #          dest_filename.parent.mkdir(parents=True, exist_ok=True)
    #          utils.symlink(src_filename, dest_filename)

    def __repr__(self):
        return f'{self.__class__.__name__}(inputs={self.inputs.__dict__}, outputs={self.outputs.__dict__})'
