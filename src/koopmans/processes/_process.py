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
        self.name: str = self.__class__.__name__.lower().replace('process', '')
        self.directory: Path | None = None

    @property
    @abstractmethod
    def _input_model(self) -> Type[BaseModel]:
        ...

    @property
    @abstractmethod
    def _output_model(self) -> Type[BaseModel]:
        ...

    def run(self):
        assert self.directory is not None, 'Process directory must be set before running'
        with utils.chdir(self.directory):
            # self._fetch_linked_files()
            self._run()

    @abstractmethod
    def _run(self):
        ...

    def __repr__(self):
        return f'{self.__class__.__name__}(inputs={self.inputs.__dict__}, outputs={self.outputs.__dict__})'

    def todict(self):
        utils.warn(
            f'Serialization of {self.__class__.__name__} is not implemented (would require rewriting the serialization module to support circular references)')
        return {}
