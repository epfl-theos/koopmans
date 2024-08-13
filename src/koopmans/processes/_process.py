"""Defines Process, a class for abstractly representing a Workflow/Calculator/other operation.

Inspired by CWL."""

import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Tuple, Type

import dill
from pydantic import BaseModel

from koopmans import calculators, utils


class Process(ABC):

    __slots__ = ['inputs', 'outputs', 'name', 'directory']

    def __init__(self, name=None, **kwargs):
        self.inputs = self._input_model(**kwargs)
        self.outputs = None
        if name is None:
            name_with_split_acroynms = re.sub(r'([a-z])([A-Z])', r'\1_\2',
                                              self.__class__.__name__.replace('Process', ''))
            self.name = re.sub(r'([A-Z])([A-Z][a-z])', r'\1_\2', name_with_split_acroynms).lower()
        else:
            self.name = name
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
            self.dump_inputs()
            self._run()
            assert self.outputs is not None, 'Process outputs must be set when running'
            self.dump_outputs()

    @abstractmethod
    def _run(self):
        ...

    def __repr__(self):
        return f'{self.__class__.__name__}(inputs={self.inputs.__dict__}, outputs={self.outputs.__dict__})'

    def dump_inputs(self):
        with open(f'{self.name}_inputs.pkl', 'wb') as f:
            dill.dump(self.inputs, f)

    def dump_outputs(self):
        with open(f'{self.name}_outputs.pkl', 'wb') as f:
            dill.dump(self.outputs, f)

    def load_outputs(self):
        if self.directory is None:
            raise ValueError('Process directory must be set before attempting to load outputs')
        with open(self.directory / f'{self.name}_outputs.pkl', 'rb') as f:
            self.outputs = dill.load(f)

    def is_complete(self):
        if self.directory is None:
            raise ValueError('Process directory must be set before checking if it is complete')
        return (self.directory / f'{self.name}_outputs.pkl').exists()
