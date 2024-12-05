"""Defines Process, a class for abstractly representing a Workflow/Calculator/other operation.

Inspired by CWL."""

import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Generic, Tuple, Type, TypeVar
from uuid import uuid4

import dill
import numpy as np
from pydantic import BaseModel

from koopmans import utils

if TYPE_CHECKING:
    from koopmans.workflows import Workflow


class IOModel(BaseModel):
    # BaseModel with an __eq__ method that compares attributes rather than memory addresses

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        for key in self.dict().keys():
            cond = (getattr(self, key) == getattr(other, key))
            if isinstance(cond, np.ndarray):
                cond = cond.all()
            if not cond:
                return False
        return True


InputModel = TypeVar('InputModel', bound=IOModel)
OutputModel = TypeVar('OutputModel', bound=IOModel)


class Process(utils.HasDirectory, ABC, Generic[InputModel, OutputModel]):

    __slots__ = utils.HasDirectory.__slots__ + ['inputs', '_outputs', 'name']

    def __init__(self, name: str | None = None, **kwargs):
        self.inputs: InputModel = self.input_model(**kwargs)
        self._outputs: OutputModel | None = None
        if name is None:
            name_with_split_acronyms = re.sub(r'([a-z])([A-Z])', r'\1_\2',
                                              self.__class__.__name__.replace('Process', ''))
            self.name = re.sub(r'([A-Z])([A-Z][a-z])', r'\1_\2', name_with_split_acronyms).lower()
        else:
            self.name = name

        # Initialize the directory information
        super().__init__()

    def run(self):
        self._pre_run()
        self._run()
        self._post_run()

    @property
    @abstractmethod
    def input_model(self) -> Type[InputModel]:
        ...

    @property
    @abstractmethod
    def output_model(self) -> Type[OutputModel]:
        ...

    @property
    def outputs(self) -> OutputModel:
        if self._outputs is None:
            raise ValueError('Process has no outputs because it has not been run yet')
        return self._outputs

    @outputs.setter
    def outputs(self, value: OutputModel):
        self._outputs = value

    def _pre_run(self):
        assert self.directory is not None, 'Process directory must be set before running'
        self.dump_inputs()
        if self.engine is None:
            raise ValueError('Process engine must be set before running')

    @abstractmethod
    def _run(self):
        ...

    def _post_run(self):
        assert self.outputs is not None, 'Process outputs must be set when running'
        self.dump_outputs()

    def __repr__(self):
        out = f'{self.__class__.__name__}(inputs={self.inputs.__dict__}'
        if self._outputs is not None:
            out += f', outputs={self.outputs.__dict__}'
        return out + ')'

    def dump_inputs(self, directory: Path | None = None):
        dst = Path(f'{self.name}_inputs.pkl')
        if directory is not None:
            dst = directory / dst
        with open(dst, 'wb') as f:
            dill.dump(self.inputs, f)

    def dump_outputs(self, directory: Path | None = None):
        dst = Path(f'{self.name}_outputs.pkl')
        if directory is not None:
            dst = directory / dst
        with open(dst, 'wb') as f:
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
