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
                raise ValueError(f'Debug: {key} differs: {getattr(self, key)} != {getattr(other, key)}')
                return False
        return True


InputModel = TypeVar('InputModel', bound=IOModel)
OutputModel = TypeVar('OutputModel', bound=IOModel)


class Process(ABC, Generic[InputModel, OutputModel]):

    __slots__ = ['inputs', '_outputs', 'name', 'directory', 'parent', 'uuid', '_base_directory']

    def __init__(self, name: str | None = None, **kwargs):
        self.inputs: InputModel = self.input_model(**kwargs)
        self._outputs: OutputModel | None = None
        if name is None:
            name_with_split_acroynms = re.sub(r'([a-z])([A-Z])', r'\1_\2',
                                              self.__class__.__name__.replace('Process', ''))
            self.name = re.sub(r'([A-Z])([A-Z][a-z])', r'\1_\2', name_with_split_acroynms).lower()
        else:
            self.name = name
        self.directory: Path | None = None
        self.parent: Workflow | None = None
        self.uuid = str(uuid4())
        self._base_directory: Path | None = None
        if self.parent is None:
            self.base_directory = Path()

    def run(self):
        assert self.directory is not None, 'Process directory must be set before running'
        with utils.chdir(self.directory):
            self.dump_inputs()
            self._run()
            assert self.outputs is not None, 'Process outputs must be set when running'
            self.dump_outputs()

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

    @abstractmethod
    def _run(self):
        ...

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

    @property
    def base_directory(self) -> Path:
        if self.parent is not None:
            return self.parent.base_directory
        else:
            if self._base_directory is None:
                raise ValueError(f'{self.__class__.__name__}.base_directory has not been set')
            return self._base_directory

    @base_directory.setter
    def base_directory(self, value: Path):
        if self.parent is not None:
            raise ValueError(f'{self.__class__.__name__}.base_directory should not be set for processes with parents')
        self._base_directory = value.resolve()

    @property
    def absolute_directory(self) -> Path:
        assert self.directory is not None
        if self.parent is None:
            return self.base_directory / self.directory
        path = self.parent.directory / self.directory

        # Recursive through the parents, adding their directories to path (these are all relative paths)...
        obj = self.parent
        while getattr(obj, 'parent', None):
            assert obj.parent is not None
            path = obj.parent.directory / path
            obj = obj.parent

        # ... until we reach the top-level parent, which should have a base_directory attribute (an absolute path)
        if not hasattr(obj, 'base_directory'):
            raise AttributeError(f'Expected `{obj.__class__.__name__}` instance to have a `base_directory` attribute')
        return obj.base_directory / path
