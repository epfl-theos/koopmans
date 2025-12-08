"""Defines Process, a class for abstractly representing a Workflow/Calculator/other operation.

Inspired by the Common Workflow Language (CWL).
"""

from __future__ import annotations

import logging
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, ClassVar, Dict, Generic, TypeVar

import dill
from pydantic import PrivateAttr, field_validator, model_validator

from koopmans.engines.engine import Engine
from koopmans.files import File
from koopmans.pydantic import BaseModel

InputModel = TypeVar('InputModel', bound=BaseModel)
OutputModel = TypeVar('OutputModel', bound=BaseModel)


class Process(BaseModel, ABC, Generic[InputModel, OutputModel]):
    """Class that defines a Process.

    This could be a Workflow, a Commandline tool, a calculator, etc.
    """

    inputs: InputModel
    input_model: ClassVar[type[InputModel]]
    output_model: ClassVar[type[OutputModel]]
    name: str
    _outputs: OutputModel | None = PrivateAttr(None)
    _linked_files: Dict[str, File] = PrivateAttr(default_factory=dict)
    _base_directory: Path | None = PrivateAttr(None)
    _directory: Path | None = PrivateAttr(None)
    engine: Engine | None = None
    parent_process: Process | None = None
    _directory_must_be_relative: bool = False

    def __init_subclass__(cls, **kwargs):
        """Define `input_model` and `output_model` automatically from the Generic specialization.

        i.e. this will define
        cls.input_model = InputModel
        cls.output_model = OutputModel
        """
        super().__init_subclass__(**kwargs)
        # detect Generic specialization
        for base in getattr(cls, "__orig_bases__", []):
            origin = getattr(base, "__origin__", None)
            if origin is Process:
                in_type, out_type = base.__args__
                cls.input_model = in_type
                cls.output_model = out_type
                import ipdb
                ipdb.set_trace()

    @model_validator(mode='before')
    @classmethod
    def base_directory_only_for_parentless_process(cls, data: dict[str, Any]) -> dict[str, Any]:
        if data.get('parent_process', None) is not None:
            data.pop('base_directory', None)
        return data

    @field_validator('name', mode='before')
    @classmethod
    def default_name(cls, name: str | None) -> str | None:
        """Generate a default name if none is provided."""
        if name is None:
            name_with_split_acronyms = re.sub(r'([a-z])([A-Z])', r'\1_\2',
                                              cls.__name__.replace('Process', ''))
            return re.sub(r'([A-Z])([A-Z][a-z])', r'\1_\2', name_with_split_acronyms).lower()
        return name

    def model_post_init(self, context: Any) -> None:
        super().model_post_init(context)
        logger = logging.getLogger(__name__)
        logger.info(f'Creating an instance of {self.__class__.__name__}')
        logger.info(f'with inputs={self.inputs.model_dump()}')

    def run(self):
        """Run the process."""
        logger = logging.getLogger(__name__)
        logger.info(f'Running {self.directory}')
        self._pre_run()
        self._run()
        self._post_run()

    @property
    def outputs(self) -> OutputModel:
        """Return the outputs of the process."""
        if self._outputs is None:
            raise ValueError('Process has no outputs because it has not been run yet')
        logger = logging.getLogger(__name__)
        logger.debug(f'Querying outputs of {self.directory}', stacklevel=2)
        return self._outputs

    @outputs.setter
    def outputs(self, value: OutputModel):
        logger = logging.getLogger(__name__)
        logger.info(f'Setting outputs of {self.directory}', stacklevel=2)
        self._outputs = value

    def _pre_run(self):
        assert self.directory is not None, 'Process directory must be set before running'
        self.dump_inputs()

        if self.engine is None:
            raise ValueError('Process engine must be set before running')

        # Link the files in self.linked_files
        for dest, src in self.linked_files.items():
            dest_file = File(self, dest)
            dest_file.symlink_to(src)

    @abstractmethod
    def _run(self):
        ...

    def _post_run(self):
        assert self._outputs is not None, 'Process outputs must be set when running'
        self.dump_outputs()

    def __repr__(self):
        out = f'{self.__class__.__name__}(inputs={self.inputs.__dict__}'
        if self._outputs is not None:
            out += f', outputs={self.outputs.__dict__}'
        return out + ')'

    def dump_inputs(self):
        """Dump the inputs of the process to the file system."""
        assert self.directory is not None
        dst = File(self, f'{self.name}_inputs.pkl')
        dst.write_bytes(dill.dumps(self.inputs))

    def dump_outputs(self):
        """Dump the outputs of the process to the file system."""
        assert self.directory is not None
        dst = File(self, f'{self.name}_outputs.pkl')
        dst.write_bytes(dill.dumps(self.outputs))

    def load_outputs(self):
        """Load the outputs of the process from the file system."""
        if self.directory is None:
            raise ValueError('Process directory must be set before attempting to load outputs')
        src = File(self, f'{self.name}_outputs.pkl')
        content = src.read_bytes()
        self.outputs = dill.loads(content)

    def is_complete(self):
        """Return True if the process has been run and the outputs have been saved."""
        if self.directory is None:
            raise ValueError('Process directory must be set before checking if it is complete')
        return (self.directory / f'{self.name}_outputs.pkl').exists()

    def __truediv__(self, other):
        assert isinstance(other, Path) or isinstance(other, str)
        return File(self, other)

    @property
    def directory(self) -> Path | None:
        """Return the directory of the process."""
        return self._directory

    @directory.setter
    def directory(self, value: Path | str | None):
        if value is None:
            return

        # Sanitize input
        if isinstance(value, str):
            value = Path(value)

        # Sanity checks
        if value.is_absolute() and self._directory_must_be_relative:
            raise ValueError(
                f'{self.__class__.__name__} directory must be a relative path '
                f'(relative to {self.__class__.__name__},base directory)'
            )

        self._directory = value

    @property
    def uid(self) -> str:
        """Return the UID of the process."""
        return str(self.directory)

    @property
    def base_directory(self) -> Path | None:
        """Return the base (i.e. root) directory of the entire workflow."""
        if self.parent_process:
            return self.parent_process.base_directory
        else:
            return self._base_directory

    @base_directory.setter
    def base_directory(self, value: Path | str | None):
        if self.parent_process is not None:
            raise ValueError('Do not directly set `base_directory` for objects with a parent_process')
        if isinstance(value, str):
            value = Path(value)
        self._base_directory = None if value is None else value.resolve()

    @property
    def absolute_directory(self) -> Path | None:
        """Return the absolute path of the directory."""
        if self.directory is None:
            return None
        if self.base_directory is None:
            return None
        return self.base_directory / self.directory

    def directory_has_been_set(self) -> bool:
        """Check if the directory has been set."""
        return self._directory is not None


class ProcessPlaceholder(Process[BaseModel, BaseModel]):

    @field_validator('parent_process', mode='after')
    @classmethod
    def must_be_none(cls, value: Any) -> None:
        if value is not None:
            raise ValueError('ProcessPlaceholder parent_process must be None')
        return None
