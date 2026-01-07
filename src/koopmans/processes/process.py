"""Defines Process, a class for abstractly representing a Workflow/Calculator/other operation.

Inspired by the Common Workflow Language (CWL).
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Generic, TypeVar

import dill
from pydantic import Field, PrivateAttr, field_validator, model_validator
from pydash import get, set_

from koopmans.engines.engine import Engine
from koopmans.files import File, LocalFile
from koopmans.pydantic import BaseModel

InputModel = TypeVar('InputModel', bound=BaseModel)
OutputModel = TypeVar('OutputModel', bound=BaseModel)


def auto_name_from_class(cls: type) -> str:
    name_with_split_acronyms = re.sub(r'([a-z])([A-Z])', r'\1_\2',
                                      cls.__name__.replace('Process', ''))
    return re.sub(r'([A-Z])([A-Z][a-z])', r'\1_\2', name_with_split_acronyms).lower()


class Process(BaseModel, Generic[InputModel, OutputModel]):
    """Class that defines a Process.

    This could be a Workflow, a Commandline tool, a calculator, etc.
    """

    inputs: InputModel
    input_model: ClassVar[type[InputModel]]
    output_model: ClassVar[type[OutputModel]]
    name: str
    _outputs: OutputModel | None = PrivateAttr(None)
    base_directory: Path = Field(default_factory=lambda: Path.cwd().resolve())
    directory: Path | None = None
    parent_process_uid: str | None = None
    engine: Engine

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        # Set a class-specific default_factory
        cls.model_fields["name"].default_factory = (
            lambda cls=cls: auto_name_from_class(cls)
        )

    @property
    def uid(self) -> str:
        """Return the UID of the process."""
        if self.directory is not None:
            return str(self.directory)
        else:
            return self.name

    @model_validator(mode='before')
    @classmethod
    def base_directory_only_for_parentless_process(cls, data: dict[str, Any]) -> dict[str, Any]:
        if get(data, 'parent_process', None) is not None:
            set_(data, 'base_directory', get(data, 'parent_process.base_directory'))
        return data

    @field_validator('inputs', mode='before')
    @classmethod
    def coerce_inputs(cls, value: Any) -> InputModel:
        """Force early coercion of inputs to the correct type."""
        try:
            input_model_cls = cls.input_model
        except AttributeError:
            raise ValueError(f'Forgot to set {cls.__name__}.input_model. Please report this bug.')

        if isinstance(value, input_model_cls):
            return value
        elif isinstance(value, BaseModel):
            return input_model_cls.model_validate(value.model_dump())
        return cls.input_model.model_validate(value)

    def model_post_init(self, context: Any) -> None:
        super().model_post_init(context)
        logger = logging.getLogger(__name__)
        logger.info(f'Creating an instance of {self.__class__.__name__}')

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
        logger.info(f'Setting outputs of {self.directory}', stacklevel=5)
        self._outputs = value

    def _pre_run(self):
        # Check directory is set
        assert self.directory is not None, 'Process directory must be set before running'

        # Dump inputs
        self.dump_inputs()

        # Copy over any inputs that correspond to files
        files_to_copy = self.inputs.fields_that_are_files()
        for inp in files_to_copy:
            self._attach_file(inp.src, inp.dest, inp.symlink, inp.overwrite, inp.recursive_symlink)

    def _attach_file(self, src: File | LocalFile, dest: Path, symlink: bool, overwrite: bool, recursive_symlink: bool) -> None:
        """Attach a file to this process.

        This involves either copying or symlinking the file to the process directory."""
        dest_file = File(parent_process_uid=self.uid, path=dest)
        if symlink:
            self.engine.link_file(src, dest_file, overwrite=overwrite, recursive=recursive_symlink)
        else:
            self.engine.copy_file(src, dest_file, exist_ok=overwrite)

    def _post_run(self):
        assert self._outputs is not None, 'Process outputs must be set when running'
        self.dump_outputs()

    def __repr__(self):
        out = f'{self.__class__.__name__}(inputs={self.inputs.__dict__}'
        if self._outputs is not None:
            out += f', outputs={self.outputs.__dict__}'
        return out + ')'

    def write_to_file(self, file: File | LocalFile | Path | str, content: str | bytes) -> None:
        """Write text to a file."""
        self.engine.write_to_file(sanitize_file_argument(file, parent_process_uid=self.uid), content)

    def read_from_file(self, file: File | LocalFile | Path | str) -> str:
        """Read text from a file."""
        return self.engine.read_file(sanitize_file_argument(file, parent_process_uid=self.uid))

    def dump_inputs(self):
        """Dump the inputs of the process to the file system."""
        self.write_to_file(f'inputs.json', self.inputs.model_dump_json(
            indent=2, exclude_defaults=True, exclude_none=True))

    def dump_outputs(self):
        """Dump the outputs of the process to the file system."""
        self.write_to_file(f'outputs.json', self.outputs.model_dump_json(
            indent=2, exclude_none=True))

    def load_outputs(self):
        """Load the outputs of the process from the file system."""
        if self.directory is None:
            raise ValueError('Process directory must be set before attempting to load outputs')
        self.outputs = self.output_model.model_validate_json(self.read_from_file(f'outputs.json'))

    def is_complete(self):
        """Return True if the process has been run and the outputs have been saved."""
        if self.directory is None:
            raise ValueError('Process directory must be set before checking if it is complete')
        return (self.directory / 'outputs.json').exists()

    def __truediv__(self, other):
        from koopmans.files import File
        assert isinstance(other, Path) or isinstance(other, str)
        return File(self, other)

    @field_validator('directory', mode='after')
    def check_directory_is_relative(cls, value: Path | None) -> Path | None:
        if value is not None and value.is_absolute():
            raise ValueError(
                f'{self.__class__.__name__} directory must be a relative path '
                f'(relative to {self.__class__.__name__}.base directory)'
            )
        return value

    @property
    def absolute_directory(self) -> Path | None:
        """Return the absolute path of the directory."""
        if self.directory is None:
            return None
        if self.base_directory is None:
            return None
        return self.base_directory / self.directory


def sanitize_file_argument(file: File | LocalFile | Path | str, parent_process_uid: str | None = None) -> File | LocalFile:
    """Sanitize a file argument to ensure it is a File or LocalFile."""
    if isinstance(file, str):
        file = Path(file)
    if isinstance(file, Path):
        if Path.is_absolute(file):
            # Interpret absolute paths as local files
            file = LocalFile(path=file)
        else:
            # Interpret relative paths as files relative to the parent process
            if parent_process_uid is None:
                raise ValueError('parent_process_uid must be provided when sanitizing a relative file path')
            file = File(parent_process_uid=parent_process_uid, path=file)
    return file
