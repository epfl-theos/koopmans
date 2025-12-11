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
from koopmans.pydantic import BaseModel

if TYPE_CHECKING:
    from koopmans.files import FileABC

InputModel = TypeVar('InputModel', bound=BaseModel)
OutputModel = TypeVar('OutputModel', bound=BaseModel)


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
    _directory: Path | None = PrivateAttr(None)
    parent_process: Process | None = None
    engine: Engine
    _directory_must_be_relative: bool = False

    @model_validator(mode='before')
    @classmethod
    def base_directory_only_for_parentless_process(cls, data: dict[str, Any]) -> dict[str, Any]:
        if get(data, 'parent_process', None) is not None:
            set_(data, 'base_directory', get(data, 'parent_process.base_directory'))
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
        logger.info(f'with inputs={self.inputs.model_dump()}')

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
        from koopmans.files import File, LinkedFile

        self.inputs.model_config['frozen'] = True

        assert self.directory is not None, 'Process directory must be set before running'
        self.dump_inputs()

        # Copy over any inputs that correspond to files
        files_to_copy = self.inputs.fields_that_are_files()
        for inp in files_to_copy:
            dest_file = File(parent_process=self, relative_path=inp.dest)
            if inp.symlink:
                dest_file.symlink_to(inp.src, overwrite=inp.overwrite, recursive=inp.recursive_symlink)
            else:
                inp.src.copy_to(dest_file, exist_ok=inp.overwrite)

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
        from koopmans.files import File
        assert self.directory is not None
        dst = File(parent_process=self, relative_path=f'{self.name}_inputs.json')
        dst.write_text(self.inputs.model_dump_json(indent=2, exclude_defaults=True, exclude_none=True))

    def dump_outputs(self):
        """Dump the outputs of the process to the file system."""
        from koopmans.files import File
        assert self.directory is not None
        dst = File(self, f'{self.name}_outputs.pkl')
        dst.write_bytes(dill.dumps(self.outputs))

    def load_outputs(self):
        """Load the outputs of the process from the file system."""
        from koopmans.files import File
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
        from koopmans.files import File
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
