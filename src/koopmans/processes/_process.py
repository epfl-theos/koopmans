"""Defines Process, a class for abstractly representing a Workflow/Calculator/other operation.

Inspired by the Common Workflow Language (CWL).
"""

import logging
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Generic, Type, TypeVar

import dill

from koopmans import utils
from koopmans.files import File
from koopmans.process_io import IOModel

InputModel = TypeVar('InputModel', bound=IOModel)
OutputModel = TypeVar('OutputModel', bound=IOModel)


class Process(utils.HasDirectory, ABC, Generic[InputModel, OutputModel]):
    """Class that defines a Process.

    This could be a Workflow, a Commandline tool, a calculator, etc.
    """

    __slots__ = utils.HasDirectory.__slots__ + ['inputs', '_outputs',
                                                'name', 'linked_files']

    input_model: Type[InputModel]
    output_model: Type[OutputModel]

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

        self.linked_files: Dict[str, File] = {}

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
        logger.info(f'Querying outputs of {self.directory}')
        return self._outputs

    @outputs.setter
    def outputs(self, value: OutputModel):
        logger = logging.getLogger(__name__)
        logger.info(f'Setting outputs of {self.directory}')
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
