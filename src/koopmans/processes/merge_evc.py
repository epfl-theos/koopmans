"""Process for running merge_evc.x."""

from typing import List

import numpy as np
from pydantic import ConfigDict

from koopmans.files import File

from ._commandlinetool import CommandLineTool
from ._process import IOModel


class MergeEVCInputs(IOModel):
    """Input model for a `MergeEVCProcess`."""

    kgrid: List[int]
    src_files: List[File]
    dest_filename: str
    model_config = ConfigDict(arbitrary_types_allowed=True)


class MergeEVCOutputs(IOModel):
    """Output model for a `MergeEVCProcess`."""

    merged_file: File
    model_config = ConfigDict(arbitrary_types_allowed=True)


class MergeEVCProcess(CommandLineTool):
    """Commandline tool for running merge_evc.x."""

    input_model = MergeEVCInputs
    output_model = MergeEVCOutputs

    @property
    def command(self):
        """Return the command that this process runs."""
        input_files = [f'input_{i}.dat' for i in range(len(self.inputs.src_files))]
        return ' '.join([f'merge_evc.x -nr {np.prod(self.inputs.kgrid)}']
                        + [f'-i {f}' for f in input_files]
                        + [f'-o {self.inputs.dest_filename}'])

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Add src_files to those that we need to link
        for i, src_file in enumerate(self.inputs.src_files):
            self.linked_files[f'input_{i}.dat'] = src_file

    def _set_outputs(self):
        self.outputs = MergeEVCOutputs(merged_file=File(self, self.inputs.dest_filename))
