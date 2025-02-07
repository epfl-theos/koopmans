from pathlib import Path
from typing import List

import numpy as np

from koopmans.files import FilePointer
from koopmans.kpoints import Kpoints
from koopmans.projections import ProjectionBlock

from ._commandlinetool import CommandLineTool
from ._process import IOModel


class MergeEVCInputs(IOModel):
    kgrid: List[int]
    src_files: List[FilePointer]
    dest_filename: str

    class Config:
        arbitrary_types_allowed = True


class MergeEVCOutputs(IOModel):
    merged_file: FilePointer

    class Config:
        arbitrary_types_allowed = True


class MergeEVCProcess(CommandLineTool):

    input_model = MergeEVCInputs
    output_model = MergeEVCOutputs

    @property
    def command(self):
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
        self.outputs = MergeEVCOutputs(merged_file=FilePointer(self, self.inputs.dest_filename))
