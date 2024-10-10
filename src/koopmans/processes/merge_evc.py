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
        return ' '.join([f'merge_evc.x -nr {np.prod(self.inputs.kgrid)}']
                        + [f'-i {f.aspath()}' for f in self.inputs.src_files]
                        + [f'-o {self.inputs.dest_filename}'])

    def _pre_run(self):
        pass

    def _post_run(self):
        self.outputs = MergeEVCOutputs(merged_file=FilePointer(self, self.inputs.dest_filename))
