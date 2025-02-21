from pathlib import Path
from typing import List

from koopmans import utils
from koopmans.commands import Command
from koopmans.files import File

from ._commandlinetool import CommandLineTool
from ._process import IOModel


class Bin2XMLInput(IOModel):
    binary: File

    class Config:
        arbitrary_types_allowed = True


class Bin2XMLOutput(IOModel):
    xml: File

    class Config:
        arbitrary_types_allowed = True


class Bin2XMLProcess(CommandLineTool):

    input_model = Bin2XMLInput
    output_model = Bin2XMLOutput

    def _pre_run(self):
        super()._pre_run()
        if not self.inputs.binary.parent_process.exists():
            raise FileNotFoundError(f'`{self.inputs.binary}` does not exist')

        # Link the input binary file to the directory of this process as input.dat
        dst = File(self, Path("input.dat"))
        dst.symlink_to(self.inputs.binary)

    @property
    def command(self):
        return Command(executable='bin2xml.x', suffix=f'input.dat output.xml')

    def _set_outputs(self):
        xml_filepointer = File(self, Path("output.xml"))
        self.outputs = self.output_model(xml=xml_filepointer)
