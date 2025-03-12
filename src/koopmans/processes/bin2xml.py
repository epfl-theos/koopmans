from pathlib import Path
from typing import List
from pydantic import ConfigDict

from koopmans import utils
from koopmans.commands import Command
from koopmans.files import File
from koopmans.process_io import IOModel

from ._commandlinetool import CommandLineTool


class Bin2XMLInput(IOModel):
    binary: File
    model_config = ConfigDict(arbitrary_types_allowed=True)


class Bin2XMLOutput(IOModel):
    xml: File
    model_config = ConfigDict(arbitrary_types_allowed=True)


class Bin2XMLProcess(CommandLineTool):

    input_model = Bin2XMLInput
    output_model = Bin2XMLOutput

    def _pre_run(self):
        super()._pre_run()
        if not self.inputs.binary.exists():
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
