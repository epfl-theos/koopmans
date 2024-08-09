from pathlib import Path
from typing import List

from pydantic import BaseModel

from koopmans import utils
from koopmans.commands import Command
from koopmans.files import FilePointer

from ._commandlinetool import CommandLineTool


class Bin2XMLInput(BaseModel):
    binary: FilePointer

    class Config:
        arbitrary_types_allowed = True


class Bin2XMLOutput(BaseModel):
    xml: FilePointer

    class Config:
        arbitrary_types_allowed = True


class Bin2XMLProcess(CommandLineTool):

    _input_model = Bin2XMLInput
    _output_model = Bin2XMLOutput

    def _pre_run(self):
        if not (self.inputs.binary.parent.directory / self.inputs.binary.name).exists():
            raise FileNotFoundError(f'{self.inputs.binary} does not exist')

        # Load the file
        binary_file_contents = utils.get_binary_content(*self.inputs.binary)

        # Write the file to disk
        utils.write_binary_content(Path("input.dat"), binary_file_contents)

    @property
    def command(self):
        return Command(executable='bin2xml.x', suffix=f'input.dat output.xml')

    def _post_run(self):
        xml_filepointer = FilePointer(self, "output.xml")
        self.outputs = self._output_model(xml=xml_filepointer)
