"""Process for running `bin2xml.x`."""

from pathlib import Path

from koopmans.files import File
from koopmans.process_io import IOModel

from ._commandlinetool import CommandLineTool


class Bin2XMLInput(IOModel):
    """Input for a `Bin2XMLProcess`."""

    binary: File


class Bin2XMLOutput(IOModel):
    """Output for a `Bin2XMLProcess`."""

    xml: File


class Bin2XMLProcess(CommandLineTool):
    """Process for running `bin2xml.x`."""

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
        """Return the command that this process runs."""
        return 'bin2xml.x input.dat output.xml'

    def _set_outputs(self):
        xml_filepointer = File(self, Path("output.xml"))
        self.outputs = self.output_model(xml=xml_filepointer)
