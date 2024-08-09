from pathlib import Path
from typing import List, Tuple

from koopmans import utils
from koopmans.files import FilePointer
from koopmans.outputs import OutputModel

from ._process import Process


class ConvertFilesFromSpin2To1InputModel(OutputModel):
    spin_2_files: List[FilePointer]
    spin_1_files: List[Path]

    class Config:
        arbitrary_types_allowed = True


class ConvertFilesOutputModel(OutputModel):
    generated_files: List[Path]


class ConvertFilesFromSpin2To1(Process):
    _input_model = ConvertFilesFromSpin2To1InputModel  # type: ignore
    _output_model = ConvertFilesOutputModel  # type: ignore

    def _run(self):

        for spin_2_file, spin_1_file in zip(self.inputs.spin_2_files, self.inputs.spin_1_files):
            contents = utils.get_binary_content(*spin_2_file)

            contents = contents.replace(b'nk="2"', b'nk="1"')
            contents = contents.replace(b'nspin="2"', b'nspin="1"')

            utils.write_binary_content(spin_1_file, contents)

        self.outputs = self._output_model(generated_files=self.inputs.spin_1_files)


class ConvertFilesFromSpin1To2InputModel(OutputModel):
    spin_1_files: List[FilePointer]
    spin_2_up_files: List[Path]
    spin_2_down_files: List[Path]

    class Config:
        arbitrary_types_allowed = True


class ConvertFilesFromSpin1To2(Process):
    _input_model = ConvertFilesFromSpin1To2InputModel  # type: ignore
    _output_model = ConvertFilesOutputModel  # type: ignore

    def _run(self):

        for spin_1_file, spin_2_up_file, spin_2_down_file in zip(self.inputs.spin_1_files,
                                                                 self.inputs.spin_2_up_files,
                                                                 self.inputs.spin_2_down_files):

            contents = utils.get_binary_content(*spin_1_file)

            contents = contents.replace(b'nk="1"', b'nk="2"')
            contents = contents.replace(b'nspin="1"', b'nspin="2"')

            utils.write_binary_content(spin_2_up_file, contents)

            contents = contents.replace(b'ik="1"', b'ik="2"')
            contents = contents.replace(b'ispin="1"', b'ispin="2"')

            utils.write_binary_content(spin_2_down_file, contents)

        self.outputs = self._output_model(generated_files=self.inputs.spin_2_up_files + self.inputs.spin_2_down_files)
