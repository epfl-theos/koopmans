from pathlib import Path
from typing import List, Tuple
from pydantic import ConfigDict

from koopmans import utils
from koopmans.files import File
from koopmans.process_io import IOModel

from ._process import Process


class ConvertFilesFromSpin2To1InputModel(IOModel):
    spin_2_files: List[File]
    spin_1_files: List[Path]
    model_config = ConfigDict(arbitrary_types_allowed=True)


class ConvertFilesOutputModel(IOModel):
    generated_files: List[File]
    model_config = ConfigDict(arbitrary_types_allowed=True)


class ConvertFilesFromSpin2To1(Process[ConvertFilesFromSpin2To1InputModel, ConvertFilesOutputModel]):
    input_model = ConvertFilesFromSpin2To1InputModel  # type: ignore
    output_model = ConvertFilesOutputModel  # type: ignore

    def _run(self):

        generated_files = []
        for spin_2_file, spin_1_file in zip(self.inputs.spin_2_files, self.inputs.spin_1_files):
            contents = spin_2_file.read_bytes()

            contents = contents.replace(b'nk="2"', b'nk="1"')
            contents = contents.replace(b'nspin="2"', b'nspin="1"')

            dst = File(self, spin_1_file)
            dst.write_bytes(contents)
            generated_files.append(dst)

        self.outputs = self.output_model(generated_files=generated_files)


class ConvertFilesFromSpin1To2InputModel(IOModel):
    spin_1_files: List[File]
    spin_2_up_files: List[Path]
    spin_2_down_files: List[Path]
    model_config = ConfigDict(arbitrary_types_allowed=True)


class ConvertFilesFromSpin1To2(Process[ConvertFilesFromSpin1To2InputModel, ConvertFilesOutputModel]):

    input_model = ConvertFilesFromSpin1To2InputModel
    output_model = ConvertFilesOutputModel

    def _run(self):

        generated_files = []

        for spin_1_file, spin_2_up_file, spin_2_down_file in zip(self.inputs.spin_1_files,
                                                                 self.inputs.spin_2_up_files,
                                                                 self.inputs.spin_2_down_files, strict=True):

            if spin_2_up_file.is_absolute() or spin_2_down_file.is_absolute():
                # Prevent writing to an absolute path because this means that the process can affect files
                # outside of its working directory
                raise ValueError(
                    f'{self.__class__.__name__} is attempting to write to a file outside of its working directory; this is not allowed')

            contents = spin_1_file.read_bytes()

            contents = contents.replace(b'nk="1"', b'nk="2"')
            contents = contents.replace(b'nspin="1"', b'nspin="2"')

            spin_up_file = File(self, spin_2_up_file)
            spin_up_file.write_bytes(contents)

            contents = contents.replace(b'ik="1"', b'ik="2"')
            contents = contents.replace(b'ispin="1"', b'ispin="2"')

            spin_down_file = File(self, spin_2_down_file)
            spin_down_file.write_bytes(contents)

            generated_files += [spin_up_file, spin_down_file]

        self.outputs = self.output_model(generated_files=generated_files)


class SwapSpinFilesInputModel(IOModel):
    read_directory: File
    model_config = ConfigDict(arbitrary_types_allowed=True)


class SwapSpinFilesOutputModel(IOModel):
    write_directory: File
    model_config = ConfigDict(arbitrary_types_allowed=True)


class SwapSpinFilesProcess(Process[SwapSpinFilesInputModel, SwapSpinFilesOutputModel]):

    input_model = SwapSpinFilesInputModel
    output_model = SwapSpinFilesOutputModel

    def _run(self):
        if not self.inputs.read_directory.exists():
            raise FileNotFoundError(f'{self.inputs.read_directory} does not exist')
        spin_up_files = list(self.inputs.read_directory.rglob('*1.*'))
        spin_down_files = list(self.inputs.read_directory.rglob('*2.*'))
        for src in self.inputs.read_directory.rglob('*'):
            if src.is_dir():
                continue

            if src.parent.exists():
                src.parent.mkdir(parents=True, exist_ok=True)

            if src in spin_up_files:
                dst = File(self, str(src.name).replace('1.', '2.'))
            elif src in spin_down_files:
                dst = File(self, str(src.name).replace('2.', '1.'))
            else:
                dst = File(self, src.name)

            try:
                dst.symlink_to(src)
            except FileExistsError:
                assert src == dst

        self.outputs = self.output_model(write_directory=self.inputs.read_directory.name)
