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
    input_model = ConvertFilesFromSpin2To1InputModel  # type: ignore
    output_model = ConvertFilesOutputModel  # type: ignore

    def _run(self):

        for spin_2_file, spin_1_file in zip(self.inputs.spin_2_files, self.inputs.spin_1_files):
            contents = utils.get_binary_content(*spin_2_file)

            contents = contents.replace(b'nk="2"', b'nk="1"')
            contents = contents.replace(b'nspin="2"', b'nspin="1"')

            utils.write_binary_content(spin_1_file, contents)

        self.outputs = self.output_model(generated_files=self.inputs.spin_1_files)


class ConvertFilesFromSpin1To2InputModel(OutputModel):
    spin_1_files: List[FilePointer]
    spin_2_up_files: List[Path]
    spin_2_down_files: List[Path]

    class Config:
        arbitrary_types_allowed = True


class ConvertFilesFromSpin1To2(Process):
    input_model = ConvertFilesFromSpin1To2InputModel  # type: ignore
    output_model = ConvertFilesOutputModel  # type: ignore

    def _run(self):

        for spin_1_file, spin_2_up_file, spin_2_down_file in zip(self.inputs.spin_1_files,
                                                                 self.inputs.spin_2_up_files,
                                                                 self.inputs.spin_2_down_files, strict=True):

            if spin_2_up_file.is_absolute() or spin_2_down_file.is_absolute():
                # Prevent writing to an absolute path because this means that the process can affect files
                # outside of its working directory
                raise ValueError(
                    f'{self.__class__.__name__} is attempting to write to a file outside of its working directory; this is not allowed')

            contents = utils.get_binary_content(*spin_1_file)

            contents = contents.replace(b'nk="1"', b'nk="2"')
            contents = contents.replace(b'nspin="1"', b'nspin="2"')

            utils.write_binary_content(spin_2_up_file, contents)

            contents = contents.replace(b'ik="1"', b'ik="2"')
            contents = contents.replace(b'ispin="1"', b'ispin="2"')

            utils.write_binary_content(spin_2_down_file, contents)

        self.outputs = self.output_model(generated_files=self.inputs.spin_2_up_files + self.inputs.spin_2_down_files)


class SwapSpinFilesInputModel(OutputModel):
    read_directory: FilePointer

    class Config:
        arbitrary_types_allowed = True


class SwapSpinFilesOutputModel(OutputModel):
    write_directory: Path


class SwapSpinFilesProcess(Process):

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

            if not src.name.parent.exists():
                src.name.parent.mkdir(parents=True, exist_ok=True)

            if src in spin_up_files:
                dst = Path(str(src.name).replace('1.', '2.'))
            elif src in spin_down_files:
                dst = Path(str(src.name).replace('2.', '1.'))
            else:
                dst = src.name

            try:
                utils.symlink(src.aspath(), dst)
            except FileExistsError:
                assert src.aspath() == dst.resolve()
                pass

        self.outputs = self.output_model(write_directory=self.inputs.read_directory.name)
