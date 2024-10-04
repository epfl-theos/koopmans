from pathlib import Path

import pytest

from koopmans import utils
from koopmans.files import AbsoluteFilePointer
from koopmans.processes.koopmans_cp import (ConvertFilesFromSpin1To2,
                                            ConvertFilesFromSpin2To1,
                                            SwapSpinFilesProcess)


def test_convert_files_from_spin_1_to_2_process(tmpdir, check_patch, datadir):
    with utils.chdir(tmpdir):
        spin_1_files = sorted([AbsoluteFilePointer(f) for f in (datadir / 'kcp' / 'spin_1_files').glob('*')])
        spin_2_up_files = sorted((datadir / 'kcp' / 'spin_2_files').glob('*1.???'))
        spin_2_down_files = sorted((datadir / 'kcp' / 'spin_2_files').glob('*2.???'))

        process = ConvertFilesFromSpin1To2(spin_1_files=spin_1_files,
                                           spin_2_up_files=spin_2_up_files,
                                           spin_2_down_files=spin_2_down_files)
        process.directory = Path()
        process.run()


def test_convert_files_from_spin_2_to_1_process(tmpdir, check_patch, datadir):
    with utils.chdir(tmpdir):
        spin_1_files = sorted((datadir / 'kcp' / 'spin_1_files').glob('*'))
        spin_2_files = sorted([AbsoluteFilePointer(f) for f in (datadir / 'kcp' / 'spin_2_files').glob('*1.???')])

        process = ConvertFilesFromSpin2To1(spin_2_files=spin_2_files,
                                           spin_1_files=spin_1_files)
        process.directory = Path()
        process.run()


def test_swap_spin_files(tmpdir, check_patch, datadir):
    with utils.chdir(tmpdir):
        process = SwapSpinFilesProcess(read_directory=AbsoluteFilePointer(datadir / 'kcp' / 'spin_2_files'))
        process.directory = Path()
        process.run()
