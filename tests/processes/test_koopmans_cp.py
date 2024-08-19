import pytest

from koopmans import utils


def test_convert_files_from_spin_1_to_2_process(tmpdir, check_patch):
    from koopmans.processes.koopmans_cp import ConvertFilesFromSpin1To2

    with utils.chdir(tmpdir):

        raise NotImplementedError("Need to get singlepoint test up and running to generate files to use in this test")


def test_convert_files_from_spin_2_to_1_process(tmpdir, check_patch):
    from koopmans.processes.koopmans_cp import ConvertFilesFromSpin2To1

    with utils.chdir(tmpdir):

        raise NotImplementedError("Need to get singlepoint test up and running to generate files to use in this test")
