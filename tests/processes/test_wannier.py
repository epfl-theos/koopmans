"""Testing the `koopmans.processes.wannier` module."""

from pathlib import Path
from typing import List

import numpy as np
import pytest  # noqa: F401

from koopmans.engines.localhost import LocalhostEngine
from koopmans.files import LocalFile
from koopmans.processes.wannier import (ExtendProcess, MergeProcess,
                                        extend_wannier_u_dis_file_content,
                                        merge_wannier_hr_file_contents)
from koopmans.utils import (chdir, parse_wannier_hr_file_contents,
                            parse_wannier_u_file_contents,
                            read_wannier_hr_file)


def test_wannierize_merge_hr_file_contents(tmp_path, datadir):
    """Test merging wannier_hr.dat files."""
    with chdir(tmp_path):
        filecontents = []
        for dir_in in sorted((datadir / 'w90').glob('occ_block*')):
            with open(dir_in / 'wann_hr.dat') as f:
                filecontents.append(f.read())

        merged_contents = merge_wannier_hr_file_contents(filecontents)

        ham, rvec, weights, num_wann = parse_wannier_hr_file_contents(merged_contents)

        ham_ref, rvec_ref, weights_ref, num_wann_ref = read_wannier_hr_file(
            LocalFile(datadir / 'w90' / 'occ' / 'wann_hr.dat'))

        assert np.allclose(ham, ham_ref)
        assert np.all(rvec == rvec_ref)
        assert weights == weights_ref
        assert num_wann == num_wann_ref


def test_merge_process(tmp_path):
    """Test MergeProcess."""
    file1_contents = 'a\nb'
    file2_contents = 'c\nd'
    engine = LocalhostEngine()

    with chdir(tmp_path):
        file1 = LocalFile('file1.txt')
        file2 = LocalFile('file2.txt')
        file3 = LocalFile('file3.txt')

        engine.write_file(file1_contents, file1)
        engine.write_file(file2_contents, file2)

        def merge_fn(all_contents: List[str]) -> str:
            return "\n".join(all_contents)

        process = MergeProcess(src_files=[file1, file2], dst_file=file3.name, merge_function=merge_fn)
        process.directory = Path()
        process.engine = engine
        process.run()

        assert file3.read_text() == file1_contents + '\n' + file2_contents


def test_extend_wannier_u_dis_file_content(tmp_path, datadir):
    """Test extending a wannier_u_dis.mat file."""
    with open(datadir / 'w90' / 'wannier90_u_dis.mat') as f:
        filecontent = f.read()

    nbnd = 10
    nwann = 8

    with chdir(tmp_path):
        contents = extend_wannier_u_dis_file_content(filecontent, nbnd=nbnd, nwann=nwann)
        umat, _, nk = parse_wannier_u_file_contents(contents)

        assert umat.shape == (nk, nwann, nbnd)


def test_extend_process(tmp_path):
    """Test ExtendProcess."""
    file1_contents = 'a\nb'
    extra_contents = 'c\nd'
    engine = LocalhostEngine()

    with chdir(tmp_path):
        file1 = LocalFile('file1.txt')
        file2 = LocalFile('file2.txt')

        engine.write_file(file1_contents, file1)

        def extend_fn(contents) -> str:
            return contents + '\n' + extra_contents

        process = ExtendProcess(src_file=file1, dst_file=file2.name, extend_function=extend_fn)
        process.directory = Path()
        process.engine = engine
        process.run()

        assert file2.read_text() == file1_contents + '\n' + extra_contents
