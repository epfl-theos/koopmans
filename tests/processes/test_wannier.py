from pathlib import Path
from typing import List

import numpy as np
import pytest

from koopmans.files import AbsoluteFilePointer
from koopmans.processes.wannier import (ExtendProcess, MergeProcess,
                                        extend_wannier_u_dis_file_content,
                                        merge_wannier_hr_file_contents)
from koopmans.utils import (chdir, parse_wannier_hr_file_contents,
                            parse_wannier_u_file_contents,
                            read_wannier_hr_file, write_content)


def test_wannierize_merge_hr_file_contents(tmp_path, datadir):
    with chdir(tmp_path):
        filecontents = []
        for dir_in in sorted((datadir / 'w90').glob('occ_block*')):
            with open(dir_in / 'wann_hr.dat') as f:
                filecontents.append(f.readlines())

        merged_contents = merge_wannier_hr_file_contents(filecontents)

        ham, rvec, weights, num_wann = parse_wannier_hr_file_contents(merged_contents)

        ham_ref, rvec_ref, weights_ref, num_wann_ref = read_wannier_hr_file(
            (None, datadir / 'w90' / 'occ' / 'wann_hr.dat'))

        assert np.allclose(ham, ham_ref)
        assert np.all(rvec == rvec_ref)
        assert weights == weights_ref
        assert num_wann == num_wann_ref


def test_merge_process(tmp_path):
    file1_contents = ['a', 'b']
    file2_contents = ['c', 'd']

    with chdir(tmp_path):
        file1 = AbsoluteFilePointer('file1.txt')
        file2 = AbsoluteFilePointer('file2.txt')
        file3 = AbsoluteFilePointer('file3.txt')

        write_content(file1.name, file1_contents)
        write_content(file2.name, file2_contents)

        def merge_fn(all_contents: List[List[str]]) -> List[str]:
            return [c for contents in all_contents for c in contents]

        process = MergeProcess(src_files=[file1, file2], dst_file=file3.name, merge_function=merge_fn)
        process.directory = Path()
        process.run()

        assert file3.read() == file1_contents + file2_contents


def test_extend_wannier_u_dis_file_content(tmp_path, datadir):
    with open(datadir / 'w90' / 'wannier90_u_dis.mat') as f:
        filecontent = f.readlines()

    nbnd = 10
    nwann = 8

    with chdir(tmp_path):
        contents = extend_wannier_u_dis_file_content(filecontent, nbnd=nbnd, nwann=nwann)
        umat, _, nk = parse_wannier_u_file_contents(contents)

        assert umat.shape == (nk, nwann, nbnd)


def test_extend_process(tmp_path):
    file1_contents = ['a', 'b']
    extra_contents = ['c', 'd']

    with chdir(tmp_path):
        file1 = AbsoluteFilePointer('file1.txt')
        file2 = AbsoluteFilePointer('file2.txt')

        write_content(file1.name, file1_contents)

        def extend_fn(contents) -> List[str]:
            return contents + extra_contents

        process = ExtendProcess(src_file=file1, dst_file=file2.name, extend_function=extend_fn)
        process.directory = Path()
        process.run()

        assert file2.read() == file1_contents + extra_contents
