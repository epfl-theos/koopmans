from pathlib import Path

import numpy as np
import pytest

from koopmans.utils import (chdir, parse_wannier_hr_file_contents,
                            read_wannier_hr_file)
from koopmans.workflows._merge_wannier import merge_wannier_hr_file_contents


def test_wannierize_merge_hr_file_contents(tmp_path, datadir):
    with chdir(tmp_path):
        filecontents = []
        for dir_in in sorted((datadir / 'w90').glob('occ_block*')):
            with open(dir_in / 'wann_hr.dat') as f:
                filecontents.append(f.readlines())

        merged_contents = merge_wannier_hr_file_contents(filecontents)

        ham, rvec, weights, num_wann = parse_wannier_hr_file_contents(merged_contents)

        ham_ref, rvec_ref, weights_ref, num_wann_ref = read_wannier_hr_file(datadir / 'w90' / 'occ' / 'wann_hr.dat')

        assert np.allclose(ham, ham_ref)
        assert np.all(rvec == rvec_ref)
        assert weights == weights_ref
        assert num_wann == num_wann_ref
