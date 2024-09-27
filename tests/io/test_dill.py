from pathlib import Path

import pytest

from koopmans.io import read_pkl, write_pkl
from koopmans.utils import chdir

base_directory = Path('/home/user/base_directory')
paths = [base_directory / 'subpath/to/directory',
         Path('/home/user/path/that/is/not/a/subpath/of/base_directory/'), Path('relative/path')]


@pytest.mark.parametrize('obj', paths)
def test_roundtrip_pkl_with_placeholder(obj, tmp_path):
    with chdir(tmp_path):
        write_pkl(obj, 'test.pkl', base_directory=base_directory)

        new_obj = read_pkl('test.pkl', base_directory=base_directory)

        assert obj == new_obj


@pytest.mark.parametrize('obj', paths)
def test_roundtrip_pkl_without_placeholder(obj, tmp_path):
    with chdir(tmp_path):
        write_pkl(obj, 'test.pkl')

        new_obj = read_pkl('test.pkl')

        assert obj == new_obj
