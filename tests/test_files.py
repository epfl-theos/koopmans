from pathlib import Path

from koopmans.calculators import KoopmansCPCalculator
from koopmans.files import AbsoluteFilePointer, FilePointer, ParentPlaceholder
from koopmans.io import read_pkl, write_pkl
from koopmans.utils import chdir
from koopmans.workflows import SinglepointWorkflow


def test_absolute_file_pointer():
    path = Path('/tmp/test.txt')
    abs_fp = AbsoluteFilePointer(path)
    assert abs_fp.aspath() == path
    assert abs_fp.parent.absolute_directory == path.parent


def test_filepointer_reduce(water, tmp_path):
    with chdir(tmp_path):
        wf = SinglepointWorkflow(**water)
        wf.directory = Path()
        calc = KoopmansCPCalculator(outdir='tmp', nspin=2, nelec=8, ndw=50, prefix='test_read_ham', **water)
        calc.parent = wf
        calc.directory = Path()
        fp = FilePointer(parent=calc, name='test.txt')

        write_pkl(fp, 'test.pkl')
        new_fp = read_pkl('test.pkl')

        assert fp == new_fp
