from pathlib import Path

from koopmans.files import File
from koopmans.io import read_pkl, write_pkl
from koopmans.utils import chdir
from koopmans.workflows import SinglepointWorkflow


def test_file_reduce(water, tmp_path):
    with chdir(tmp_path):
        wf = SinglepointWorkflow(**water)
        wf.directory = Path()
        calc = wf.new_calculator('kcp', outdir='tmp', nspin=2, nelec=8, ndw=50, prefix='test_read_ham')
        calc.directory = Path()
        f = File(parent=calc, name='test.txt')

        write_pkl(f, 'test.pkl')
        new_f = read_pkl('test.pkl')

        assert f == new_f
