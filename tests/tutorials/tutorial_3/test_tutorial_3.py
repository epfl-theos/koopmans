import pytest

from koopmans import base_directory
from koopmans.io import read, write_json
from koopmans.utils import chdir

tutorial_dir = base_directory / 'tutorials' / 'tutorial_3'

w90_override = {
    "calculator_parameters": {
        "w90": {
            "occ": {
                "projections_blocks": [
                    [{"site": "Zn", "ang_mtm": "l=0"}],
                    [{"site": "Zn", "ang_mtm": "l=1"}],
                    [{"site": "O", "ang_mtm": "l=0"}],
                    [{"site": "Zn", "ang_mtm": "l=2"},
                     {"site": "O", "ang_mtm": "l=1"}]
                ]},
            "emp": {
                "dis_froz_max": 14.5,
                "dis_win_max": 17.0,
                "projections": [
                    {"site": "Zn", "ang_mtm": "l=0"}
                ]
            }
        }
    }
}


@pytest.mark.tutorials
def test_zno_pdos(tutorial_patch):
    with chdir(tutorial_dir):
        wf = read('zno.json')
        wf.run()


@pytest.mark.tutorials
def test_replot_dft_bandstructure():
    with chdir(tutorial_dir):
        exec(open('replot_dft_bandstructure.py', 'r').read())


@pytest.mark.tutorials
def test_zno_wannierize(tutorial_patch):
    with chdir(tutorial_dir):
        wf = read('zno.json', override={'workflow': {'task': 'wannierize'}, **w90_override})
        wf.run()


@pytest.mark.tutorials
def test_zno_ki(tutorial_patch):
    with chdir(tutorial_dir):
        wf = read('zno.json', override={'workflow': {'task': 'singlepoint'}, **w90_override})
        write_json(wf, 'test.json')
        wf.run()


@pytest.mark.tutorials
def test_zno_plot_bandstructures():
    with chdir(tutorial_dir):
        exec(open('plot_bandstructures.py', 'r').read())
