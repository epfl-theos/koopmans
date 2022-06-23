from pathlib import Path
from koopmans import utils
from koopmans.qei_to_json import qei_to_json


def test_pwi_to_json(tmpdir, datadir):
    with utils.chdir(tmpdir):
        basedir = datadir / 'qei_to_json'
        qei_to_json(basedir / 'example.pwi', basedir / 'example.json')


def test_cpi_to_json(tmpdir, datadir):
    with utils.chdir(tmpdir):
        basedir = datadir / 'qei_to_json'
        qei_to_json(basedir / 'example.cpi', basedir / 'example.json')
