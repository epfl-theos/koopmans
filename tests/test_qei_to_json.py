"""Test koopmans.qei_to_json."""

from koopmans import utils
from koopmans.qei_to_json import qei_to_json


def test_pwi_to_json(tmpdir, datadir):
    """Test the conversion of .pwi to .json."""
    with utils.chdir(tmpdir):
        qei_to_json(datadir / 'qei_to_json' / 'example.pwi', 'example.json')


def test_cpi_to_json(tmpdir, datadir):
    """Test the conversion of .cpi to .json."""
    with utils.chdir(tmpdir):
        qei_to_json(datadir / 'qei_to_json' / 'example.cpi', 'example.json')
