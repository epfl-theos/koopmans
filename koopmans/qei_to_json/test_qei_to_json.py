from pathlib import Path
from koopmans.qei_to_json import qei_to_json


def test_pwi_to_json():
    cwd = Path(__file__).parent
    qei_to_json(cwd / 'example.pwi', cwd / 'example.json')


def test_cpi_to_json():
    cwd = Path(__file__).parent
    qei_to_json(cwd / 'example.cpi', cwd / 'example.json')
