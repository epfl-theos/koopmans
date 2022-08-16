import pytest

from koopmans.calculators import UnfoldAndInterpolateCalculator


def test_read_input(silicon, datadir):
    calc = UnfoldAndInterpolateCalculator(atoms=silicon['atoms'], directory=datadir / 'ui')
    calc.prefix = 'ki'
    calc.read_input()

    # Verify a few miscellaneous keywords
    assert calc.parameters.smooth_int_factor == [2, 2, 2]
    assert calc.parameters.do_dos is True
    assert calc.parameters.kgrid == [2, 2, 2]
    assert calc.parameters.plotting.Emin == -10


def test_read_results(silicon, datadir):
    calc = UnfoldAndInterpolateCalculator(atoms=silicon['atoms'], directory=datadir / 'ui')
    calc.prefix = 'ki'
    calc.read_results()

    # Verify the results (.uio files only confirm if the job is done)
    assert 'job done' in calc.results
    assert calc.results['job done']
