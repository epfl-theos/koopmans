'''

Tests for koopmans.io._calculators

'''

import pytest

from koopmans.calculators import (Calc, KoopmansCPCalculator,
                                  KoopmansHamCalculator,
                                  KoopmansScreenCalculator, PhCalculator,
                                  PW2WannierCalculator, PWCalculator,
                                  UnfoldAndInterpolateCalculator,
                                  Wann2KCCalculator, Wann2KCPCalculator,
                                  Wannier90Calculator)
from koopmans.io._calculators import read_calculator

example_files = {KoopmansCPCalculator: 'kcp/ki_final',
                 KoopmansHamCalculator: ['kcw/kc.khi', 'kcw/kc.kho'],
                 KoopmansScreenCalculator: ['kcw/kc.ksi', 'kcw/kc.kso'],
                 PW2WannierCalculator: 'pw2wannier/pw2wan',
                 PhCalculator: 'ph/eps',
                 PWCalculator: 'pw/scf',
                 UnfoldAndInterpolateCalculator: 'ui/ki',
                 Wann2KCCalculator: ['kcw/kc.w2ki', 'kcw/kc.w2ko'],
                 Wann2KCPCalculator: 'wann2kcp/w2kcp',
                 Wannier90Calculator: 'w90/wann'}


@pytest.mark.parametrize('calc_class,files', example_files.items())
def test_read_calculator(calc_class, files, datadir):
    if isinstance(files, str):
        path = datadir / files
    else:
        path = [datadir / x for x in files]

    calc: Calc = read_calculator(path)
    assert isinstance(calc, calc_class)
