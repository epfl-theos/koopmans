from pathlib import Path

import numpy as np
import pytest
from ase_koopmans.dft.kpoints import BandPath
from ase_koopmans.io.wannier90 import read_wannier90_out

from koopmans.engines import LocalhostEngine
from koopmans.files import LocalFile
from koopmans.processes import ui
from koopmans.settings import (PlotSettingsDict,
                               UnfoldAndInterpolateSettingsDict)
from koopmans.utils import chdir
from tests.helpers.patches import benchmark_filename


def test_ui_si(silicon, tmp_path, datadir, check_patch):
    with chdir(tmp_path):
        atoms = silicon['atoms']
        parameters = UnfoldAndInterpolateSettingsDict(smooth_int_factor=2, w90_calc='pc', do_map=True, use_ws_distance=True,
                                                      do_dos=True, kgrid=[2, 2, 2], kpath=atoms.cell.bandpath().interpolate('GL', density=10))

        # Load the centers and spreads from file
        with open(datadir / "ui" / "wann.wout") as f:
            w90_calc = next(read_wannier90_out(f)).calc

        proc = ui.UnfoldAndInterpolateProcess(atoms=atoms,
                                              parameters=parameters,
                                              centers=np.array(w90_calc.results['centers']),
                                              spreads=w90_calc.results['spreads'],
                                              kc_ham_file=LocalFile(datadir / "ui" / "kc_ham.dat"),
                                              dft_ham_file=LocalFile(datadir / "ui" / "dft_ham.dat"),
                                              dft_smooth_ham_file=LocalFile(
                                                  datadir / "ui" / "smooth_dft_ham.dat"),
                                              plotting_parameters=PlotSettingsDict(degauss=0.05, nstep=1000, Emin=-10, Emax=4))
        proc.directory = Path()
        proc.engine = LocalhostEngine()
        proc.run()
