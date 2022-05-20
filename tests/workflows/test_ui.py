import pytest
from koopmans import workflows
from koopmans.utils import chdir


def test_ui_si(silicon, tmp_path, sys2file, datadir):
    with chdir(tmp_path):
        parameters = {'ui': {'smooth_int_factor': 2,
                             "kc_ham_file": datadir / "ui" / "kc_ham.dat",
                             "w90_seedname": datadir / "ui" / "wann",
                             "dft_ham_file": datadir / "ui" / "dft_ham.dat",
                             "dft_smooth_ham_file": datadir / "ui" / "smooth_dft_ham.dat",
                             "w90_calc": "pc",
                             "do_map": True,
                             "use_ws_distance": True,
                             "do_dos": True
                             },
                      "plot": {"degauss": 0.05,
                               "nstep": 1000,
                               "Emin": -10,
                               "Emax": 4
                               }
                      }
        wf = workflows.SingleUnfoldAndInterpolateWorkflow(atoms=silicon['atoms'],
                                                          kgrid=[2, 2, 2],
                                                          kpath="GL",
                                                          master_calc_params=parameters)
        wf.run()
