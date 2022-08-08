import pytest

from ase.spectrum.band_structure import BandStructure
from koopmans import workflows
from koopmans.io import read_kwf as read_encoded_json
from koopmans.testing import benchmark_filename
from koopmans.utils import chdir


def test_ui_si(silicon, tmp_path, sys2file, datadir, ui_patch):
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
                                                          kpoints=silicon['kpoints'],
                                                          name='si_ui',
                                                          calculator_parameters=parameters)
        wf.kpoints.set_path('GL', silicon['atoms'].cell)
        wf.run()

        # Check the results. Usually result-checking is automated within the MockCalc class, but here...
        # 1) we don't want to mock the calculators
        # 2) we can afford to be a lot more strict, because it's just python we're running so the results should be
        #    identical
        calc = wf.calculations[-1]
        results = calc.results
        # Save the calculator as an encoded json in the benchmarks directory
        with open(benchmark_filename(calc), 'r') as fd:
            calc_ref = read_encoded_json(fd)
        for key, result in results.items():
            # Don't compare walltime
            if key == 'walltime':
                continue

            assert result == calc_ref.results[key]
