"""

Contains the default settings for KI/KIPZ calculations

Written by Edward Linscott Jan 2020

"""

from koopmans.io import warn
from koopmans.calculators.cp import CP_calc
from koopmans.calculators.pw import PW_calc
from koopmans.calculators.wannier90 import W90_calc
from koopmans.calculators.pw2wannier import PW2Wannier_calc

defaults = {CP_calc: {'calculation': 'cp',
                      'iprint': 1,
                      'outdir': '../TMP-CP/',
                      'prefix': 'kc',
                      'verbosity': 'low',
                      'disk_io': 'high',
                      'do_ee': True,
                      'do_wf_cmplx': True,
                      'electron_dynamics': 'cg',
                      'nspin': 2,
                      'ortho_para': 1,
                      'passop': 2.0,
                      'ion_dynamics': 'none',
                      'ion_nstepe': 5,
                      'ion_radius(1)': 1.0,
                      'ion_radius(2)': 1.0,
                      'ion_radius(3)': 1.0,
                      'ion_radius(4)': 1.0,
                      'do_innerloop_cg': True,
                      'innerloop_cg_nreset': 20,
                      'innerloop_cg_nsd': 2,
                      'innerloop_init_n': 3,
                      'innerloop_nmax': 100,
                      'hartree_only_sic': False,
                      'conv_thr': '1.0e-9*nelec',
                      'esic_conv_thr': '1.0e-9*nelec'}}


def load_defaults(calc):
    calc_class = calc.__class__

    if calc_class not in defaults:
        warn(f'No defaults found for {calc_class} calculator')
        return

    for key, value in defaults[calc_class].items():
        if getattr(calc, key, value) not in [None, value]:
            # If a setting has already been set, keep that value but print a warning
            warn(
                f'Suggested value for {key} is being overwritten; do this with caution')
        else:
            setattr(calc, key, value)
    return
