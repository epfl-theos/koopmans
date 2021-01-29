"""

Contains the default settings for KI/KIPZ calculations

Written by Edward Linscott Jan 2020

"""

from koopmans.utils import warn
from koopmans.calculators.kcp import KCP_calc
from koopmans.calculators.pw import PW_calc
from koopmans.calculators.wannier90 import W90_calc
from koopmans.calculators.pw2wannier import PW2Wannier_calc

defaults = {KCP_calc: {
            'calculation': 'cp',
            'outdir': './TMP-CP/',
            'iprint': 1,
            'prefix': 'kc',
            'verbosity': 'low',
            'disk_io': 'high',
            'write_hr': False,
            'do_ee': True,
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
            'esic_conv_thr': '1.0e-9*nelec'},
            PW_calc: {
            'calculation': 'scf',
            'outdir': './TMP-PW/',
            'prefix': 'kc',
            'conv_thr': '2.0e-9*nelec'},
            W90_calc: {
            'num_iter': 10000,
            'conv_tol': 1.e-10,
            'conv_window': 5,
            'write_hr': True,
            'guiding_centres': True,
            'gamma_only': False},
            PW2Wannier_calc: {
            'outdir': '../TMP-PW/',
            'prefix': 'kc',
            'seedname': 'wann',
            'wan_mode': 'standalone'}}
