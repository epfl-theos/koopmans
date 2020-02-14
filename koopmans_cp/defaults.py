"""

Contains the default settings for KI/KIPZ calculations

Written by Edward Linscott Jan 2020

"""

from koopmans_cp.io import warn

defaults = {'calculation':         'cp',
            'iprint':              1,
            'outdir':              './TMP-CP/',
            'prefix':              'kc',
            'verbosity':           'low',
            'disk_io':             'high',
            'do_ee':               True,
            'do_wf_cmplx':         True,
            'electron_dynamics':   'cg',
            'nspin':               2,
            'ortho_para':          1,
            'passop':              2.0,
            'which_compensation':  'tcc',
            'ion_dynamics':        'none',
            'ion_nstepe':          5,
            'ion_radius(1)':       1.0,
            'ion_radius(2)':       1.0,
            'ion_radius(3)':       1.0,
            'ion_radius(4)':       1.0,
            'do_innerloop_cg':     True,
            'innerloop_cg_nreset': 7,
            'innerloop_cg_nsd':    2,
            'innerloop_init_n':    3,
            'innerloop_nmax':      300,
            'hartree_only_sic':    False}


def load_defaults(calc):
    for key, value in defaults.items():
        if getattr(calc, key, None) is not None:
            # If a setting has already been set, keep that value but print a warning
            warn(f'Default value for {key} is being overwritten')
        else:
            setattr(calc, key, value)
    return calc
