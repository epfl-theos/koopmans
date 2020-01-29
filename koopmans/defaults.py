"""

Contains the default settings for KI/KIPZ calculations

Written by Edward Linscott Jan 2020

"""

defaults = {'calculation':         'cp',
            'iprint':              1,
            'outdir':              './TMP-CP/',
            'prefix':              'pbe',
            'verbosity':           'low',
            'disk_io':             'high',
            'do_ee':               True,
            'do_wf_cmplx':         True,
            'ecutrho':             260.0,
            'ecutwfc':             65.0,
            'ibrav':               0,
            'electron_dynamics':   'cg',
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
