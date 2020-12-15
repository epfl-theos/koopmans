"""

pw2wannier calculator module for python_KI

Written by Edward Linscott Sep 2020

"""

from koopmans.utils import warn
from ase.io import pw2wannier as p2w_io
from ase.calculators.pw2wannier import PW2Wannier
from koopmans.calculators.generic import QE_calc


class PW2Wannier_calc(QE_calc):
    # Link to relevant ase io module
    _io = p2w_io

    # Define the appropriate file extensions
    ext_in = '.p2wi'
    ext_out = '.p2wo'

    # Adding all wannier90 keywords as decorated properties of the W90_calc class.
    # This means one can set and get wannier90 keywords as self.<keyword> but
    # internally they are stored as self._settings['keyword'] rather than
    # self.<keyword>
    _recognised_keywords = ['outdir', 'prefix', 'seedname', 'write_mmn',
                            'write_amn', 'write_uHu', 'uHu_formatted', 'write_unk', 'reduce_unk',
                            'wan_mode', 'wannier_plot', 'wannier_plot_list', 'split_evc_file']

    for k in _recognised_keywords:
        # We need to use these make_get/set functions so that get/set_k are
        # evaluated immediately (otherwise we run into late binding and 'k'
        # is not defined when get/set_k are called)
        def make_get(key):
            def get_k(self):
                # Return 'None' rather than an error if the keyword has not
                # been defined
                return self._settings.get(key, None)
            return get_k

        def make_set(key):
            def set_k(self, value):
                self._settings[key] = value
            return set_k

        get_k = make_get(k)
        set_k = make_set(k)
        locals()[k] = property(get_k, set_k)

    def __init__(self, *args, **kwargs):
        self._ase_calc_class = PW2Wannier
        self.settings_to_not_parse = ['wan_mode']
        super().__init__(*args, **kwargs)

    def _update_settings_dict(self):
        if 'inputpp' not in self._ase_calc.parameters:
            self._ase_calc.parameters['inputpp'] = {}
        self._settings = self._ase_calc.parameters['inputpp']

    def is_converged(self):
        warn("is_converged is not properly implemented")
        return True

    def is_complete(self):
        return self.results['job done']
