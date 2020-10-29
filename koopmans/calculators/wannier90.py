"""

wannier90 calculator module for python_KI

Written by Edward Linscott Sep 2020

"""

from ase.io import wannier90 as w90_io
from koopmans.calculators.calculator import QE_calc

class W90_calc(QE_calc):
    # Link to relevant ase io module
    _io = w90_io

    # Define the appropriate file extensions
    ext_in = '.win'

    # Adding all wannier90 keywords as decorated properties of the W90_calc class.
    # This means one can set and get wannier90 keywords as self.<keyword> but
    # internally they are stored as self._settings['keyword'] rather than 
    # self.<keyword>
    _recognised_keywords = ['atoms_frac', 'bands_plot', 'dis_froz_max', 
        'dis_num_iter', 'dis_win_max', 'guiding_centres', 'kpoint_path', 'kpoints',
        'mp_grid', 'num_bands', 'num_iter', 'num_print_cycles', 'num_wann', 
        'projections', 'unit_cell_cart']

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

    def _update_settings_dict(self):
       self._settings = self._ase_calc.parameters

    def calculate(self):
       self._ase_calc.calculate()

