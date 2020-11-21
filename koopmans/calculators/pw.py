"""

pw calculator module for python_KI

Written by Edward Linscott Sep 2020

"""

from ase.io import espresso_cp as pw_io
from ase.units import create_units
from koopmans.calculators.calculator import QE_calc


class PW_calc(QE_calc):
    # Subclass of QE_calc for performing calculations with pw.x

    # Point to the appropriate ASE IO module
    _io = pw_io

    # Define the appropriate file extensions
    ext_in = '.pwi'
    ext_out = '.pwo'

    # Adding all pw.x keywords as decorated properties of the PW_calc class.
    # This means one can set and get pw.x keywords as self.<keyword> but
    # internally they are stored as self._settings['keyword'] rather than
    # self.<keyword>
    _recognised_keywords = []

    for keywords in _io.KEYS.values():
        for k in keywords:
            _recognised_keywords.append(k)

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

    def is_complete(self):
        return self.results.get('job done', False)

    def is_converged(self):
        return self.results.get('energy', None) is not None
