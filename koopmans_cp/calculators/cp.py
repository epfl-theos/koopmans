"""

cp calculator module for python_KI

Written by Edward Linscott Sep 2020

"""

from ase.io import espresso_cp as cp_io
from ase.units import create_units
from koopmans_cp.calc import QE_calc

# Quantum ESPRESSO uses CODATA 2006 internally
units = create_units('2006')

class CP_calc(QE_calc):
    # Subclass of QE_calc for performing calculations with cp.x
     
    # Point to the appropriate ASE IO module
    _io = cp_io

    # Define the appropriate file extensions
    ext_in = '.cpi'
    ext_out = '.cpo'

    # Adding all cp.x keywords as decorated properties of the CP_calc class.
    # This means one can set and get cp.x keywords as self.<keyword> but
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
        return self.results['job_done']

    def is_converged(self):
        # Checks convergence of the calculation
        if self.conv_thr is None:
            raise ValueError('Cannot check convergence when "conv_thr" is not set')
        return self._ase_is_converged()

    def _ase_is_converged(self):
        if 'convergence' not in self.results:
            raise ValueError('Could not locate calculation details to check convergence')

        # Check convergence for both filled and empty, allowing for the possibility
        # of do_outerloop(_empty) = False meaning the calculation is immpediately
        # 'converged'
        do_outers = [self.do_outerloop, self.do_outerloop_empty]
        convergence_data = self.results['convergence'].values()
        converged = []
        for do_outer, convergence in zip(do_outers, convergence_data):
            if not do_outer:
                converged.append(True)
            else:
                converged.append(convergence[-1]['delta_E'] < self.conv_thr*units.Hartree)
        return all(converged)
