"""

pw calculator module for python_KI

Written by Edward Linscott Sep 2020

"""

from ase.io import espresso as pw_io
from ase.calculators.espresso import Espresso
from ase.units import create_units
from koopmans import io
from koopmans.calculators.generic import QE_calc


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

    def __init__(self, *args, **kwargs):
        self._ase_calc_class = Espresso
        self.settings_to_not_parse = ['pseudo_dir', 'assume_isolated']
        super().__init__(*args, **kwargs)
        self.results_for_qc = ['energy']

    @property
    def calc(self):
        # First, update the param block
        self._ase_calc.parameters['input_data'] = self.construct_namelist()

        return self._ase_calc

    @calc.setter
    def calc(self, value):
        self._ase_calc = value

    def is_complete(self):
        return self.results.get('job done', False)

    def is_converged(self):
        return self.results.get('energy', None) is not None

    def _ase_calculate(self):
        # Before running the calculation, update the keywords for the ASE calculator object
        self._ase_calc.parameters['input_data'] = self.construct_namelist()
        super()._ase_calculate()

    def construct_namelist(self):
        # Returns a namelist of settings, grouped by their Quantum Espresso headings
        return pw_io.construct_namelist(**self._settings, warn=True)

    def _update_settings_dict(self):
        # Updates self._settings based on self._ase_calc
        self._settings = {}
        for namelist in self._ase_calc.parameters.get('input_data', {}).values():
            for key, val in namelist.items():
                self._settings[key] = val

    @property
    def nelec(self):
        # PW does not directly have an "nelec" keyword, but it is useful for setting
        # extensive properties such as energy convergence thresholds
        return io.nelec_from_pseudos(self.calc)

    @nelec.setter
    def nelec(self, value):
        raise ValueError('You tried to set PW_calc.nelec, which is invalid')
