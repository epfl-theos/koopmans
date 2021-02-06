"""
Unfolding & interpolating calculator for python_KI

Originally written by Riccardo De Gennaro as the standalone 'unfolding and interpolate' code
Integrated within python_KI by Edward Linscott Jan 2021

"""

import os
from time import time
from ase.atoms import Atoms
from ase.calculators.calculator import FileIOCalculator
from koopmans.calculators.generic import QE_calc
from koopmans import utils


class UI_calc(QE_calc):

    from ._io import parse_w90, parse_hr, parse_phases, print_centers, write_results, write_bands, write_dos, \
        write_input_file, read_input_file, read_output_file
    from ._interpolate import interpolate, calc_bands, correct_phase, calc_dos
    from ._settings import load_defaults, valid_settings

    ext_in = '.uii'
    ext_out = '.uio'

    # Adding all UI keywords as decorated properties of the W90_calc class.
    # This means one can set and get wannier90 keywords as self.<keyword> but
    # internally they are stored as self._settings['keyword'] rather than
    # self.<keyword>
    _recognised_keywords = [s.name for s in valid_settings]

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

    # Subclass of QE_calc for performing unfolding and interpolation
    def __init__(self, calc=None, qe_files=[], skip_qc=False, **kwargs):
        self._recognised_keywords = []
        self._ase_calc_class = None
        self.settings_to_not_parse = ['w90_seedname', 'kc_ham_file', 'dft_ham_file', 'dft_smooth_ham_file']

        super().__init__(calc, qe_files, skip_qc, **kwargs)

        self.results_for_qc = []

    def calculate(self):

        for setting in self.mandatory_settings:
            assert getattr(self, setting, None), f'The mandatory key {setting} has not been provided'

        if self.name is None:
            self.name = 'ui'

        if self.directory is None:
            self.directory = '.'

        self.write_input_file()

        start = time()
        reset = time()

        with utils.chdir(self.directory):
            with open(f'{self.name}{self.ext_out}', 'w') as f_out:
                self.f_out = f_out

                self.f_out.write('\nUNFOLDING & INTERPOLATION\n\n')

                """
                 1) Parse data:
                    - calc parameters from the JSON file
                    - other parameters from W90 output file
                    - hamiltonian(s)
                    - WFs phases read from file wf_phases.dat
                """

                self.parse_w90()
                self.parse_hr()
                self.parse_phases()

                self.f_out.write(f'\tParsing input in:{time() - reset:25.3f} sec\n')
                reset = time()

                """
                 2) Core of the unfolding and interpolation code:
                    - build the map |i> ---> |Rn>
                    - calc interpolated (if needed) bands
                    - calc DOS (if needed)
                """

                self.interpolate(reset)

                reset = time()

                """
                 3) Print out the results:
                    - bands into 'bands_interpolated.dat' file
                    - DOS into 'dos_interpolated.dat' file
                """

                self.write_results()

                self.f_out.write(f'\tPrinting output in: {time() - reset:24.3f} sec\n')

                self.f_out.write(f'\n\tTotal time: {time() - start:32.3f} sec\n')
                self.f_out.write('\nALL DONE\n\n')

                self.calc.results['job done'] = True

    def check_code_is_installed(self):
        # This calculator is entirely python-based, so we don't need to check for an installed binary
        return True

    def is_complete(self):
        return self.calc.results.get('job done', False)

    @property
    def preprocessing_flags(self):
        return None

    @preprocessing_flags.setter
    def preprocessing_flags(self, value):
        if value != ['']:
            raise ValueError('You should not try and set the preprocessing flags for the UI calculator')

    @property
    def do_smooth_interpolation(self):
        return any([f > 1 for f in self.smooth_int_factor])
