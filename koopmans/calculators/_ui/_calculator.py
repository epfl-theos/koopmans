"""

The calculator class defining the Unfolding & interpolating calculator

"""

from typing import List
from time import time
import numpy as np
from .._utils import ExtendedCalculator
from ._settings import valid_settings
from koopmans import utils


class UnfoldAndInterpolateCalculator(ExtendedCalculator):
    # Subclass of ExtendedCalculator for performing unfolding and interpolation
    from ._io import parse_w90, parse_hr, parse_phases, print_centers, write_results, write_bands, write_dos, \
        write_input_file, read_input_file, read_output_file, read_bands
    from ._interpolate import interpolate, calc_bands, correct_phase, calc_dos
    from ._settings import load_defaults, valid_settings

    # UnfoldAndInterpolateCalculator does not use ASE
    _io = None

    ext_in = '.uii'
    ext_out = '.uio'

    # Adding all UI keywords as decorated properties of the UI calc class
    _valid_settings: List[str] = [s.name for s in valid_settings]
    _settings_that_are_paths = ['w90_seedname', 'kc_ham_file', 'dft_ham_file', 'dft_smooth_ham_file']

    def __init__(self, calc=None, qe_files=[], skip_qc=False, **kwargs):
        self._ase_calc_class = None
        super().__init__(calc, qe_files, skip_qc, **kwargs)

        # If we were reading generating this object from files, look for bands, too
        if qe_files and any([self.ext_out in f for f in qe_files]):
            self.read_bands()

        self.results_for_qc = ['band structure', 'dos']

    def calculate(self):

        if self.kpath is None:
            utils.warn('"kpath" missing in input, the energies will be calculated on a commensurate Monkhorst-Pack '
                       'mesh')

        for setting in self.mandatory_settings:
            assert getattr(self, setting, None), f'The mandatory key {setting} has not been provided'

        if self.name is None:
            self.name = 'ui'

        if self.directory is None:
            self.directory = '.'

        self._calculate()

    def _calculate(self):
        # The core of the calculation machinery is separated into self._calculate() to allow for monkeypatching
        # during testing

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

                self.write_results(directory='.')

                self.f_out.write(f'\tPrinting output in: {time() - reset:24.3f} sec\n')

                walltime = time() - start
                self.calc.results['walltime'] = walltime
                self.f_out.write(f'\n\tTotal time: {walltime:32.3f} sec\n')
                self.f_out.write('\nALL DONE\n\n')

                self.calc.results['job done'] = True

                # Unlink the output file
                delattr(self, 'f_out')

    def check_code_is_installed(self):
        # This calculator is entirely python-based, so we don't need to check for an installed binary
        return True

    def is_complete(self):
        return self.calc.results.get('job done', False)

    @property
    def do_smooth_interpolation(self):
        return any([f > 1 for f in self.smooth_int_factor])

    def get_k_point_weights(self):
        return np.ones(len(self.kpath.kpts))

    def get_number_of_spins(self):
        return 1

    def get_eigenvalues(self, kpt=None, spin=0):
        if spin != 0:
            raise NotImplementedError(
                f'Unfolding and interpolating calculator is not implemented for spin-polarised systems')

        if 'band structure' not in self.results:
            raise ValueError('You must first calculate the band structure before you try to access the KS eigenvalues')

        if kpt is None:
            return self.results['band structure'].energies[spin, :]
        else:
            return self.results['band structure'].energies[spin, kpt]

    def get_fermi_level(self):
        return 0

    @property
    def Emin(self):
        if self._Emin is None:
            return np.min(self.get_eigenvalues())
        else:
            return self._Emin

    @Emin.setter
    def Emin(self, value):
        self._Emin = value

    @property
    def Emax(self):
        if self._Emax is None:
            return np.max(self.get_eigenvalues())
        else:
            return self._Emax

    @Emax.setter
    def Emax(self, value):
        self._Emax = value

    @property
    def at(self):
        # basis vectors of direct lattice (in PC alat units)
        return self.calc.atoms.cell / self.alat

    @at.setter
    def at(self, value):
        self.calc.atoms.cell = value * self.alat

    @property
    def alat(self):
        return self.alat_sc / self.sc_dim[0]

    @property
    def bg(self):
        # basis vectors of reciprocal lattice (in PC 2pi/alat units)
        return np.linalg.inv(self.at).transpose()

    def todict(self):
        dct = super().todict()
        del dct['valid_settings']
        return dct
