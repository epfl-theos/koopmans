"""

The calculator class defining the Unfolding & interpolating calculator

"""

from typing import List
from time import time
import numpy as np
from ase.dft.kpoints import BandPath
from .._utils import ExtendedCalculator
from ._settings import valid_settings, UISettingsDictWithChecks
from koopmans import utils


class UnfoldAndInterpolateCalculator(ExtendedCalculator):
    # Subclass of ExtendedCalculator for performing unfolding and interpolation
    from ._io import parse_w90, parse_hr, parse_phases, print_centers, write_results, write_bands, write_dos, \
        write_input_file, read_input_file, read_output_file, read_bands
    from ._interpolate import interpolate, calc_bands, correct_phase, calc_dos

    def __init__(self, calc=None, qe_files=[], skip_qc=False, **kwargs):

        self.parameters = UISettingsDictWithChecks(settings=valid_settings,
                                                   are_paths=['w90_seedname', 'kc_ham_file',
                                                              'dft_ham_file', 'dft_smooth_ham_file'],
                                                   to_not_parse=[],
                                                   physicals=['alat_sc', 'degauss', 'Emin', 'Emax'])
        # Link to the corresponding ASE Calculator (it does not use ASE)
        self._ase_calc_class = None

        # Define the appropriate file extensions
        self.ext_in = '.uii'
        self.ext_out = '.uio'

        super().__init__(calc, qe_files, skip_qc, **kwargs)

        # If we were reading generating this object from files, look for bands, too
        if qe_files and any([self.ext_out in f for f in qe_files]):
            self.read_bands()

        self.results_for_qc = ['band structure', 'dos']

        if self.parameters.kpath is None:
            # By default, use ASE's default bandpath for this cell (see
            # https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#brillouin-zone-data)
            if any(self.atoms.pbc):
                self.parameters.kpath = self.atoms.cell.get_bravais_lattice().bandpath()
            else:
                self.parameters.kpath = 'G'
        if isinstance(self.parameters.kpath, str):
            # If kpath is provided as a string, convert it to a BandPath first
            utils.read_kpath(self, self.parameters.kpath)
        assert isinstance(self.parameters.kpath, BandPath)

        if self.parameters.do_smooth_interpolation:
            assert self.parameters.dft_ham_file, 'Missing file_hr_coarse for smooth interpolation'
            assert self.parameters.dft_smooth_ham_file, 'Missing dft_smooth_ham_file for smooth interpolation'

    def calculate(self):
        # Check mandatory settings
        for mandatory_setting in ['w90_seedname', 'kc_ham_file', 'alat_sc', 'sc_dim']:
            if mandatory_setting not in self.parameters:
                raise ValueError(f'You must provide the "{mandatory_setting}" setting for a UI calculation')

        if self.parameters.kpath is None:
            utils.warn('"kpath" missing in input, the energies will be calculated on a commensurate Monkhorst-Pack '
                       'mesh')

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
                self.results['walltime'] = walltime
                self.f_out.write(f'\n\tTotal time: {walltime:32.3f} sec\n')
                self.f_out.write('\nALL DONE\n\n')

                self.calc.results['job done'] = True

                # Unlink the output file
                delattr(self, 'f_out')

    def check_code_is_installed(self):
        # This calculator is entirely python-based, so we don't need to check for an installed binary
        return True

    def is_complete(self):
        return self.results.get('job done', False)

    def get_k_point_weights(self):
        return np.ones(len(self.parameters.kpath.kpts))

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

    def todict(self):
        dct = super().todict()
        del dct['valid_settings']
        return dct
