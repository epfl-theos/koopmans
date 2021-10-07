"""

The calculator class defining the Unfolding & interpolating calculator

"""

from time import time
from ase.calculators.calculator import Calculator
from ase.units import Hartree
import numpy as np
from numpy.typing import ArrayLike
from typing import Union, List, Optional
from pathlib import Path
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.dft.kpoints import BandPath
from .._utils import CalculatorExt, CalculatorABC, sanitise_filenames
from ._atoms import UIAtoms
from ._io import parse_w90, parse_hr, parse_phases, print_centers, write_results, write_bands, write_dos, \
    write_input, read_input, read_results, read_bands
from ._interpolate import interpolate, calc_bands, correct_phase, calc_dos
from koopmans import utils
from koopmans.settings import UnfoldAndInterpolateSettingsDict


class UnfoldAndInterpolateCalculator(CalculatorExt, Calculator, CalculatorABC):
    # Subclass of CalculatorExt for performing unfolding and interpolation, using the base ASE calculator 'Calculator'

    ext_in = '.uii'
    ext_out = '.uio'
    results_for_qc = ['band structure', 'dos']

    def __init__(self, atoms: Atoms, *args, **kwargs):
        self.parameters = UnfoldAndInterpolateSettingsDict()

        # Initialise first with the base ASE calculator, and then with the calculator extensions
        Calculator.__init__(self, atoms=atoms, *args, **kwargs)
        CalculatorExt.__init__(self, *args, **kwargs)

        # Ensure that self.atoms is a UIAtoms and not just a Atoms object
        if not isinstance(atoms, UIAtoms):
            atoms = UIAtoms.fromatoms(atoms=atoms, supercell_matrix=np.diag(self.parameters.kpts))
        self.atoms = atoms
        self.atoms.calc = self

        # Intermediate variables
        self.centers: ArrayLike = []
        self.spreads: ArrayLike = []
        self.phases: ArrayLike = []
        self.hr: ArrayLike = []
        self.hr_smooth: ArrayLike = []
        self.hk: ArrayLike = []
        self.Rvec: ArrayLike = []
        self.Rsmooth: ArrayLike = []
        self.wRs: ArrayLike = []

    # Adding all the various functions defined in ._io
    parse_w90 = parse_w90
    parse_hr = parse_hr
    parse_phases = parse_phases
    print_centers = print_centers
    write_results = write_results
    write_bands = write_bands
    write_dos = write_dos
    write_input = write_input
    read_input = read_input
    read_results = read_results
    read_bands = read_bands
    # Adding various functions defined in ._interpolate
    interpolate = interpolate
    calc_bands = calc_bands
    calc_dos = calc_dos
    correct_phase = correct_phase

    @classmethod
    def fromfile(cls, filenames: Union[str, Path, List[str], List[Path]]) -> 'UnfoldAndInterpolateCalculator':
        calc = CalculatorABC.fromfile(filenames)

        sanitised_filenames = sanitise_filenames(filenames, cls.ext_in, cls.ext_out)

        # If we were reading generating this object from files, look for bands, too
        if any([f.suffix == calc.ext_out for f in sanitised_filenames]):
            calc.read_bands()

        return calc

    def calculate(self):
        # Check mandatory settings
        for mandatory_setting in ['w90_seedname', 'kc_ham_file']:
            if mandatory_setting not in self.parameters:
                raise ValueError(f'You must provide the "{mandatory_setting}" setting for a UI calculation')

        # Check we have the requisite files
        if self.parameters.do_smooth_interpolation:
            assert self.parameters.dft_ham_file.is_file(), 'Missing file_hr_coarse for smooth interpolation'
            assert self.parameters.dft_smooth_ham_file.is_file(), 'Missing dft_smooth_ham_file for smooth interpolation'

        if self.prefix is None:
            self.prefix = 'ui'

        if self.directory is None:
            self.directory = '.'

        self._calculate()

    def _calculate(self):
        # The core of the calculation machinery is separated into self._calculate() to allow for monkeypatching
        # during testing

        self.write_input(self.atoms)

        start = time()
        reset = time()

        with utils.chdir(self.directory):
            with open(f'{self.prefix}{self.ext_out}', 'w') as f_out:
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

                self.results['job done'] = True

                # Unlink the output file
                delattr(self, 'f_out')

    def check_code_is_installed(self):
        # This calculator is entirely python-based, so we don't need to check for an installed binary
        return True

    def is_complete(self):
        return self.results.get('job done', False)

    def is_converged(self):
        return True

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
