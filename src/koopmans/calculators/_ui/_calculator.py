"""

The calculator class defining the Unfolding & interpolating calculator

"""

import copy
import json
import os
from datetime import datetime
from pathlib import Path
from time import time
from typing import List, Optional, Union

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.dft.dos import DOS
from ase.geometry.cell import crystal_structure_from_cell
from ase.spectrum.band_structure import BandStructure
from numpy.typing import ArrayLike, NDArray

from koopmans import utils
from koopmans.kpoints import Kpoints, kpath_to_dict
from koopmans.settings import (PlotSettingsDict,
                               UnfoldAndInterpolateSettingsDict)

from .._utils import CalculatorABC, CalculatorExt, sanitize_filenames
from ._atoms import UIAtoms
from ._utils import crys_to_cart, extract_hr, latt_vect


class UnfoldAndInterpolateCalculator(CalculatorExt, Calculator, CalculatorABC):
    # Subclass of CalculatorExt for performing unfolding and interpolation, using the base ASE calculator 'Calculator'

    ext_in = '.uii'
    ext_out = '.uio'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        self.parameters = UnfoldAndInterpolateSettingsDict()

        # Initialize first with the base ASE calculator, and then with the calculator extensions
        Calculator.__init__(self, atoms=atoms, *args, **kwargs)
        CalculatorExt.__init__(self, *args, **kwargs)

        # Ensure that self.atoms is a UIAtoms and not just a Atoms object
        if not isinstance(atoms, UIAtoms) and self.parameters.kgrid:
            atoms = UIAtoms.fromatoms(atoms=atoms, supercell_matrix=np.diag(self.parameters.kgrid))
        self.atoms = atoms
        self.atoms.calc = self

        # Intermediate variables
        self.centers: NDArray[np.float_] = np.array([])
        self.spreads: List[float] = []
        self.phases: List[complex] = []
        self.hr: NDArray[np.complex_] = np.array([])
        self.hr_coarse: NDArray[np.complex_] = np.array([])
        self.hr_smooth: NDArray[np.complex_] = np.array([])
        self.hk: NDArray[np.complex_] = np.array([])
        self.Rvec: NDArray[np.int_] = np.array([])
        self.Rsmooth: NDArray[np.int_] = np.array([])
        self.wRs: List[int] = []

        # Does not have a command (but we still want self.command to be defined)
        self.command = None

    @classmethod
    def fromfile(cls, filenames: Union[str, Path, List[str], List[Path]]) -> 'UnfoldAndInterpolateCalculator':
        sanitized_filenames = sanitize_filenames(filenames, cls.ext_in, cls.ext_out)

        calc = super(UnfoldAndInterpolateCalculator, cls).fromfile(sanitized_filenames)

        # If we were reading generating this object from files, look for bands, too
        if any([f.suffix == calc.ext_out for f in sanitized_filenames]):
            calc.read_bands()

        return calc

    def calculate(self):
        # Check mandatory settings
        for mandatory_setting in ['w90_seedname', 'kc_ham_file']:
            if mandatory_setting not in self.parameters:
                raise ValueError(f'You must provide the "{mandatory_setting}" setting for a UI calculation')

        # Check we have the requisite files
        if self.parameters.do_smooth_interpolation:
            if not self.parameters.dft_ham_file.is_file():
                raise FileNotFoundError(f'Cannot find {self.parameters.dft_ham_file} for smooth interpolation')
            if not self.parameters.dft_smooth_ham_file.is_file():
                raise FileNotFoundError(f'Cannot find {self.parameters.dft_smooth_ham_file} for smooth interpolation')

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

                self.f_out.write('UNFOLDING & INTERPOLATION\n\n')

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

                self.f_out.write(f'\tParsing input in:{time() - reset:27.3f} sec\n')
                reset = time()

                """
                 2) Core of the unfolding and interpolation code:
                    - build the map |i> ---> |Rn>
                    - calc interpolated bands (if needed)
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
                self.f_out.write('\nALL DONE\n')

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
                f'Unfolding and interpolating calculator is not implemented for spin-polarized systems')

        if 'band structure' not in self.results:
            raise ValueError('You must first calculate the band structure before you try to access the KS eigenvalues')

        if kpt is None:
            return self.results['band structure'].energies[spin, :]
        else:
            return self.results['band structure'].energies[spin, kpt]

    def get_fermi_level(self):
        return 0

    def parse_w90(self) -> None:
        '''
        centers : centers of WFs (in PC crystal units)
        spreads : spreads of WFs (in Ang^2)
        '''

        if len(self.centers) > 0 and len(self.spreads) > 0:
            num_wann = len(self.centers)
        else:
            with open(self.parameters.w90_seedname.with_suffix('.wout'), 'r') as ifile:
                lines = ifile.readlines()

            centers = []
            self.spreads = []
            count = 0

            for line in lines:
                if 'Number of Wannier Functions' in line:
                    num_wann = int(line.split()[6])
                if count > 0 and count <= num_wann:
                    start = line.find('(')
                    end = line.find(')')
                    centers.append(np.array(line[start + 1:end].replace(',', ' ').split(),
                                            dtype=float))
                    self.spreads.append(float(line.split()[-1]))
                    count += 1
                if 'Final State' in line:
                    count += 1

            self.centers = np.array(centers)

        self.Rvec = latt_vect(*self.parameters.kgrid)

        if self.parameters.w90_input_sc:
            self.parameters.num_wann_sc = num_wann
            self.parameters.num_wann = num_wann // np.prod(self.parameters.kgrid)
        else:
            self.parameters.num_wann = num_wann
            self.parameters.num_wann_sc = num_wann * np.prod(self.parameters.kgrid)

        self.centers /= np.linalg.norm(self.atoms.cell[0])
        self.centers = crys_to_cart(self.centers, self.atoms.acell.reciprocal(), -1)

        # generate the centers and spreads of all the other (R/=0) WFs
        self.centers = np.concatenate([self.centers + rvec for rvec in self.Rvec])
        self.spreads *= len(self.Rvec)

        return

    def parse_hr(self) -> None:
        """
        parse_hr reads the hamiltonian file passed as argument and it sets it up
        as self.hr. It also reads in the coarse and smooth Hamiltonians, if smooth
        interploation is being performed

        There is only one possible file format: the Wannier90 formatting. kcp files
        now also have this format. kc_occ_old and kc_emp_old have been deprecated.

        """

        # Read the Hamiltonian
        hr, rvect, _, nrpts = utils.read_wannier_hr_file(self.parameters.kc_ham_file)

        # Reshape the hamiltonian and convert it to a numpy array
        if nrpts == 1:
            assert len(hr) == self.parameters.num_wann_sc**2, \
                f'Wrong number of matrix elements ({len(hr)}) for the input hamiltonian'
            self.hr = np.array(hr, dtype=complex).reshape(self.parameters.num_wann_sc, self.parameters.num_wann_sc)
        else:
            assert len(hr) == nrpts * self.parameters.num_wann**2, \
                f'Wrong number of matrix elements ({len(hr)}) for the input hamiltonian'
            self.hr = np.array(hr, dtype=complex).reshape(nrpts, self.parameters.num_wann, self.parameters.num_wann)
            self.hr = extract_hr(self.hr, rvect, *self.parameters.kgrid)
            self.hr = self.hr.reshape(self.parameters.num_wann_sc, self.parameters.num_wann)

        # Reading the two Hamiltonians for the smooth interpolation method
        if self.parameters.do_smooth_interpolation:
            # The coarse Hamiltonian
            hr_coarse, rvect, _, nrpts = utils.read_wannier_hr_file(self.parameters.dft_ham_file)
            if nrpts == 1:
                assert len(hr_coarse) == self.parameters.num_wann_sc**2, \
                    f'Wrong number of matrix elements for hr_coarse {len(hr_coarse)}'
                self.hr_coarse = np.array(self.hr_coarse, dtype=complex)
                self.hr_coarse = self.hr_coarse.reshape(self.parameters.num_wann_sc, self.parameters.num_wann_sc)
                self.hr_coarse = self.hr_coarse[:, :self.parameters.num_wann]
            else:
                assert len(hr_coarse) == nrpts * \
                    self.parameters.num_wann**2, f'Wrong number of matrix elements for hr_coarse {len(hr_coarse)}'
                self.hr_coarse = np.array(hr_coarse, dtype=complex)
                self.hr_coarse = self.hr_coarse.reshape(nrpts, self.parameters.num_wann, self.parameters.num_wann)
                self.hr_coarse = extract_hr(self.hr_coarse, rvect, *self.parameters.kgrid)
                self.hr_coarse = self.hr_coarse.reshape(self.parameters.num_wann_sc, self.parameters.num_wann)

            # The smooth Hamiltonian
            hr_smooth, self.Rsmooth, self.wRs, nrpts = utils.read_wannier_hr_file(self.parameters.dft_smooth_ham_file)
            assert len(hr_smooth) == nrpts * \
                self.parameters.num_wann**2, f'Wrong number of matrix elements for hr_smooth {len(self.hr_smooth)}'
            self.hr_smooth = np.array(hr_smooth, dtype=complex)
            self.hr_smooth = self.hr_smooth.reshape(nrpts, self.parameters.num_wann, self.parameters.num_wann)

        return

    def parse_phases(self) -> None:
        """
        parse_phases gets the phases of WFs from the file 'wf_phases.dat'. If the file
                     is not found a warning is print out and the WFs phases are ignored.
        """

        try:
            with open('wf_phases.dat', 'r') as ifile:
                lines = ifile.readlines()
            self.phases = [float(l.split()[0]) + float(l.split()[1]) * 1j for l in lines]
        except FileNotFoundError:
            if self.parameters.w90_input_sc:
                utils.warn('file "wf_phases.dat" not found; phases are ignored')
            self.phases = []
        return

    def print_centers(self, centers: NDArray[np.float_] = np.array([])) -> None:
        """
        print_centers simply prints out the centers in the following Xcrysden-readable format:

                      X  0.000  0.000  0.000
                      X  0.000  0.000  0.000
                      *    *      *      *
                      *    *      *      *
                      *    *      *      *
                      X  0.000  0.000  0.000

        """

        if len(centers) == 0:
            centers = self.centers

        for n in range(self.parameters.num_wann_sc):
            self.f_out.write(' X' + ''.join([f'  {x:10.6f}' for x in centers[n]]) + '\n')

        return

    def write_results(self, directory: Optional[Path] = None) -> None:
        """
        write_results calls write_bands and write_dos if the DOS was calculated
        """
        if directory is None:
            directory = self.directory

        self.write_bands(directory)

        if self.parameters.do_dos:
            self.write_dos(directory)

        return

    def write_bands(self, directory=None) -> None:
        """
        write_bands prints the interpolated bands, in the QE format, in a file called
                    'bands_interpolated.dat'.
                    (see PP/src/bands.f90 around line 574 for the linearized path)
        """

        if directory is None:
            directory = self.directory

        kvec = []
        for kpt in self.parameters.kpath.kpts:
            kvec.append(crys_to_cart(kpt, self.atoms.acell.reciprocal(), +1))

        kx: List[float] = [0.0]
        for ik in range(1, len(kvec)):
            dxmod = float(np.linalg.norm(kvec[ik] - kvec[ik - 1]))
            if ik == 1:
                dxmod_save = dxmod
            if dxmod > 5 * dxmod_save:
                kx.append(kx[ik - 1])
            elif dxmod > 1.e-4:
                kx.append(kx[ik - 1] + dxmod)
                dxmod_save = dxmod
            else:
                kx.append(kx[ik - 1] + dxmod)

        bs = self.results['band structure'].energies
        for energies_spin, label in zip(bs, ['up', 'down']):
            fname = 'bands_interpolated'
            if bs.shape[0] == 2:
                fname += f'_spin_{label}'

            with open(f'{directory}/{fname}.dat', 'w') as ofile:
                ofile.write('# Written at ' + datetime.now().isoformat(timespec='seconds'))

                for energies in energies_spin.transpose():
                    assert len(kx) == len(energies)
                    for k, energy in zip(kx, energies):
                        ofile.write(f'\n{k:16.8f}{energy:16.8f}')
                    ofile.write('\n')

        return

    def read_bands(self, directory: Optional[Path] = None) -> None:
        """
        read_bands reads the interpolated bands, in the QE format, in a file called
                   'bands_interpolated.dat'
                   (see PP/src/bands.f90 around line 574 for the linearized path)

                   This function also then regenerates the DOS based off the bandstructure
        """

        if directory is None:
            directory = self.directory

        energies: List[List[List[float]]] = []
        for suffix in ['', '_spin_up', '_spin_down']:
            band_file = directory / f'bands_interpolated{suffix}.dat'
            if os.path.isfile(band_file):
                energies.append([[]])
                with open(band_file, 'r') as f:
                    flines = f.readlines()
                for line in flines[1:]:
                    splitline = line.strip().split()
                    if len(splitline) == 0:
                        energies[-1].append([])
                    else:
                        energies[-1][-1].append(float(splitline[-1]))

        if len(energies) > 0:
            self.results['band structure'] = BandStructure(
                path=self.parameters.kpath, energies=np.transpose(energies, (0, 2, 1)))

        self.calc_dos()

    def write_dos(self, directory=None) -> None:
        """
        write_dos prints the DOS in a file called 'dos_interpolated.dat', in a format (E , DOS(E))

        """
        if directory is None:
            directory = self.directory

        with open(f'{directory}/dos_interpolated.dat', 'w') as ofile:
            ofile.write('# Written at ' + datetime.now().isoformat(timespec='seconds'))
            dos = self.results['dos']
            for e, d in zip(dos.get_energies(), dos.get_dos()):
                ofile.write('\n{:10.4f}{:12.6f}'.format(e, d))
            ofile.write('\n')

        return

    def write_input(self, atoms: Atoms) -> None:
        """
        write_input writes out a JSON file containing the settings used for the calculation. This "input" file is
        never actually used in a standard calculation, but it is useful for debugging
        """

        with utils.chdir(self.directory):
            with open(f'{self.prefix}{self.ext_in}', 'w') as fd:
                settings = copy.deepcopy(self.parameters.data)

                # Remove the kpoints information from the settings dict
                kgrid = settings.pop('kgrid')
                kpath = settings.pop('kpath')

                # Remove the plot parameters from the settings dict
                plotting = settings.pop('plotting')

                # Converting Paths to JSON-serialisable strings
                for k in self.parameters.are_paths:
                    if k in settings:
                        settings[k] = str(settings[k])

                # Store all the settings in one big dictionary
                bigdct = {"workflow": {"task": "ui"}, "ui": settings}

                # Provide the bandpath information in the form of a string
                bigdct['kpoints'] = {'grid': kgrid, **kpath_to_dict(kpath)}
                # The cell is stored elsewhere
                bigdct['kpoints'].pop('cell')

                # Provide the plot information
                bigdct['plotting'] = {k: v for k, v in plotting.items()}

                # We also need to provide a cell so the explicit kpath can be reconstructed from the string alone
                bigdct['atoms'] = {'cell_parameters': utils.construct_cell_parameters_block(atoms)}

                json.dump(bigdct, fd, indent=2)

    def read_input(self, input_file: Optional[Path] = None) -> None:
        """
        read_input reads in the settings from a JSON-formatted input file and loads them onto this calculator (useful
        for restarting)
        """

        if input_file is None:
            input_file = self.directory / (self.prefix + self.ext_in)

        with open(input_file, 'r') as fd:

            # Load the input file
            bigdct = json.load(fd)

        assert bigdct['workflow']['task'] == 'ui', 'UI input file should have "task": "ui"'

        # Update the parameters
        self.parameters = bigdct['ui']

        # Update plot parameters
        self.parameters.plotting = PlotSettingsDict(**bigdct.get('plotting', {}))

        # Load the cell and kpts if they are provided
        if 'atoms' in bigdct:
            utils.read_cell_parameters(self.atoms, bigdct['atoms'].get('cell_parameters', {}))

        kpoint_block = bigdct.get('kpoints', {})
        if kpoint_block:
            kpts = Kpoints(**kpoint_block, cell=self.atoms.cell)
            self.parameters.kgrid = kpts.grid
            self.parameters.kpath = kpts.path

        return

    def read_results(self) -> None:
        """
        read_results parses a .uio file for the solitary result 'job done'
        """

        output_file = self.directory / (self.prefix + self.ext_out)

        assert output_file.is_file()

        # Check the calculation is done
        with open(output_file, 'r') as f:
            flines = f.readlines()
        self.results = {'job done': any(['ALL DONE' in line for line in flines])}

        return

    def interpolate(self, start_time) -> None:
        """
        interpolate is the main program in this module and it calls consecutively
                    the three independent functions:
                    - map_wannier
                    - calc_bands
                    - calc_dos

        """

        # Step 1: map the WFs
        if self.parameters.do_map:
            self.map_wannier()
            self.f_out.write(f'\tBuilding the map |i> --> |Rn> in:{time()-start_time:11.3f} sec\n')
        reset = time()

        # Step 2: calculate the electronic bands along kpath
        self.calc_bands()
        self.f_out.write(f'\tCalculating bands in: {time()-reset:22.3f} sec\n')
        reset = time()

        # Step 3: calculate the density-of-states
        if self.parameters.do_dos:
            self.calc_dos()
            self.f_out.write(f'\tCalculating DOS in: {time()-reset:24.3f} sec\n')

        return

    def map_wannier(self) -> None:
        """
        map_wannier builds the map |i> --> |Rn> between the WFs in the SC and in the PC.
        """

        centers = []
        spreads = []
        index = []

        # here we identify the WFs within the R=0 cell
        self.centers /= self.parameters.kgrid
        self.centers -= np.floor(self.centers)
        self.centers *= self.parameters.kgrid
        for n in range(self.parameters.num_wann_sc):
            if all([x - 1 < 1.e-3 for x in self.centers[n]]):
                centers.append(self.centers[n])
                spreads.append(self.spreads[n])
                index.append(n)

        # check on the WFs found in the R=0 cell
        assert len(centers) == self.parameters.num_wann, 'Did not find the right number of WFs in the R=0 cell'

        # here we identify with |Rn> the WFs in the rest of the SC, by comparing centers and spreads
        # the WFs are now ordered as (R0,1),(R0,2),...,(R0,n),(R1,1),...
        for rvect in self.Rvec[1:]:
            count = 0
            for m in range(self.parameters.num_wann):
                for n in range(self.parameters.num_wann_sc):
                    if all(abs(self.centers[n] - centers[m] - rvect) < 1.e-3) and \
                       abs(self.spreads[n] - spreads[m]) < 1.e-3:
                        centers.append(self.centers[n])
                        spreads.append(self.spreads[n])
                        index.append(n)
                        count += 1
            assert count == self.parameters.num_wann, f'Found {count} WFs in the {rvect} cell'

        # permute phases and Hamiltonian matrix elements in order to follow the new order of WFs
        hr = [self.hr[i, j] for i in index for j in index]
        if self.phases:
            self.phases = [self.phases[i] for i in index]

        self.centers = np.array(centers, dtype=float)
        self.spreads = spreads
        self.hr = np.array(hr, dtype=complex).reshape(self.parameters.num_wann_sc, self.parameters.num_wann_sc)

        return

    def calc_bands(self) -> None:
        """
        calc_bands interpolates the k-space hamiltonian along the input path, by Fourier
                   transforming the Wannier hamiltonian H(R). The function generates two
                   new attributes:
                   - self.hk containing H(k) for any k-vector in the input path
                   - self.results['band structure'] containing the interpolated electronic energies

        """

        # when smooth interpolation is on, we remove the DFT part from hr
        hr = self.hr[:, :self.parameters.num_wann]
        if self.parameters.do_smooth_interpolation:
            hr = hr - self.hr_coarse
        hr = hr.reshape(len(self.Rvec), self.parameters.num_wann, self.parameters.num_wann)

        # renormalize H(R) on the WF phases
        if self.phases:
            hr = np.conjugate(self.phases) * (hr.transpose() * self.phases).transpose()

        # calculate phase and phase correction
        # phi:      (Nkpath, NR)
        # phi_corr: (Nkpath, NR, num_wann, num_wann)
        phi = np.exp(2j * np.pi * np.dot(self.parameters.kpath.kpts, self.Rvec.transpose()))
        phi_corr = self.correct_phase()

        # interpolate H(k)
        hk = np.transpose(np.sum(phi * np.transpose(hr * phi_corr, axes=(2, 3, 0, 1)), axis=3), axes=(2, 0, 1))
        if self.parameters.do_smooth_interpolation:
            phi = np.exp(2j * np.pi * np.dot(self.parameters.kpath.kpts, self.Rsmooth.transpose()))
            hr_smooth = np.transpose(self.hr_smooth, axes=(2, 1, 0)) / self.wRs
            hk += np.dot(phi, np.transpose(hr_smooth, axes=(1, 2, 0)))

        bands = np.linalg.eigvalsh(hk)
        self.hk = hk
        self.results['band structure'] = BandStructure(self.parameters.kpath, [bands])

        return

    def calc_dos(self) -> None:
        """
        calc_dos calculates the density of states using the DOS function from ASE
        """

        self.results['dos'] = DOS(self, width=self.parameters.plotting.degauss, window=(
            self.parameters.plotting.Emin, self.parameters.plotting.Emax),
            npts=self.parameters.plotting.nstep + 1)

        return

    def correct_phase(self) -> NDArray[np.complex_]:
        """
        correct_phase calculate the correct phase factor to put in the Fourier transform
                      to get the interpolated k-space hamiltonian. The correction consists
                      of finding the right distance, i.e. the right R-vector, considering
                      also the BVK boundary conditions.
                      if use_ws_distance=True, the function accounts also for the intracell
                      distance between Wannier functions, otherwise only the intercell
                      distances are considered.

           IMPORTANT: the vectors must all be in crystal units otherwise the distances are
                      not properly evaluated.
        """

        if self.parameters.use_ws_distance:
            # create an array containing all the distances between reference (R=0) WFs and all the other WFs:
            # 1) accounting for their positions within the unit cell
            wf_dist = np.concatenate([self.centers] * self.parameters.num_wann) \
                - np.concatenate([[c] * self.parameters.num_wann_sc for c in self.centers[:self.parameters.num_wann]])

        else:
            # 2) considering only the distance between the unit cells they belong to
            wf_dist = np.array(np.concatenate([[rvec] * self.parameters.num_wann for rvec in self.Rvec]).tolist()
                               * self.parameters.num_wann)

        # supercell lattice vectors
        Tvec = [np.array((i, j, k)) * self.parameters.kgrid for i in range(-1, 2)
                for j in range(-1, 2) for k in range(-1, 2)]
        Tlist = []
        for dist in wf_dist:
            distance = crys_to_cart(dist + np.array(Tvec), self.atoms.acell, +1)
            norm = np.linalg.norm(distance, axis=1)
            Tlist.append(np.where(norm - norm.min() < 1.e-3)[0])

        phase = np.zeros((len(self.parameters.kpath.kpts), len(Tlist)), dtype=complex)
        for i, t_index in enumerate(Tlist):
            for ik, kvect in enumerate(self.parameters.kpath.kpts):
                for it in t_index:
                    phase[ik, i] += np.exp(2j * np.pi * np.dot(kvect, Tvec[it]))
                phase[ik, i] /= len(t_index)

        phase = phase.reshape(len(self.parameters.kpath.kpts), self.parameters.num_wann, len(self.Rvec),
                              self.parameters.num_wann)
        phase = np.transpose(phase, axes=(0, 2, 3, 1))

        return phase
