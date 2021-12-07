"""

The calculator class defining the Unfolding & interpolating calculator

"""

import os
import copy
import json
from time import time
import numpy as np
from numpy.typing import ArrayLike
from typing import Union, List, Optional
from pathlib import Path
from datetime import datetime
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.dft.dos import DOS
from ase.spectrum.band_structure import BandStructure
from koopmans import utils
from koopmans.settings import UnfoldAndInterpolateSettingsDict
from .._utils import CalculatorExt, CalculatorABC, sanitise_filenames
from ._utils import crys_to_cart, extract_hr
from ._atoms import UIAtoms


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
        if not isinstance(atoms, UIAtoms) and self.parameters.kgrid:
            atoms = UIAtoms.fromatoms(atoms=atoms, supercell_matrix=np.diag(self.parameters.kgrid))
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

    @classmethod
    def fromfile(cls, filenames: Union[str, Path, List[str], List[Path]]) -> 'UnfoldAndInterpolateCalculator':
        sanitised_filenames = sanitise_filenames(filenames, cls.ext_in, cls.ext_out)

        calc = super(UnfoldAndInterpolateCalculator, cls).fromfile(sanitised_filenames)

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

    def parse_w90(self):
        '''
        centers : centers of WFs (in PC crystal units)
        spreads : spreads of WFs (in Ang^2)
        '''

        if self.centers and self.spreads:
            num_wann = len(self.centers)
        else:
            with open(self.parameters.w90_seedname.with_suffix('.wout'), 'r') as ifile:
                lines = ifile.readlines()

            self.centers = []
            self.spreads = []
            count = 0

            for line in lines:
                if 'Number of Wannier Functions' in line:
                    num_wann = int(line.split()[6])
                if count > 0 and count <= num_wann:
                    start = line.find('(')
                    end = line.find(')')
                    self.centers.append(np.array(line[start + 1:end].replace(',', ' ').split(),
                                                 dtype=float))
                    self.spreads.append(float(line.split()[-1]))
                    count += 1
                if 'Final State' in line:
                    count += 1

        self.Rvec = np.array([[x, y, z] for x in range(self.parameters.kgrid[0])
                              for y in range(self.parameters.kgrid[1]) for z in range(self.parameters.kgrid[2])])

        if self.parameters.w90_input_sc:
            self.parameters.num_wann_sc = num_wann
            self.parameters.num_wann = num_wann // np.prod(self.parameters.kgrid)

        else:
            self.parameters.num_wann = num_wann
            self.parameters.num_wann_sc = num_wann * np.prod(self.parameters.kgrid)

        for n in range(num_wann):
            self.centers[n] = self.centers[n] / np.linalg.norm(self.atoms.cell[0])
            self.centers[n] = crys_to_cart(self.centers[n], self.atoms.acell.reciprocal(), -1)

        # generate the centers and spreads of all the other (R/=0) WFs
        if not self.parameters.w90_input_sc:
            for rvect in self.Rvec[1:]:
                for n in range(self.parameters.num_wann):
                    self.centers.append(self.centers[n] + rvect)
                    self.spreads.append(self.spreads[n])

        return

    def parse_hr(self):
        """
        parse_hr reads the hamiltonian file passed as argument and it sets it up
        as self.hr. It also reads in the coarse and smooth Hamiltonians, if smooth
        interploation is being performed

        There is only one possible file format: the Wannier90 formatting. kcp files
        now also have this format. kc_occ_old and kc_emp_old have been deprecated.

        """

        # Read the Hamiltonian
        hr, rvect, _, nrpts = utils.read_hr_file(self.parameters.kc_ham_file)

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
            hr_coarse, rvect, _, nrpts = utils.read_hr_file(self.parameters.dft_ham_file)
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
            hr_smooth, self.Rsmooth, self.wRs, nrpts = utils.read_hr_file(self.parameters.dft_smooth_ham_file)
            assert len(hr_smooth) == nrpts * \
                self.parameters.num_wann**2, f'Wrong number of matrix elements for hr_smooth {len(self.hr_smooth)}'
            self.hr_smooth = np.array(hr_smooth, dtype=complex)
            self.hr_smooth = self.hr_smooth.reshape(nrpts, self.parameters.num_wann, self.parameters.num_wann)

        return

    def parse_phases(self):
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
            self.phases = [1 for _ in range(self.parameters.num_wann_sc)]
        return

    def print_centers(self, centers=None):
        """
        print_centers simply prints out the centers in the following Xcrysden-readable format:

                      X  0.000  0.000  0.000
                      X  0.000  0.000  0.000
                      *    *      *      *
                      *    *      *      *
                      *    *      *      *
                      X  0.000  0.000  0.000

        """

        if centers is None:
            centers = self.centers

        for n in range(self.parameters.num_wann_sc):
            self.f_out.write(' X' + ''.join([f'  {x:10.6f}' for x in centers[n]]) + '\n')

        return

    def write_results(self, directory: Optional[Path] = None):
        """
        write_results calls write_bands and write_dos if the DOS was calculated
        """
        if directory is None:
            directory = self.directory

        self.write_bands(directory)

        if self.parameters.do_dos:
            self.write_dos(directory)

        return

    def write_bands(self, directory=None):
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

        kx = [0.]
        for ik in range(1, len(kvec)):
            dxmod = np.linalg.norm(kvec[ik] - kvec[ik - 1])
            if ik == 1:
                dxmod_save = dxmod

            if dxmod > 5 * dxmod_save:
                kx.append(kx[ik - 1])

            elif dxmod > 1.e-4:
                kx.append(kx[ik - 1] + dxmod)
                dxmod_save = dxmod

            else:
                kx.append(kx[ik - 1] + dxmod)

        with open(f'{directory}/bands_interpolated.dat', 'w') as ofile:
            ofile.write('# Written at ' + datetime.now().isoformat(timespec='seconds'))

            for energies in self.results['band structure'].energies[0].transpose():
                assert len(kx) == len(energies)
                for k, energy in zip(kx, energies):
                    ofile.write(f'\n{k:16.8f}{energy:16.8f}')
                ofile.write('\n')

        return

    def read_bands(self, directory=None):
        """
        read_bands reads the interpolated bands, in the QE format, in a file called
                   'bands_interpolated.dat'
                   (see PP/src/bands.f90 around line 574 for the linearized path)

                   This function also then regenerates the DOS based off the bandstructure
        """

        if directory is None:
            directory = self.directory

        band_file = f'{directory}/bands_interpolated.dat'
        bands = [[]]
        if os.path.isfile(band_file):
            with open(band_file, 'r') as f:
                flines = f.readlines()
            for line in flines[1:]:
                splitline = line.strip().split()
                if len(splitline) == 0:
                    bands.append([])
                else:
                    bands[-1].append(float(splitline[-1]))
            self.results['band structure'] = BandStructure(path=self.parameters.kpath, energies=[np.transpose(bands)])
            self.calc_dos()

    def write_dos(self, directory=None):
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

    def write_input(self, atoms: Atoms):
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

                # Converting Paths to JSON-serialisable strings
                for k in self.parameters.are_paths:
                    if k in settings:
                        settings[k] = str(settings[k])

                # Store all the settings in one big dictionary
                bigdct = {"workflow": {"task": "ui"}, "ui": settings}

                # Provide the bandpath information in the form of a string
                bigdct['setup'] = {'kpoints': {'kpath': kpath.path, 'kgrid': kgrid}}

                # We also need to provide a cell so the explicit kpath can be reconstructed from the string alone
                bigdct['setup']['cell_parameters'] = utils.construct_cell_parameters_block(atoms)

                json.dump(bigdct, fd, indent=2)

    def read_input(self, input_file: Optional[Path] = None):
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

        # Load the cell and kpts if they are provided
        if 'setup' in bigdct:
            cell = utils.read_cell_parameters(None, bigdct['setup'].get('cell_parameters', {}))
            if cell:
                self.atoms.cell = cell
            kpoint_block = bigdct['setup'].get('kpoints', {})
            if kpoint_block:
                self.parameters.kgrid = kpoint_block['kgrid']
                utils.read_kpath(self, kpoint_block['kpath'])

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

    def interpolate(self, start_time):
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
            self.f_out.write(f'\tBuilding the map |i> --> |Rn> in:\t{time()-start_time:.3f} sec\n')
        reset = time()

        # Step 2: calculate the electronic bands along kpath
        self.calc_bands()
        self.f_out.write(f'\tCalculating bands in: {time()-reset:22.3f} sec\n')
        reset = time()

        # Step 3 (optional) : calculate the density-of-states
        if self.parameters.do_dos:
            self.calc_dos()
            self.f_out.write(f'\tCalculating DOS in: {time()-reset:24.3f} sec\n')

        return

    def map_wannier(self):
        """
        map_wannier builds the map |i> --> |Rn> between the WFs in the SC and in the PC.
        """

        centers = []
        spreads = []
        index = []
        num_wann = 0

        # here we identify the WFs within the R=0 cell
        for n in range(self.parameters.num_wann_sc):
            # shift the WFs within the SC
            self.centers[n] = self.centers[n] / self.parameters.kgrid
            self.centers[n] = self.centers[n] - np.floor(self.centers[n])

            # converting centers from crystal units of SC to crystal units of PC
            self.centers[n] = self.centers[n] * self.parameters.kgrid

            if all([x - 1 < 1.e-3 for x in self.centers[n]]):
                centers.append(self.centers[n])
                spreads.append(self.spreads[n])
                index.append(n)
                num_wann += 1

        # check on the WFs found in the R=0 cell
        assert num_wann == self.parameters.num_wann, 'Did not find the right number of WFs in the R=0 cell'

        # here we identify with |Rn> the WFs in the rest of the SC
        # the WFs are now ordered as (R0,1),(R0,2),...,(R0,n),(R1,1),...
        for rvect in self.Rvec[1:]:
            count = 0
            for m in range(self.parameters.num_wann):
                for n in range(self.parameters.num_wann_sc):
                    wf_dist = self.centers[n] - centers[m]
                    if (np.linalg.norm(wf_dist - rvect) < 1.e-3) and (abs(self.spreads[n] - spreads[m] < 1.e-3)):
                        centers.append(self.centers[n])
                        spreads.append(self.spreads[n])
                        index.append(n)
                        count += 1

            assert count == self.parameters.num_wann, f'Found {count} WFs in the {rvect} cell'

        # redefine phases and hr in order to follow the new order of WFs
        phases = []
        hr = np.zeros((self.parameters.num_wann_sc, self.parameters.num_wann_sc), dtype=complex)
        for n in range(self.parameters.num_wann_sc):
            phases.append(self.phases[index[n]])
            for m in range(self.parameters.num_wann_sc):
                hr[m, n] = self.hr[index[m], index[n]]

        self.centers = centers
        self.spreads = spreads
        self.hr = hr
        self.phases = phases

        return

    def calc_bands(self):
        """
        calc_bands interpolates the k-space hamiltonian along the input path, by Fourier
                   transforming the Wannier hamiltonian H(R). The function generates two
                   new attributes:
                   - self.hk containing H(k) for any k-vector in the input path
                   - self.results['band structure'] containing the interpolated electronic energies

        """

        # when smooth interpolation is on, we remove the DFT part from hr
        hr = np.array(self.hr[:, :self.parameters.num_wann], dtype=complex)
        if self.parameters.do_smooth_interpolation:
            hr = hr - self.hr_coarse

        # renormalize H(R) on the WF phases
        for m in range(self.parameters.num_wann_sc):
            for n in range(self.parameters.num_wann):
                hr[m, n] = self.phases[m].conjugate() * hr[m, n] * self.phases[n]

        # here we build the interpolated H(k)
        hk = np.zeros((len(self.parameters.kpath.kpts), self.parameters.num_wann,
                      self.parameters.num_wann), dtype=complex)
        bands = []
        self.f_out.write(f"\n\t\tTotal number of k-points: {len(self.parameters.kpath.kpts):6d}\n\n")
        for ik, kpt in enumerate(self.parameters.kpath.kpts):
            self.f_out.write(f"\t\t      calculating point # {ik+1}\n")
            for m in range(self.parameters.num_wann):
                for n in range(self.parameters.num_wann):
                    for ir in range(len(self.Rvec)):

                        mm = ir * self.parameters.num_wann + m
                        phase = self.correct_phase(self.centers[n], self.centers[mm], self.Rvec[ir], kpt)

                        hk[ik, m, n] = hk[ik, m, n] + np.exp(1j * 2 * np.pi * np.dot(kpt, self.Rvec[ir])) * \
                            phase * hr[mm, n]

                    if self.parameters.do_smooth_interpolation:
                        hk_smooth = 0
                        for jr in range(len(self.Rsmooth)):
                            Rs = self.Rsmooth[jr]
                            hk_smooth = hk_smooth + np.exp(1j * 2 * np.pi * np.dot(kpt, Rs)) * \
                                self.hr_smooth[jr, m, n] / self.wRs[jr]
                        hk[ik, m, n] = hk[ik, m, n] + hk_smooth

            bands.append(np.linalg.eigvalsh(hk[ik]))

        self.f_out.write('\n')

        self.hk = hk
        self.results['band structure'] = BandStructure(self.parameters.kpath, [bands])

        return

    def calc_dos(self) -> None:
        """
        calc_dos calculates the density of states using a gaussian smearing. The DOS is saved
                 as a list [ [E1, DOS(E1)], [E2, DOS[E2]], ... , [En, DOS(En)] ]
        """

        if self.parameters.Emin is None:
            self.parameters.Emin = np.min(self.get_eigenvalues() - 0.1)
        if self.parameters.Emax is None:
            self.parameters.Emax = np.max(self.get_eigenvalues() + 0.1)

        self.results['dos'] = DOS(self, width=self.parameters.degauss, window=(
            self.parameters.Emin, self.parameters.Emax), npts=self.parameters.nstep + 1)

        return

    def correct_phase(self, center_ref, center, rvect, kvect):
        """
        correct_phase calculate the correct phase factor to put in the Fourier transform
                      to get the interpolated k-space hamiltonian. The correction consists
                      of finding the right distance, i.e. the right R-vector, considering
                      also the BVK boundary conditions.
                      if use_ws_distance=True, the function accounts also for the intracell
                      distance between Wannier functions, otherwise only the intercell
                      distances are considered.

           IMPORTANT: the input vectors (center_ref, center, Rvec and kvec) must all
                      be in crystal units otherwise the distances are not properly evaluated.
        """

        if self.parameters.use_ws_distance:
            wf_dist = crys_to_cart(center - center_ref, self.atoms.acell, +1)
        else:
            wf_dist = crys_to_cart(rvect, self.atoms.acell, +1)

        dist_min = np.linalg.norm(wf_dist)
        Tlist = []

        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    tvect = np.array([i, j, k]) * self.parameters.kgrid
                    Tvec = crys_to_cart(tvect, self.atoms.acell, +1)
                    dist = np.linalg.norm(wf_dist + Tvec)

                    if (abs(dist - dist_min) < 1.e-3):
                        #
                        # an equidistant replica is found
                        #
                        Tlist.append(tvect)

                    elif (dist < dist_min):
                        #
                        # a new nearest replica is found
                        # reset dist_min and reinitialize Tlist
                        #
                        dist_min = dist
                        Tlist = [tvect]

                    else:
                        #
                        # this replica is rejected
                        #
                        continue

        phase = sum(np.exp(1j * 2 * np.pi * np.dot(Tlist, kvect))) / len(Tlist)

        return phase
