"""

The process class defining the Unfolding & interpolating step

"""

import copy
import json
import os
from datetime import datetime
from pathlib import Path
from time import time
from typing import List, Optional, Tuple, Type, Union

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.dft.dos import DOS
from ase.geometry.cell import crystal_structure_from_cell
from ase.spectrum.band_structure import BandStructure
from numpy.typing import ArrayLike, NDArray
from pydantic import BaseModel, ConfigDict, Field

from koopmans import calculators, utils
from koopmans.kpoints import Kpoints, kpath_to_dict
from koopmans.settings import (PlotSettingsDict,
                               UnfoldAndInterpolateSettingsDict)

from .._process import Process
from ._atoms import UIAtoms
from ._utils import crys_to_cart, extract_hr, latt_vect


class UnfoldAndInterpolateInputs(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)

    atoms: UIAtoms = \
        Field(..., description='An ASE Atoms object augmented with supercell information')

    parameters: UnfoldAndInterpolateSettingsDict = \
        Field(..., description='The parameters for the unfolding and interpolation process')

    centers: NDArray[np.float64] = \
        Field(..., description="centers of the Wannier functions (in primitive cell crystal units)")

    spreads: List[float] = \
        Field(..., description="spreads of Wannier functions (in Ang^2)")

    kc_ham_file: Tuple[utils.HasDirectoryAttr, Path] = \
        Field(..., description='The Koopmans Hamiltonian')

    dft_ham_file: Tuple[utils.HasDirectoryAttr, Path] = \
        Field(..., description='The DFT Hamiltonian')

    dft_smooth_ham_file: Optional[Tuple[utils.HasDirectoryAttr, Path]] = \
        Field(default=None, description='The DFT Hamiltonian, evaluated on the smooth grid')

    wf_phases_file: Optional[Tuple[utils.HasDirectoryAttr, Path]] = \
        Field(default=None, description='The file containing the phases of the Wannier functions')

    plotting_parameters: PlotSettingsDict = \
        Field(default=PlotSettingsDict(), description='The plotting parameters (used when generating the DOS)')


class UnfoldAndInterpolateOutputs(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)

    band_structure: BandStructure
    dos: DOS | None = None


class UnfoldAndInterpolateProcess(Process):
    __slots__ = ['inputs', 'outputs', '_centers', '_spreads', '_Rvec', '_hr',
                 '_hr_smooth', '_hr_coarse', '_phases', '_wRs', '_Rsmooth', '_hk']

    def __init__(self, atoms: Atoms, centers, parameters, **kwargs):

        # Update the parameters field based on the centers and spreads objects
        num_wann = len(centers)
        if parameters.w90_input_sc:
            parameters.num_wann_sc = num_wann
            parameters.num_wann = num_wann // np.prod(parameters.kgrid)
        else:
            parameters.num_wann = num_wann
            parameters.num_wann_sc = num_wann * np.prod(parameters.kgrid)

        # Convert atoms from an Atoms to a UIAtoms object
        ui_atoms = UIAtoms.fromatoms(atoms, supercell_matrix=np.diag(parameters.kgrid))

        super().__init__(atoms=ui_atoms, centers=centers, parameters=parameters, **kwargs)

    @property
    def _input_model(self) -> Type[UnfoldAndInterpolateInputs]:
        return UnfoldAndInterpolateInputs

    @property
    def _output_model(self) -> Type[UnfoldAndInterpolateOutputs]:
        return UnfoldAndInterpolateOutputs

    def _run(self):
        self._centers = self.inputs.centers
        self._spreads = self.inputs.spreads
        self._Rvec = None
        self._hr = None
        self._hr_smooth = None
        self._hr_coarse = None
        self._phases = None
        self._wRs = None
        self._Rsmooth = None

        self.write_input(self.inputs.atoms)

        """
         1) Parse data:
            - calc parameters from the JSON file
            - other parameters from W90 output file
            - hamiltonian(s)
            - WFs phases read from file wf_phases.dat
        """

        # Generate the centers and spreads the non-primitive-cell (R/=0) Wannier functions
        self._Rvec = latt_vect(*self.inputs.parameters.kgrid)
        self._centers /= np.linalg.norm(self.inputs.atoms.cell[0])
        self._centers = crys_to_cart(self._centers, self.inputs.atoms.acell.reciprocal(), -1)

        self._centers = np.concatenate([self._centers + rvec for rvec in self._Rvec])
        self._spreads *= len(self._Rvec)

        self.parse_hr()
        self.parse_phases()

        """
         2) Core of the unfolding and interpolation code:
            - build the map |i> ---> |Rn>
            - calc interpolated bands (if needed)
            - calc DOS (if needed)
        """

        self.interpolate()

        """
         3) Print out the results:
            - bands into 'bands_interpolated.dat' file
            - DOS into 'dos_interpolated.dat' file
        """

        self.write_results()

    def parse_w90(self) -> None:
        '''
        centers : centers of WFs (in PC crystal units)
        spreads : spreads of WFs (in Ang^2)
        '''

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
        hr, rvect, _, nrpts = utils.read_wannier_hr_file(self.inputs.kc_ham_file)

        # Reshape the hamiltonian and convert it to a numpy array
        if nrpts == 1:
            assert len(hr) == self.inputs.parameters.num_wann_sc**2, \
                f'Wrong number of matrix elements ({len(hr)}) for the input hamiltonian'
            self._hr = np.array(hr, dtype=complex).reshape(
                self.inputs.parameters.num_wann_sc, self.inputs.parameters.num_wann_sc)
        else:
            assert len(hr) == nrpts * self.inputs.parameters.num_wann**2, \
                f'Wrong number of matrix elements ({len(hr)}) for the input hamiltonian'
            self._hr = np.array(hr, dtype=complex).reshape(
                nrpts, self.inputs.parameters.num_wann, self.inputs.parameters.num_wann)
            self._hr = extract_hr(self._hr, rvect, *self.inputs.parameters.kgrid)
            self._hr = self._hr.reshape(self.inputs.parameters.num_wann_sc, self.inputs.parameters.num_wann)

        # Reading the two Hamiltonians for the smooth interpolation method
        if self.inputs.parameters.do_smooth_interpolation:
            # The coarse Hamiltonian
            hr_coarse, rvect, _, nrpts = utils.read_wannier_hr_file(self.inputs.dft_ham_file)
            if nrpts == 1:
                assert len(hr_coarse) == self.inputs.parameters.num_wann_sc**2, \
                    f'Wrong number of matrix elements for hr_coarse {len(hr_coarse)}'
                self._hr_coarse = np.array(self._hr_coarse, dtype=complex)
                self._hr_coarse = self._hr_coarse.reshape(
                    self.inputs.parameters.num_wann_sc, self.inputs.parameters.num_wann_sc)
                self._hr_coarse = self._hr_coarse[:, :self.inputs.parameters.num_wann]
            else:
                assert len(hr_coarse) == nrpts * \
                    self.inputs.parameters.num_wann**2, f'Wrong number of matrix elements for hr_coarse {len(hr_coarse)}'
                self._hr_coarse = np.array(hr_coarse, dtype=complex)
                self._hr_coarse = self._hr_coarse.reshape(
                    nrpts, self.inputs.parameters.num_wann, self.inputs.parameters.num_wann)
                self._hr_coarse = extract_hr(self._hr_coarse, rvect, *self.inputs.parameters.kgrid)
                self._hr_coarse = self._hr_coarse.reshape(
                    self.inputs.parameters.num_wann_sc, self.inputs.parameters.num_wann)

            # The smooth Hamiltonian
            hr_smooth, self._Rsmooth, self._wRs, nrpts = utils.read_wannier_hr_file(self.inputs.dft_smooth_ham_file)
            assert len(hr_smooth) == nrpts * \
                self.inputs.parameters.num_wann**2, f'Wrong number of matrix elements for hr_smooth {len(self._hr_smooth)}'
            self._hr_smooth = np.array(hr_smooth, dtype=complex)
            self._hr_smooth = self._hr_smooth.reshape(
                nrpts, self.inputs.parameters.num_wann, self.inputs.parameters.num_wann)

        return

    def parse_phases(self) -> None:
        """
        parse_phases gets the phases of WFs from the file 'wf_phases.dat'. If the file
                     is not found a warning is print out and the WFs phases are ignored.
        """

        try:
            with open('wf_phases.dat', 'r') as ifile:
                lines = ifile.readlines()
            self._phases = [float(l.split()[0]) + float(l.split()[1]) * 1j for l in lines]
        except FileNotFoundError:
            if self.inputs.parameters.w90_input_sc:
                utils.warn('file "wf_phases.dat" not found; phases are ignored')
            self._phases = []
        return

    def write_results(self, directory: Optional[Path] = None) -> None:
        """
        write_results calls write_bands and write_dos if the DOS was calculated
        """
        if directory is None:
            directory = self.directory

        self.write_bands(directory)

        if self.inputs.parameters.do_dos:
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
        for kpt in self.inputs.parameters.kpath.kpts:
            kvec.append(crys_to_cart(kpt, self.inputs.atoms.acell.reciprocal(), +1))

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

        bs = self.outputs.band_structure.energies
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

    def write_dos(self, directory=None) -> None:
        """
        write_dos prints the DOS in a file called 'dos_interpolated.dat', in a format (E , DOS(E))

        """
        if directory is None:
            directory = self.directory

        with open(f'{directory}/dos_interpolated.dat', 'w') as ofile:
            ofile.write('# Written at ' + datetime.now().isoformat(timespec='seconds'))
            dos = self.outputs.dos
            for e, d in zip(dos.get_energies(), dos.get_dos()):
                ofile.write('\n{:10.4f}{:12.6f}'.format(e, d))
            ofile.write('\n')

        return

    def write_input(self, atoms: Atoms) -> None:
        """
        write_input writes out a JSON file containing the settings used for the calculation. This "input" file is
        never actually used in a standard calculation, but it is useful for debugging
        """

        with open(f'{self.name}.json', 'w') as fd:
            settings = copy.deepcopy(self.inputs.parameters.data)

            # Remove the kpoints information from the settings dict
            kgrid = settings.pop('kgrid')
            kpath = settings.pop('kpath')

            # Converting Paths to JSON-serialisable strings
            for k in self.inputs.parameters.are_paths:
                if k in settings:
                    settings[k] = str(settings[k])

            # Store all the settings in one big dictionary
            bigdct = {"workflow": {"task": "ui"}, "ui": settings}

            # Provide the bandpath information in the form of a string
            bigdct['kpoints'] = {'grid': kgrid, **kpath_to_dict(kpath)}
            # The cell is stored elsewhere
            bigdct['kpoints'].pop('cell')

            # Provide the plot information
            bigdct['plotting'] = {k: v for k, v in self.inputs.plotting_parameters.data.items()}

            # We also need to provide a cell so the explicit kpath can be reconstructed from the string alone
            bigdct['atoms'] = {'cell_parameters': utils.construct_cell_parameters_block(atoms)}

            json.dump(bigdct, fd, indent=2)

    def interpolate(self):
        """
        interpolate is the main program in this module and it calls consecutively
                    the three independent functions:
                    - map_wannier
                    - calc_bands
                    - generate_dos

        """

        # Step 1: map the WFs
        if self.inputs.parameters.do_map:
            self.map_wannier()

        # Step 2: calculate the electronic bands along kpath
        bs = self.calc_bands()

        # Step 3: calculate the density-of-states
        dos = generate_dos(bs, self.inputs.plotting_parameters) if self.inputs.parameters.do_dos else None

        # Store the outputs
        self.outputs = self._output_model(band_structure=bs, dos=dos)

        return

    def map_wannier(self) -> None:
        """
        map_wannier builds the map |i> --> |Rn> between the WFs in the SC and in the PC.
        """

        centers = []
        spreads = []
        index = []

        # here we identify the WFs within the R=0 cell
        self._centers /= self.inputs.parameters.kgrid
        self._centers -= np.floor(self._centers)
        self._centers *= self.inputs.parameters.kgrid
        for n in range(self.inputs.parameters.num_wann_sc):
            if all([x - 1 < 1.e-3 for x in self._centers[n]]):
                centers.append(self._centers[n])
                spreads.append(self._spreads[n])
                index.append(n)

        # check on the WFs found in the R=0 cell
        assert len(centers) == self.inputs.parameters.num_wann, 'Did not find the right number of WFs in the R=0 cell'

        # here we identify with |Rn> the WFs in the rest of the SC, by comparing centers and spreads
        # the WFs are now ordered as (R0,1),(R0,2),...,(R0,n),(R1,1),...
        for rvect in self._Rvec[1:]:
            count = 0
            for m in range(self.inputs.parameters.num_wann):
                for n in range(self.inputs.parameters.num_wann_sc):
                    if all(abs(self._centers[n] - centers[m] - rvect) < 1.e-3) and \
                       abs(self._spreads[n] - spreads[m]) < 1.e-3:
                        centers.append(self._centers[n])
                        spreads.append(self._spreads[n])
                        index.append(n)
                        count += 1
            assert count == self.inputs.parameters.num_wann, f'Found {count} WFs in the {rvect} cell'

        # permute phases and Hamiltonian matrix elements in order to follow the new order of WFs
        hr = [self._hr[i, j] for i in index for j in index]
        if self._phases:
            self._phases = [self._phases[i] for i in index]

        self._centers = np.array(centers, dtype=float)
        self._spreads = spreads
        self._hr = np.array(hr, dtype=complex).reshape(
            self.inputs.parameters.num_wann_sc, self.inputs.parameters.num_wann_sc)

        return

    def calc_bands(self) -> BandStructure:
        """
        calc_bands interpolates the k-space hamiltonian along the input path, by Fourier
                   transforming the Wannier hamiltonian H(R). The function generates two
                   new attributes:
                   - self.hk containing H(k) for any k-vector in the input path
                   - self.results['band structure'] containing the interpolated electronic energies

        """

        # when smooth interpolation is on, we remove the DFT part from hr
        hr = self._hr[:, :self.inputs.parameters.num_wann]
        if self.inputs.parameters.do_smooth_interpolation:
            hr = hr - self._hr_coarse
        hr = hr.reshape(len(self._Rvec), self.inputs.parameters.num_wann, self.inputs.parameters.num_wann)

        # renormalize H(R) on the WF phases
        if self._phases:
            hr = np.conjugate(self._phases) * (hr.transpose() * self._phases).transpose()

        # calculate phase and phase correction
        # phi:      (Nkpath, NR)
        # phi_corr: (Nkpath, NR, num_wann, num_wann)
        phi = np.exp(2j * np.pi * np.dot(self.inputs.parameters.kpath.kpts, self._Rvec.transpose()))
        phi_corr = self.correct_phase()

        # interpolate H(k)
        hk = np.transpose(np.sum(phi * np.transpose(hr * phi_corr, axes=(2, 3, 0, 1)), axis=3), axes=(2, 0, 1))
        if self.inputs.parameters.do_smooth_interpolation:
            phi = np.exp(2j * np.pi * np.dot(self.inputs.parameters.kpath.kpts, self._Rsmooth.transpose()))
            hr_smooth = np.transpose(self._hr_smooth, axes=(2, 1, 0)) / self._wRs
            hk += np.dot(phi, np.transpose(hr_smooth, axes=(1, 2, 0)))

        bands = np.linalg.eigvalsh(hk)
        self._hk = hk
        return BandStructure(self.inputs.parameters.kpath, [bands])

    def correct_phase(self) -> NDArray[np.complex128]:
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

        if self.inputs.parameters.use_ws_distance:
            # create an array containing all the distances between reference (R=0) WFs and all the other WFs:
            # 1) accounting for their positions within the unit cell
            wf_dist = np.concatenate([self._centers] * self.inputs.parameters.num_wann) \
                - np.concatenate([[c] * self.inputs.parameters.num_wann_sc for c in self._centers[:self.inputs.parameters.num_wann]])

        else:
            # 2) considering only the distance between the unit cells they belong to
            wf_dist = np.array(np.concatenate([[rvec] * self.inputs.parameters.num_wann for rvec in self._Rvec]).tolist()
                               * self.inputs.parameters.num_wann)

        # supercell lattice vectors
        Tvec = [np.array((i, j, k)) * self.inputs.parameters.kgrid for i in range(-1, 2)
                for j in range(-1, 2) for k in range(-1, 2)]
        Tlist = []
        for dist in wf_dist:
            distance = crys_to_cart(dist + np.array(Tvec), self.inputs.atoms.acell, +1)
            norm = np.linalg.norm(distance, axis=1)
            Tlist.append(np.where(norm - norm.min() < 1.e-3)[0])

        phase = np.zeros((len(self.inputs.parameters.kpath.kpts), len(Tlist)), dtype=complex)
        for i, t_index in enumerate(Tlist):
            for ik, kvect in enumerate(self.inputs.parameters.kpath.kpts):
                for it in t_index:
                    phase[ik, i] += np.exp(2j * np.pi * np.dot(kvect, Tvec[it]))
                phase[ik, i] /= len(t_index)

        phase = phase.reshape(len(self.inputs.parameters.kpath.kpts), self.inputs.parameters.num_wann, len(self._Rvec),
                              self.inputs.parameters.num_wann)
        phase = np.transpose(phase, axes=(0, 2, 3, 1))

        return phase


def generate_dos(band_structure: BandStructure, plotting_parameters: PlotSettingsDict) -> DOS:
    """
    Generate the density of states using the DOS function from ASE
    """

    dos = DOS(None,
              width=plotting_parameters.degauss,
              window=(plotting_parameters.Emin, plotting_parameters.Emax),
              npts=plotting_parameters.nstep + 1,
              w_k=np.ones(len(band_structure.path.kpts)),
              nspins=1,
              e_skn=np.array([band_structure.energies[0, :]]))

    return dos
