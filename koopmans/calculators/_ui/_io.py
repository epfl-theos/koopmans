'''
I/O functions associated with the UI calculator

'''

import os
import numpy as np
import copy
import json
from datetime import datetime
from ase.atoms import Atoms
from ase.calculators.calculator import FileIOCalculator
from ase.spectrum.band_structure import BandStructure
from koopmans import utils
from ._utils import latt_vect, crys_to_cart, extract_hr


def parse_w90(self):
    '''
    centers : centers of WFs (in PC crystal units)
    spreads : spreads of WFs (in Ang^2)
    '''

    with open(self.w90_seedname + '.wout', 'r') as ifile:
        lines = ifile.readlines()

    at = np.empty((3, 3), dtype=float)
    self.centers = []
    self.spreads = []
    count = 0

    for line in lines:
        if 'Number of Wannier Functions' in line:
            num_wann = int(line.split()[6])
        if ' a_1 ' in line:
            at[0, :] = line.split()[1:]
        if ' a_2 ' in line:
            at[1, :] = line.split()[1:]
        if ' a_3 ' in line:
            at[2, :] = line.split()[1:]
        if count > 0 and count <= num_wann:
            start = line.find('(')
            end = line.find(')')
            self.centers.append(np.array(line[start + 1:end].replace(',', ' ').split(),
                                         dtype=float))
            self.spreads.append(float(line.split()[-1]))
            count += 1
        if 'Final State' in line:
            count += 1

    self.Rvec = latt_vect(*self.sc_dim)
    self.at = at / self.alat

    if self.w90_input_sc:
        self.num_wann_sc = num_wann
        self.num_wann = num_wann // np.prod(self.sc_dim)

        # converting at and bg from the SC to the PC
        self.at_sc = copy.copy(self.at)
        self.at = np.array([x / n for x, n in zip(self.at, self.sc_dim)])
        self.bg_sc = copy.copy(self.bg)
        self.bg = np.array([x * n for x, n in zip(self.bg, self.sc_dim)])

    else:
        self.num_wann = num_wann
        self.num_wann_sc = num_wann * np.prod(self.sc_dim)
        self.at_sc = np.array([x * n for x, n in zip(self.at, self.sc_dim)])
        self.bg_sc = np.array([x / n for x, n in zip(self.bg, self.sc_dim)])

    for n in range(num_wann):
        self.centers[n] = self.centers[n] / self.alat
        self.centers[n] = crys_to_cart(self.centers[n], self.bg, -1)

    # generate the centers and spreads of all the other (R/=0) WFs
    if not self.w90_input_sc:
        for rvect in self.Rvec[1:]:
            for n in range(self.num_wann):
                self.centers.append(self.centers[n] + rvect)
                self.spreads.append(self.spreads[n])

    return


def parse_hr(self):
    """
    parse_hr reads the hamiltonian file passed as argument and it sets it up
             as attribute self.hr

    there are 3 possible types of file:
      - w90 for a Wannier90 type of file (now also CP-koopmans files have this format)
      - kc_occ_old for the old type of CP hamiltonian for occ states (hamiltonian1.xml)
      - kc_emp_old for the old type of CP hamiltonian for emp states (hamiltonian_emp.dat)

    NB: kc_emp_old must be called 'hamiltonian_emp.dat' otherwise the code may crash
        or misread the matrix elements. if the file name is different the code
        should be updated.

    """

    with open(self.kc_ham_file, 'r') as ifile:
        lines = ifile.readlines()

    if 'written on' in lines[0].lower():
        hr_type = 'w90'
    elif 'xml version' in lines[0]:
        hr_type = 'kc_occ_old'
    elif self.kc_ham_file[-19:] == 'hamiltonian_emp.dat':
        hr_type = 'kc_emp_old'
    else:
        raise ValueError('Hamiltonian file not recognised')

    self.hr = []

    if hr_type == 'w90':

        nrpts = int(lines[2].split()[0])
        if (nrpts == 1):
            single_R = True
        else:
            single_R = False

        if single_R:
            assert int(lines[1].split()[0]) == self.num_wann_sc, 'In parse_hr inconsistency in num_wann'
        else:
            assert int(lines[1].split()[0]) != self.num_wann, 'In parse_hr inconsistency in num_wann'

        if not single_R:
            rvect = []
        lines_to_skip = 3 + int(nrpts / 15)
        if nrpts % 15 > 0:
            lines_to_skip += 1
        counter = 0

        for line in lines[lines_to_skip:]:
            self.hr.append(float(line.split()[5]) + 1j * float(line.split()[6]))
            counter += 1
            if not single_R and counter == self.num_wann**2:
                rvect.append(np.array(line.split()[0:3], dtype=int))
                counter = 0

    if hr_type == 'kc_occ_old':
        for line in lines[5:-2]:
            self.hr.append(line.split()[0])

    if hr_type == 'kc_emp_old':
        for line in lines:
            self.hr.append(line.split()[0])

    if hr_type == 'w90' and not single_R:
        assert len(self.hr) == nrpts * self.num_wann**2, \
            f'Wrong number of matrix elements ({len(self.hr)}) for the input hamiltonian'
        self.hr = np.array(self.hr, dtype=complex).reshape(nrpts, self.num_wann, self.num_wann)
        self.hr = extract_hr(self.hr, rvect, self.sc_dim[0], self.sc_dim[1], self.sc_dim[2])
        self.hr = self.hr.reshape(self.num_wann_sc, self.num_wann)
    else:
        assert len(self.hr) == self.num_wann_sc**2, \
            f'Wrong number of matrix elements ({len(self.hr)}) for the input hamiltonian'
        self.hr = np.array(self.hr, dtype=complex).reshape(self.num_wann_sc, self.num_wann_sc)

    # conversion to eV (hamiltonian from CP Koopmans code is in Hartree)
    if hr_type == 'kc_occ_old' or hr_type == 'kc_emp_old':
        self.hr = self.hr * utils.units.Hartree

    # check the hermiticity of the hamiltonian (except for H_R(m,n))
    if not (hr_type == 'w90' and not single_R):
        for m in range(self.num_wann):
            for n in range(self.num_wann):
                assert self.hr[m, n] - self.hr[n, m].conjugate() <= 1.e-6, \
                    f'Hamiltonian matrix (index {m}, {n}) not hermitian'

    # reading the 2 hamiltonians for the smooth interpolation method
    if self.do_smooth_interpolation:
        self.hr_coarse = []
        self.hr_smooth = []

        # parsing hr_coarse
        with open(self.dft_ham_file, 'r') as ifile:
            lines = ifile.readlines()

        nrpts = int(lines[2].split()[0])
        if nrpts == 1:
            single_R = True
        else:
            single_R = False

        if single_R:
            assert int(lines[1].split()[0]) == self.num_wann_sc, 'In parse_hr inconsistency in num_wann in hr_coarse'
        else:
            assert int(lines[1].split()[0]) == self.num_wann, 'In parse_hr inconsistency in num_wann in hr_coarse'

        if not single_R:
            rvect = []
        lines_to_skip = 3 + int(nrpts / 15)
        if nrpts % 15 > 0:
            lines_to_skip += 1
        counter = 0

        for line in lines[lines_to_skip:]:
            self.hr_coarse.append(float(line.split()[5]) + 1j * float(line.split()[6]))
            counter += 1
            if not single_R and counter == self.num_wann**2:
                rvect.append(np.array(line.split()[0:3], dtype=int))
                counter = 0

        if single_R:
            assert len(self.hr_coarse) == self.num_wann_sc**2, \
                f'Wrong number of matrix elements for hr_coarse {len(self.hr_coarse)}'
            self.hr_coarse = np.array(self.hr_coarse, dtype=complex)
            self.hr_coarse = self.hr_coarse.reshape(self.num_wann_sc, self.num_wann_sc)
            self.hr_coarse = self.hr_coarse[:, :self.num_wann]
        else:
            assert len(self.hr_coarse) == nrpts * \
                self.num_wann**2, f'Wrong number of matrix elements for hr_coarse {len(self.hr_coarse)}'
            self.hr_coarse = np.array(self.hr_coarse, dtype=complex)
            self.hr_coarse = self.hr_coarse.reshape(nrpts, self.num_wann, self.num_wann)
            self.hr_coarse = extract_hr(self.hr_coarse, rvect, self.sc_dim[0], self.sc_dim[1], self.sc_dim[2])
            self.hr_coarse = self.hr_coarse.reshape(self.num_wann_sc, self.num_wann)

        # parsing hr_smooth
        with open(self.dft_smooth_ham_file, 'r') as ifile:
            lines = ifile.readlines()

        assert int(lines[1].split()[0]) == self.num_wann, 'In parse_hr inconsistency in num_wann in hr_smooth'

        weights = []
        rvect = []
        nrpts = int(lines[2].split()[0])
        lines_to_skip = 3 + int(nrpts / 15)
        if nrpts % 15 > 0:
            lines_to_skip += 1
        counter = 0

        for line in lines[3:lines_to_skip]:
            for n in range(len(line.split())):
                weights.append(int(line.split()[n]))

        for line in lines[lines_to_skip:]:
            self.hr_smooth.append(float(line.split()[5]) + 1j * float(line.split()[6]))
            counter += 1
            if counter == self.num_wann**2:
                rvect.append(np.array(line.split()[0:3], dtype=int))
                counter = 0

        assert len(self.hr_smooth) == nrpts * \
            self.num_wann**2, f'Wrong number of matrix elements for hr_smooth {len(self.hr_smooth)}'
        self.wRs = weights
        self.Rsmooth = rvect
        self.hr_smooth = np.array(self.hr_smooth, dtype=complex)
        self.hr_smooth = self.hr_smooth.reshape(nrpts, self.num_wann, self.num_wann)

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
        if self.w90_input_sc:
            utils.warn('file "wf_phases.dat" not found; phases are ignored')
        self.phases = [1 for _ in range(self.num_wann_sc)]
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

    for n in range(self.num_wann_sc):
        self.f_out.write(' X' + ''.join([f'  {x:10.6f}' for x in centers[n]]) + '\n')

    return


def write_results(self, directory=None):
    """
    write_results calls write_bands and write_dos if the DOS was calculated
    """
    if directory is None:
        directory = self.directory

    self.write_bands(directory)

    if self.do_dos:
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
    for kpt in self.kpath.kpts:
        kvec.append(crys_to_cart(kpt, self.bg, +1))

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
        self.results['band structure'] = BandStructure(path=self.kpath, energies=[np.transpose(bands)])
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


def write_input_file(self):
    """
    write_input_file writes out a JSON file containing the settings used for the calculation. This "input" file is
    never actually used in a standard calculation, but it is useful for debugging
    """

    with utils.chdir(self.directory):
        with open(f'{self.name}{self.ext_in}', 'w') as fd:
            settings = copy.deepcopy(self.calc.parameters)

            # Convert alat_sc to Bohr
            settings['alat_sc'] *= utils.units.Bohr

            # Remove the ASE BandPath object
            kpath = settings.pop('kpath')

            # Store all the settings in one big dictionary
            bigdct = {"workflow": {"task": "ui"}, "ui": settings}

            # Provide the bandpath information in the form of a string
            bigdct['setup'] = {'kpoints': {'kpath': kpath.path}}

            # We also need to provide a cell so the explicit kpath can be reconstructed from the string alone
            bigdct['setup']['cell_parameters'] = utils.construct_cell_parameters_block(self.calc)

            json.dump(bigdct, fd, indent=2)


def read_input_file(self, input_file=None):
    """
    read_input_file reads in the settings from a JSON-formatted input file and returns an ASE calculator (useful for
    restarting)
    """

    if input_file is None:
        input_file = self.directory + '/' + self.name + self.ext_in

    with open(input_file, 'r') as fd:

        # Load the input file
        bigdct = json.load(fd)

        assert bigdct['workflow']['task'] == 'ui', 'UI input file should have "task": "ui"'

        # Create an empty dummy ASE calculator required by read_ui_dict
        atoms = Atoms(calculator=FileIOCalculator())
        atoms.calc.atoms = atoms

        # Load the cell if it is provided
        if 'setup' in bigdct:
            cell = utils.read_cell_parameters(None, bigdct['setup'].get('cell_parameters', {}))
            if cell:
                atoms.cell = cell
            kpath = bigdct['setup'].get('kpoints', {}).get('kpath', None)
            if kpath:
                utils.read_kpath(atoms.calc, kpath)

        # Parse the UI dict
        return read_ui_dict(bigdct['ui'], atoms)


def read_output_file(self, output_file=None):
    """
    read_output_file reads in an output file and returns an ASE calculator object with the solitary result 'job done'
    """

    if output_file is None:
        output_file = f'{self.directory}/{self.name}{self.ext_out}'

    # Create an empty dummy ASE calculator
    atoms = Atoms(calculator=FileIOCalculator())
    atoms.calc.atoms = atoms
    calc = atoms.calc
    calc.directory = output_file.rsplit('/', 1)[0]

    assert os.path.isfile(output_file)

    # Check the calculation is done
    with open(output_file, 'r') as f:
        flines = f.readlines()
    calc.results = {'job done': any(['ALL DONE' in line for line in flines])}

    return calc

def read_ui_dict(dct, generic_atoms):
    # Use a generic calculator with no command
    atoms = copy.deepcopy(generic_atoms)
    calc = atoms.calc
    calc.atoms.calc = calc

    calc.command = ''

    # Overwrite the parameters with the provided JSON dict
    calc.parameters = utils.parse_dict(dct)

    # Use kpath specified in the setup block if it is not present in the ui dict
    setup_kpath = generic_atoms.calc.parameters.get('kpath', None)
    if 'kpath' not in dct and setup_kpath:
        calc.parameters['kpath'] = setup_kpath

    # Convert units of alat_sc
    if 'alat_sc' in calc.parameters:
        calc.parameters['alat_sc'] *= utils.units.Bohr

    return calc