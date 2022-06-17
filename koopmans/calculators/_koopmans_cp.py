"""

kcp calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

from __future__ import annotations
import os
import copy
import math
import numpy as np
import pickle
from pathlib import Path
from scipy.linalg import block_diag
from typing import Optional, List, Union
from pandas.core.series import Series
import xml.etree.ElementTree as ET
from ase import Atoms
from ase.calculators.espresso import Espresso_kcp
from koopmans import utils, settings, pseudopotentials, bands
from koopmans.commands import ParallelCommand
from ._utils import CalculatorExt, CalculatorABC, bin_directory, CalculatorCanEnforceSpinSym


def read_ham_file(filename: Path) -> np.ndarray:
    # Read a single hamiltonian XML file
    if not filename.exists():
        raise FileExistsError(f'{filename} does not exist')

    with open(filename, 'r') as fd:
        tree = ET.parse(fd)
    ham_xml = tree.getroot()

    length = int(math.sqrt(int(ham_xml.attrib['size'])))

    assert ham_xml.text is not None, f'{filename} is empty'

    ham_array = np.array([complex(*[float(x) for x in line.split(',')])
                          for line in ham_xml.text.strip().split('\n')], dtype=complex) * utils.units.Hartree

    return ham_array.reshape((length, length))


class KoopmansCPCalculator(CalculatorCanEnforceSpinSym, CalculatorExt, Espresso_kcp, CalculatorABC):
    # Subclass of CalculatorExt for performing calculations with kcp.x
    ext_in = '.cpi'
    ext_out = '.cpo'

    def __init__(self, atoms: Atoms, skip_qc: bool = False, alphas: Optional[List[List[float]]] = None,
                 filling: Optional[List[List[bool]]] = None, **kwargs):
        # Define the valid parameters
        self.parameters = settings.KoopmansCPSettingsDict()


        # Initialise first using the ASE parent and then CalculatorExt
        Espresso_kcp.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, skip_qc, **kwargs)

        # Add nelec, nelup, neldw if they are missing
        if 'nelec' not in self.parameters and 'pseudopotentials' in self.parameters:
            self.parameters.nelec = pseudopotentials.nelec_from_pseudos(
                self.atoms, self.pseudopotentials, self.parameters.pseudo_dir)
        if 'nelec' in self.parameters:
            if 'nelup' not in self.parameters:
                self.parameters.nelup = self.parameters.nelec // 2
            if 'neldw' not in self.parameters:
                self.parameters.neldw = self.parameters.nelec // 2

        if not isinstance(self.command, ParallelCommand):
            self.command = ParallelCommand(os.environ.get(
                'ASE_ESPRESSO_KCP_COMMAND', str(bin_directory) + os.path.sep + self.command))

        if alphas is not None:
            self.alphas = alphas
        if filling is not None:
            self.filling = filling
        self.results_for_qc = ['energy', 'homo_energy', 'lumo_energy']

        # Give the calculator an attribute to keep track of which band has been held fixed, for calculations where
        # fixed_state = .true.. N.B. this differs from self.parameters.fixed_band in the case of empty orbitals (see
        # koopmans.workflows._koopmans_dscf.py for more details)
        self.fixed_band: Optional[bands.Band] = None


    def calculate(self):
        # kcp.x imposes nelup >= neldw, so if we try to run a calcualtion with neldw > nelup, swap the spin channels
        if self.parameters.nspin == 2:
            spin_channels_are_swapped = self.parameters.nelup < self.parameters.neldw
        else:
            spin_channels_are_swapped = False

        # Swap the spin channels
        if spin_channels_are_swapped:
            self._swap_spin_channels()

        # Write out screening parameters to file
        if self.parameters.get('do_orbdep', False):
            self.write_alphas()

        super().calculate()

        # Check spin-up and spin-down eigenvalues match
        if 'eigenvalues' in self.results and self.parameters.do_outerloop \
                and self.parameters.nspin == 2 and self.parameters.tot_magnetization == 0 \
                and not self.parameters.fixed_state and len(self.results['eigenvalues']) > 0:
            rms_eigenval_difference = np.sqrt(np.mean(np.diff(self.results['eigenvalues'], axis=0)**2))
            if rms_eigenval_difference > 0.05:
                utils.warn('Spin-up and spin-down eigenvalues differ substantially')

        # Swap the spin channels back
        if spin_channels_are_swapped:
            self._swap_spin_channels()

    def _swap_spin_channels(self):
        # Parameters
        if self.parameters.fixed_band is not None and self.parameters.fixed_state:
            if self.parameters.nbnd is None:
                if self.parameters.nspin == 2:
                    nbup = self.parameters.nelup
                    nbdw = self.parameters.neldw
                else:
                    nbup = self.parameters.nelec // 2
                    nbdw = self.parameters.nelec // 2
            else:
                nbup = self.parameters.nbnd
                nbdw = self.parameters.nbnd
            if self.parameters.fixed_band > nbup:
                # The fixed band was spin-down
                self.parameters.fixed_band -= nbup
            else:
                # The fixed band was spin-up
                self.parameters.fixed_band += nbdw
        self.parameters.nelup, self.parameters.neldw = self.parameters.neldw, self.parameters.nelup
        self.parameters.tot_magnetization *= -1

        # alphas and filling
        self.alphas = self.alphas[::-1]
        self.filling = self.filling[::-1]

        # Results
        if 'orbital_data' in self.results:
            self.results['orbital_data'] = {k: v[::-1] for k, v in self.results['orbital_data'].items()}
        for key in ['eigenvalues', 'lambda']:
            if key in self.results:
                self.results[key] = self.results[key][::-1]

        # Input and output files
        for nd in [self.parameters.ndr, self.parameters.ndw]:
            outdir = self.parameters.outdir / f'{self.parameters.prefix}_{nd}.save/K00001'

            for fpath_1 in outdir.glob('*1.*'):
                # Swap the two files around
                fpath_tmp = fpath_1.parent / fpath_1.name.replace('1', 'tmp')
                fpath_2 = fpath_1.parent / fpath_1.name.replace('1', '2')

                if not fpath_2.exists():
                    raise FileNotFoundError(
                        'Error in {self.__class__.__name__}._swap_spin_channels: I expected {fpath_2} to exist')

                fpath_1.replace(fpath_tmp)
                fpath_2.replace(fpath_1)
                fpath_tmp.replace(fpath_2)

    def is_complete(self) -> bool:
        return self.results.get('job_done', False)

    def is_converged(self) -> bool:
        # Checks convergence of the calculation
        if 'conv_thr' not in self.parameters:
            raise ValueError('Cannot check convergence when "conv_thr" is not set')

        if 'convergence' not in self.results:
            raise ValueError('Could not locate calculation details to check convergence')

        # Check convergence for both filled and empty, allowing for the possibility
        # of do_outerloop(_empty) = False meaning the calculation is immediately
        # 'converged'
        do_outers = [self.parameters.do_outerloop, self.parameters.do_outerloop_empty]
        convergence_data = self.results['convergence'].values()
        converged = []
        for do_outer, convergence in zip(do_outers, convergence_data):
            if not do_outer:
                converged.append(True)
            elif len(convergence) == 0:
                return False
            else:
                converged.append(
                    convergence[-1]['delta_E'] < self.parameters.conv_thr * utils.units.Hartree)
        return all(converged)

    def read_ham_xml_files(self, bare=False) -> List[np.ndarray]:
        # Reads all expected hamiltonian XML files
        ham_dir = self.parameters.outdir / f'{self.parameters.prefix}_{self.parameters.ndw}.save/K00001'
        ham_matrix: List[np.ndarray] = []

        for ispin in range(1, self.parameters.nspin + 1):
            # Construct the filename
            filename = 'hamiltonian'
            if bare:
                filename += '0'
            if self.parameters.nspin > 1:
                filename += str(ispin)
            filename += '.xml'

            # Work out the shape of the arrays we expect (important for padded arrays)
            if self.parameters.nspin == 2:
                if ispin == 1:
                    n_filled = self.parameters.nelup
                else:
                    n_filled = self.parameters.neldw
            else:
                n_filled = self.parameters.nelec // 2
            n_empty = self.parameters.get('nbnd', n_filled) - n_filled

            # Read the hamiltonian
            ham_filled = read_ham_file(ham_dir / filename)[:n_filled, :n_filled]

            if self.has_empty_states():
                # Construct the filename
                filename = 'hamiltonian'
                if bare:
                    filename += '0'
                filename += '_emp'
                if self.parameters.nspin > 1:
                    filename += str(ispin)
                filename += '.xml'

                # Read the hamiltonian
                ham_empty = read_ham_file(ham_dir / filename)[:n_empty, :n_empty]
                ham = block_diag(ham_filled, ham_empty)
            else:
                ham = ham_filled

            # Store the hamiltonian
            ham_matrix.append(ham)

        return ham_matrix

    def _ham_pkl_file(self, bare: bool = False) -> Path:
        if bare:
            suffix = '.bare_ham.pkl'
        else:
            suffix = '.ham.pkl'
        return self.directory / (self.prefix + suffix)

    def read_ham_pkl_files(self, bare: bool = False) -> List[np.ndarray]:
        with open(self._ham_pkl_file(bare), 'rb') as fd:
            ham_matrix = pickle.load(fd)
        return ham_matrix

    def write_ham_pkl_files(self, ham_matrix: List[np.ndarray], bare: bool = False) -> None:
        with open(self._ham_pkl_file(bare), 'wb') as fd:
            pickle.dump(ham_matrix, fd)

    def read_ham_files(self, bare: bool = False) -> List[np.ndarray]:
        # While the hamiltonian information is stored in xml files inside the outdir of the corresponding calculations,
        # we want a workflow to be able to be reconstructed even if these outdirs have been deleted. This means that we
        # need to store the contents of these xml files elsewhere. We do these as python-readable pickle files

        if self._ham_pkl_file(bare).exists():
            ham_matrix = self.read_ham_pkl_files(bare)
        else:
            ham_matrix = self.read_ham_xml_files(bare)
            self.write_ham_pkl_files(ham_matrix, bare)

        return ham_matrix

    def read_results(self):
        super().read_results()

        self.results['lambda'] = self.read_ham_files()
        if self.parameters.do_bare_eigs:
            self.results['bare lambda'] = self.read_ham_files(bare=True)

    @property
    def alphas(self) -> List[List[float]]:
        if not hasattr(self, '_alphas'):
            raise AttributeError(f'{self}.alphas has not been initialised')
        return self._alphas

    @alphas.setter
    def alphas(self, val: Union[Series, List[List[float]]]):

        if isinstance(val, Series):
            val = val.to_numpy()

        assert len(val) == self.parameters.nspin, \
            f'Dimensions of {self.__class__.__name__}.alphas must match nspin = {self.parameters.nspin}'

        self._alphas = val

    @property
    def filling(self) -> List[List[bool]]:
        # Filling is indexed by [i_spin, i_orbital]
        # Filling is written in this way such that we can calculate it automatically,
        # but if it is explicitly set then we will always use that value instead
        if not hasattr(self, '_filling'):
            filling = []

            # Work out how many filled and empty bands we will have for each spin channel
            if self.parameters.nspin == 2:
                nel_list = [self.parameters.nelup, self.parameters.neldw]
            else:
                nel_list = [self.parameters.nelec // 2]

            # Generate the filling list
            for nel in nel_list:
                if 'nbnd' in self.parameters:
                    nemp = self.parameters.nbnd - nel
                else:
                    nemp = 0
                filling.append([True for _ in range(nel)] + [False for _ in range(nemp)])

            self._filling = filling
        return self._filling

    @filling.setter
    def filling(self, val: List[List[bool]]):
        assert len(val) == self.parameters.nspin, \
            f'Dimensions of {self.__class__.__name__}.filling must match nspin = {self.parameters.nspin}'
        self._filling = val

    def write_alphas(self):
        '''
        Generates file_alpharef.txt and file_alpharef_empty.txt

        '''

        if not self.parameters.do_orbdep or not self.parameters.odd_nkscalfact:
            return

        flat_alphas = [a for sublist in self.alphas for a in sublist]
        flat_filling = [f for sublist in self.filling for f in sublist]
        utils.write_alpha_file(self.directory, flat_alphas, flat_filling)

    def read_alphas(self) -> List[List[float]]:
        '''
        Reads in file_alpharef.txt and file_alpharef_empty.txt from this calculation's directory

        Output:
           alphas -- a list of alpha values (1 per orbital)
        '''

        if not self.parameters.do_orbdep or not self.parameters.odd_nkscalfact:
            return [[]]

        flat_alphas = utils.read_alpha_file(self.directory)

        return convert_flat_alphas_for_kcp(flat_alphas, self.parameters)

    # The following functions enable DOS generation via ase.dft.dos.DOS(<KoopmansCPCalculator object>)
    def get_k_point_weights(self):
        return [1]

    def get_number_of_spins(self):
        return 1

    def get_eigenvalues(self, kpt=None, spin=0):
        if 'eigenvalues' not in self.results:
            raise ValueError('You must first perform a calculation before you try to access the KS eigenvalues')

        if kpt is None:
            return [self.results['eigenvalues'][spin]]
        elif kpt == 0:
            return self.results['eigenvalues'][spin]
        else:
            raise ValueError(f'{self.__class__.__name__} does not have k-point-resolved KS eigenvalues')

    def get_fermi_level(self):
        return 0

    def has_empty_states(self, spin: Optional[int] = None) -> bool:
        if 'nbnd' not in self.parameters:
            return False
        if self.parameters.nspin == 1:
            nel = self.parameters.nelec // 2
        elif spin == 0:
            nel = self.parameters.nelup
        elif spin == 1:
            nel = self.parameters.neldw
        elif 'nelup' in self.parameters and 'neldw' in self.parameters:
            return self.has_empty_states(spin=0) or self.has_empty_states(spin=1)
        else:
            nel = self.parameters.nelec // 2
        return self.parameters.nbnd > nel

    def nspin1_dummy_calculator(self) -> KoopmansCPCalculator:
        calc = copy.deepcopy(self)
        calc.prefix += '_nspin1_dummy'
        calc.parameters.do_outerloop = False
        calc.parameters.do_outerloop_empty = False
        calc.parameters.nspin = 1
        if hasattr(calc, 'alphas'):
            calc.alphas = [calc.alphas[0]]
        if hasattr(calc, 'filling'):
            calc.filling = [calc.filling[0]]
        calc.parameters.nelup = None
        calc.parameters.neldw = None
        calc.parameters.tot_magnetization = None
        calc.parameters.ndw, calc.parameters.ndr = 98, 98
        calc.parameters.restart_mode = 'from_scratch'
        return calc

    def nspin1_calculator(self) -> KoopmansCPCalculator:
        calc = copy.deepcopy(self)
        calc.prefix += '_nspin1'
        calc.parameters.nspin = 1
        calc.parameters.nelup = None
        calc.parameters.neldw = None
        if hasattr(calc, 'alphas'):
            calc.alphas = [calc.alphas[0]]
        if hasattr(calc, 'filling'):
            calc.filling = [calc.filling[0]]
        calc.parameters.tot_magnetization = None
        calc.parameters.ndw, calc.parameters.ndr = 98, 98
        return calc

    def nspin2_dummy_calculator(self) -> KoopmansCPCalculator:
        calc = copy.deepcopy(self)
        calc.prefix += '_nspin2_dummy'
        calc.parameters.restart_mode = 'from_scratch'
        calc.parameters.do_outerloop = False
        calc.parameters.do_outerloop_empty = False
        calc.parameters.ndw = 99
        return calc

    def prepare_to_read_nspin1(self):
        self.prefix += '_nspin2'
        self.parameters.restart_mode = 'restart'
        self.parameters.ndr = 99

    @property
    def from_scratch(self):
        return self.parameters.restart_mode == 'from_scratch'

    def convert_wavefunction_2to1(self):
        nspin2_tmpdir = self.parameters.outdir / f'{self.parameters.prefix}_{self.parameters.ndr}.save/K00001'
        nspin1_tmpdir = self.parameters.outdir / f'{self.parameters.prefix}_98.save/K00001'

        for directory in [nspin2_tmpdir, nspin1_tmpdir]:
            if not directory.is_dir():
                raise OSError(f'{directory} not found')

        for wfile in ['evc0.dat', 'evc0_empty1.dat', 'evcm.dat', 'evc.dat', 'evcm.dat', 'hamiltonian.xml',
                      'eigenval.xml', 'evc_empty1.dat', 'lambda01.dat', 'lambdam1.dat']:
            if '1.' in wfile:
                prefix, suffix = wfile.split('1.')
            else:
                prefix, suffix = wfile.split('.')

            file_out = nspin1_tmpdir / wfile
            file_in = nspin2_tmpdir / f'{prefix}1.{suffix}'

            if file_in.is_file():

                with open(file_in, 'rb') as fd:
                    contents = fd.read()

                contents = contents.replace(b'nk="2"', b'nk="1"')
                contents = contents.replace(b'nspin="2"', b'nspin="1"')

                with open(file_out, 'wb') as fd:
                    fd.write(contents)

    def convert_wavefunction_1to2(self):
        nspin1_tmpdir = self.parameters.outdir / f'{self.parameters.prefix}_98.save/K00001'
        nspin2_tmpdir = self.parameters.outdir / f'{self.parameters.prefix}_99.save/K00001'

        for directory in [nspin2_tmpdir, nspin1_tmpdir]:
            if not directory.is_dir():
                raise OSError(f'{directory} not found')

        for wfile in ['evc0.dat', 'evc0_empty1.dat', 'evcm.dat', 'evc.dat', 'evcm.dat', 'hamiltonian.xml',
                      'eigenval.xml', 'evc_empty1.dat', 'lambda01.dat']:
            if '1.' in wfile:
                prefix, suffix = wfile.split('1.')
            else:
                prefix, suffix = wfile.split('.')

            file_in = nspin1_tmpdir / wfile

            if file_in.is_file():
                with open(file_in, 'rb') as fd:
                    contents = fd.read()

                contents = contents.replace(b'nk="1"', b'nk="2"')
                contents = contents.replace(b'nspin="1"', b'nspin="2"')

                file_out = nspin2_tmpdir / f'{prefix}1.{suffix}'
                with open(file_out, 'wb') as fd:
                    fd.write(contents)

                contents = contents.replace(b'ik="1"', b'ik="2"')
                contents = contents.replace(b'ispin="1"', b'ispin="2"')

                file_out = nspin2_tmpdir / f'{prefix}2.{suffix}'
                with open(file_out, 'wb') as fd:
                    fd.write(contents)


def convert_flat_alphas_for_kcp(flat_alphas: List[float],
                                parameters: settings.KoopmansCPSettingsDict) -> List[List[float]]:
    # Read alpha file returns a flat list ordered by filled spin up, filled spin down, empty spin up, empty spin down
    # Here we reorder this into a nested list indexed by [i_spin][i_orbital]
    if parameters.nspin == 2:
        nbnd = len(flat_alphas) // 2
        alphas = [flat_alphas[:parameters.nelup]
                  + flat_alphas[parameters.nelec:-(nbnd - parameters.neldw)],
                  flat_alphas[parameters.nelup:parameters.nelec]
                  + flat_alphas[-(nbnd - parameters.neldw):]]
    else:
        alphas = [flat_alphas]
    return alphas
