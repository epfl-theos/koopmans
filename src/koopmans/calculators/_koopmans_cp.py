"""kcp calculator module for koopmans."""

from __future__ import annotations

import copy
import math
import os
import pickle
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
from ase_koopmans import Atoms
from ase_koopmans.calculators.espresso import Espresso_kcp
from pandas.core.series import Series
from scipy.linalg import block_diag

from koopmans import bands, pseudopotentials, settings, utils
from koopmans.cell import cell_follows_qe_conventions, cell_to_parameters
from koopmans.commands import ParallelCommand
from koopmans.files import File

from ._calculator import (CalculatorABC, CalculatorCanEnforceSpinSym,
                          CalculatorExt)


def allowed(nr: int) -> bool:
    """Return whether i is a good fft grid number."""
    if nr < 1:
        return False
    mr = nr
    factor = [2, 3, 5, 7, 11]
    allowed = False
    pwr = [0, 0, 0, 0, 0]
    for i, fac in enumerate(factor):
        maxpwr = int(np.log(mr) / np.log(fac))
        for _ in range(maxpwr):
            if mr == 1:
                break
            if mr % fac == 0:
                mr //= fac
                pwr[i] += 1
    if mr != 1:
        allowed = False
    else:
        allowed = (pwr[3] == 0) and (pwr[4] == 0)
    return allowed


def good_fft(nr: int) -> int:
    """Return good grid dimension (optimal for the FFT)."""
    nfftx = 2049
    new = nr
    while allowed(new) is False and (new <= nfftx):
        new = new + 1
    nr = new
    return nr


def read_ham_file(filename: Path) -> np.ndarray[Any, np.dtype[np.complex128]]:
    """Read a single hamiltonian XML file."""
    if not filename.exists():
        raise FileExistsError(f'`{filename}` does not exist')

    with open(filename, 'r') as fd:
        tree = ET.parse(fd)
    ham_xml = tree.getroot()

    length = int(math.sqrt(int(ham_xml.attrib['size'])))

    assert ham_xml.text is not None, f'{filename} is empty'

    ham_array = np.array([complex(*[float(x) for x in line.split(',')])
                          for line in ham_xml.text.strip().split('\n')], dtype=np.complex128) * utils.units.Hartree

    return ham_array.reshape((length, length))


class KoopmansCPCalculator(CalculatorCanEnforceSpinSym, CalculatorExt, Espresso_kcp, CalculatorABC):
    """Subclass of CalculatorExt for performing calculations with kcp.x."""

    ext_in = '.cpi'
    ext_out = '.cpo'

    def __init__(self, atoms: Atoms, alphas: Optional[List[List[float]]] = None,
                 filling: Optional[List[List[bool]]] = None, **kwargs):

        self.parent_process = None

        # Define the valid parameters
        self.parameters = settings.KoopmansCPSettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        Espresso_kcp.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, **kwargs)

        # Add nelec, nelup, neldw if they are missing
        if 'nelec' not in self.parameters and 'pseudopotentials' in self.parameters:
            self.parameters.nelec = pseudopotentials.nelec_from_pseudos(self.atoms, self.parameters.pseudopotentials)
        if 'nelec' in self.parameters:
            if 'nelup' not in self.parameters:
                self.parameters.nelup = self.parameters.nelec // 2
            if 'neldw' not in self.parameters:
                self.parameters.neldw = self.parameters.nelec // 2

        if not isinstance(self.command, ParallelCommand):
            self.command = ParallelCommand(os.environ.get('ASE_ESPRESSO_KCP_COMMAND', self.command))

        if alphas is not None:
            self.alphas = alphas
        if filling is not None:
            self.filling = filling

        # Give the calculator an attribute to keep track of which band has been held fixed, for calculations where
        # fixed_state = .true.. N.B. this differs from self.parameters.fixed_band in the case of empty orbitals (see
        # koopmans.workflows._koopmans_dscf.py for more details)
        self.fixed_band: Optional[bands.Band] = None

        # Create a private attribute to keep track of whether the spin channels have been swapped
        self._spin_channels_are_swapped: bool = False

    def _pre_calculate(self):
        # kcp.x imposes nelup >= neldw, so if we try to run a calcualtion with neldw > nelup, swap the spin channels
        # This needs to happen before super()._pre_calculate() because this is when the linked files are fetched
        if self.parameters.nspin == 2:
            self._spin_channels_are_swapped = self.parameters.nelup < self.parameters.neldw
            # Swap the spin channels if required
            if self._spin_channels_are_swapped:
                self._swap_spin_channels()

        super()._pre_calculate()

        # Write out screening parameters to file
        if self.parameters.get('do_orbdep', False):
            self.write_alphas()

        # Update ibrav and celldms
        if cell_follows_qe_conventions(self.atoms.cell):
            self.parameters.update(**cell_to_parameters(self.atoms.cell))
        else:
            self.parameters.ibrav = 0

        # Autogenerate the nr keywords
        self._autogenerate_nr()

    def _post_calculate(self):

        super()._post_calculate()

        # Check spin-up and spin-down eigenvalues match
        if 'eigenvalues' in self.results and self.parameters.do_outerloop \
                and self.parameters.nspin == 2 and self.parameters.tot_magnetization == 0 \
                and not self.parameters.fixed_state and len(self.results['eigenvalues']) > 0:
            rms_eigenval_difference = np.sqrt(np.mean(np.diff(self.results['eigenvalues'], axis=0)**2))
            if rms_eigenval_difference > 0.05:
                utils.warn('Spin-up and spin-down eigenvalues differ substantially')

        # Swap the spin channels back
        if self._spin_channels_are_swapped:
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

        # Linked files

        # First, recurse over the linked files and replace directories with their contents
        for key in list(self.linked_files.keys()):
            dst_filepointer, sym, rsym, force = self.linked_files[key]
            if dst_filepointer.is_dir():
                self.linked_files.pop(key)
                for subfile in dst_filepointer.rglob('*'):
                    if subfile.is_dir():
                        continue
                    self.linked_files[str(subfile.name)] = (subfile, rsym, False, force)

        # Now iterate over the linked files and swap any pairs of spin up/down files
        for dst_file in self.linked_files:
            if dst_file.endswith('1.dat'):
                other_dst_file = dst_file.replace('1.dat', '2.dat')
                if other_dst_file not in self.linked_files:
                    raise FileNotFoundError(f'Expected {other_dst_file} to be linked to the {self.prefix} calculator')
                # Switch the links
                self.linked_files[dst_file], self.linked_files[other_dst_file] = self.linked_files[other_dst_file], \
                    self.linked_files[dst_file]

        # alphas and filling
        self.alphas = self.alphas[::-1]
        self.filling = self.filling[::-1]

        # Results
        if 'orbital_data' in self.results:
            self.results['orbital_data'] = {k: v[::-1] for k, v in self.results['orbital_data'].items()}
        for key in ['eigenvalues', 'lambda']:
            if key in self.results:
                self.results[key] = self.results[key][::-1]

    def _autogenerate_nr(self):
        """Autogenerate the nr parameters.

        For norm-conserving pseudopotentials the small box grid (nr1b, nr2b, nr3b) is needed in case the pseudo has
        non-linear core corrections. This function automatically defines this small box using a conservative guess.
        """
        has_nlcc = False
        for p in self.parameters.pseudopotentials.values():
            upf = pseudopotentials.read_pseudo_file(self.directory / self.parameters.pseudo_dir / p)
            if upf['header']['core_correction']:
                has_nlcc = True
        if has_nlcc and (self.parameters.nr1b is None or self.parameters.nr2b is None or self.parameters.nr3b is None):
            # Extract alat (in Bohr)
            if cell_follows_qe_conventions(self.atoms.cell):
                params = cell_to_parameters(self.atoms.cell)
                alat = params['celldms'][1]
            else:
                alat = np.linalg.norm(self.atoms.cell[0]) / utils.units.Bohr

            # Define reduced lattice vectors ("at" in espresso)
            at = self.atoms.cell / (alat * utils.units.Bohr)

            # nr1 = int ( sqrt (gcutm) * sqrt (at(1, 1)**2 + at(2, 1)**2 + at(3, 1)**2) ) + 1
            ecutrho = self.parameters.get('ecutrho', 4 * self.parameters.ecutwfc)
            [nr1, nr2, nr3] = [2 * int(np.sqrt(ecutrho) / (2.0 * np.pi / alat) * np.linalg.norm(vec) + 1) for vec in at]

            # set good_fft dimensions
            nr1 = good_fft(nr1)
            nr2 = good_fft(nr2)
            nr3 = good_fft(nr3)

            # At this stage nr should match with the one generated by CP or PW for the charge density. This can be a
            # very safe choice for nrb, but most probably too conservative. If we have access to the pseudopotential
            # cutoff radius we can define a more reasonable one as nr1b = nr1 * (2*rc/L1) where rc is the cutoff
            # radius used to generate the PP and L1 is the dimension of the simulation Box.
            # N.B. we assume here 3 Bohr is a safe choice; all the rc in the DOJO pseudos are <= 2.6:
            rc_safe = 3.0
            [nr1b, nr2b, nr3b] = [int(nr * 2 * rc_safe / (np.linalg.norm(vec) * alat))
                                  for vec, nr in zip(at, [nr1, nr2, nr3])]

            self.parameters.nr1b = good_fft(nr1b)
            self.parameters.nr2b = good_fft(nr2b)
            self.parameters.nr3b = good_fft(nr3b)

            utils.warn('Small box parameters `nrb` not provided in input: these will be automatically set to safe '
                       'default values. These values can probably be decreased, but this would require convergence '
                       f'tests. Estimated real mesh dimension `(nr1, nr2, nr3) = {nr1} {nr2} {nr3}`. Small box '
                       f'mesh dimension `(nr1b, nr2b, nr3b) = {self.parameters.nr1b} {self.parameters.nr2b} '
                       f'{self.parameters.nr3b}`.')

    def is_complete(self) -> bool:
        """Check if the calculation is complete."""
        return self.results.get('job_done', False)

    def is_converged(self) -> bool:
        """Check the convergence of the calculation."""
        if 'conv_thr' not in self.parameters:
            raise ValueError('Cannot check convergence when `conv_thr` is not set')

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
        """Read all the expected Hamiltonian XML files."""
        ham_dir = self.directory / self.parameters.outdir / \
            f'{self.parameters.prefix}_{self.parameters.ndw}.save/K00001'
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
        """Return the path to the pickle file containing the Hamiltonian."""
        assert self.directory is not None
        if bare:
            suffix = '.bare_ham.pkl'
        else:
            suffix = '.ham.pkl'
        return self.directory / (self.prefix + suffix)

    def read_ham_pkl_files(self, bare: bool = False) -> List[np.ndarray]:
        """Read the Hamiltonian from pickle files."""
        with open(self._ham_pkl_file(bare), 'rb') as fd:
            ham_matrix = pickle.load(fd)
        return ham_matrix

    def write_ham_pkl_files(self, ham_matrix: List[np.ndarray], bare: bool = False) -> None:
        """Write the Hamiltonian to pickle files."""
        with open(self._ham_pkl_file(bare), 'wb') as fd:
            pickle.dump(ham_matrix, fd)

    def read_ham_files(self, bare: bool = False) -> List[np.ndarray]:
        """Read the hamiltonian files from the calculation.

        While the hamiltonian information is stored in xml files inside the outdir of the corresponding calculations,
        we want a workflow to be able to be reconstructed even if these outdirs have been deleted. This means that we
        need to store the contents of these xml files elsewhere. We do these as python-readable pickle files
        """
        if self._ham_pkl_file(bare).exists():
            ham_matrix = self.read_ham_pkl_files(bare)
        else:
            ham_matrix = self.read_ham_xml_files(bare)
            self.write_ham_pkl_files(ham_matrix, bare)

        return ham_matrix

    def read_results(self):
        """Read the results of the calculation from the output file."""
        super().read_results()

        self.results['lambda'] = self.read_ham_files()
        if self.parameters.do_bare_eigs:
            self.results['bare lambda'] = self.read_ham_files(bare=True)

    @property
    def alphas(self) -> List[List[float]]:
        """Return the screening parameters for the calculation."""
        if not hasattr(self, '_alphas'):
            raise AttributeError(f'`{self}.alphas` has not been initialized')
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
        """Return the filling of each orbital the calculation.

        Filling is indexed by [i_spin, i_orbital]
        Filling is written in this way such that we can calculate it automatically,
        but if it is explicitly set then we will always use that value instead
        """
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
        """Generate file_alpharef.txt and file_alpharef_empty.txt."""
        if not self.parameters.do_orbdep or not self.parameters.odd_nkscalfact:
            return

        flat_alphas = [a for sublist in self.alphas for a in sublist]
        flat_filling = [f for sublist in self.filling for f in sublist]
        utils.write_alpha_file(self, flat_alphas, flat_filling)

    def read_alphas(self) -> List[List[float]]:
        """Read in file_alpharef.txt and file_alpharef_empty.txt from this calculation's directory.

        Output:
           alphas -- a list of alpha values (1 per orbital)
        """
        if not self.parameters.do_orbdep or not self.parameters.odd_nkscalfact:
            return [[]]

        assert self.directory is not None
        flat_alphas = utils.read_alpha_file(self)

        assert isinstance(self.parameters, settings.KoopmansCPSettingsDict)

        return convert_flat_alphas_for_kcp(flat_alphas, self.parameters)

    def get_k_point_weights(self):
        """Return the k-point weights for the calculation.

        This is a dummy function because it is assumed we have the Gamma-point only. We need to define this function
        to enable the DOS generation via ase.dft.dos.DOS(<KoopmansCPCalculator object>).
        """
        return [1]

    def get_number_of_spins(self):
        """Return the number of spins in the calculation."""
        return 1

    def get_eigenvalues(self, kpt=None, spin=0):
        """Return the eigenvalues of the calculation."""
        if 'eigenvalues' not in self.results:
            raise ValueError('You must first perform a calculation before you try to access the KS eigenvalues')

        if kpt is None:
            return [self.results['eigenvalues'][spin]]
        elif kpt == 0:
            return self.results['eigenvalues'][spin]
        else:
            raise ValueError(f'`{self.__class__.__name__}` does not have k-point-resolved KS eigenvalues')

    def get_fermi_level(self):
        """Return the Fermi level."""
        return 0

    def has_empty_states(self, spin: Optional[int] = None) -> bool:
        """Return True if the calculation has empty states."""
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
        """Create a copy of the calculator that is set up to run a dummy nspin=1 calculation."""
        self.parent_process, parent_process = None, self.parent_process
        self.linked_files, linked_files = {}, self.linked_files
        calc = copy.deepcopy(self)
        calc.linked_files = linked_files
        self.linked_files = linked_files
        self.parent_process = parent_process
        calc.parent_process = parent_process
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
        """Create a copy of the calculator that is set up to run a nspin=1 calculation."""
        self.parent_process, parent_process = None, self.parent_process
        self.linked_files, linked_files = {}, self.linked_files
        calc = copy.deepcopy(self)
        calc.linked_files = linked_files
        self.linked_files = linked_files
        calc.parent_process = parent_process
        self.parent_process = parent_process
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
        """Create a copy of the calculator that is set up to run a nspin=2 calculation from scratch."""
        calc = copy.deepcopy(self)
        calc.prefix += '_nspin2_dummy'
        calc.parameters.restart_mode = 'from_scratch'
        calc.parameters.do_outerloop = False
        calc.parameters.do_outerloop_empty = False
        calc.parameters.ndw = 99
        return calc

    def prepare_to_read_nspin1(self):
        """Set up the calculation to read from a nspin=1 calculation."""
        self.prefix += '_nspin2'
        self.parameters.restart_mode = 'restart'
        self.parameters.ndr = 99

    @property
    def from_scratch(self):
        """Return True if the calculation is running from scratch."""
        return self.parameters.restart_mode == 'from_scratch'

    @property
    def files_to_convert_with_spin2_to_spin1(self) -> Dict[str, List[File] | List[Path]]:
        """Return a list of files that need to be converted when going from spin 2 to spin 1."""
        nspin_2_files = []
        nspin_1_files = []
        for f in ['evc0.dat', 'evc0_empty1.dat', 'evcm.dat', 'evc.dat', 'evcm.dat', 'hamiltonian.xml',
                  'eigenval.xml', 'evc_empty1.dat', 'lambda01.dat', 'lambdam1.dat']:
            if '1.' in f:
                prefix, suffix = f.split('1.')
            else:
                prefix, suffix = f.split('.')

            nspin_2_file = self.read_directory / 'K00001' / f'{prefix}1.{suffix}'
            if nspin_2_file.exists():
                nspin_2_files.append(nspin_2_file)
                nspin_1_files.append(Path(f))

        return {'spin_2_files': nspin_2_files, 'spin_1_files': nspin_1_files}

    @property
    def files_to_convert_with_spin1_to_spin2(self):
        """Return a list of files that need to be converted when going from spin 1 to spin 2."""
        nspin_1_files = []
        nspin_2up_files = []
        nspin_2dw_files = []

        for nspin_1_file in ['evc0.dat', 'evc0_empty1.dat', 'evcm.dat', 'evc.dat', 'evcm.dat', 'hamiltonian.xml',
                             'eigenval.xml', 'evc_empty1.dat', 'lambda01.dat']:

            if '1.' in nspin_1_file:
                prefix, suffix = nspin_1_file.split('1.')
            else:
                prefix, suffix = nspin_1_file.split('.')

            nspin_1_file = self.read_directory / 'K00001' / nspin_1_file
            if nspin_1_file.exists():
                nspin_1_files.append(nspin_1_file)
                nspin_2up_files.append(f'{prefix}1.{suffix}')
                nspin_2dw_files.append(f'{prefix}2.{suffix}')

        return {'spin_1_files': nspin_1_files,
                'spin_2_up_files': nspin_2up_files,
                'spin_2_down_files': nspin_2dw_files}

    @property
    def read_directory(self) -> File:
        """Return the directory where the calculation reads from."""
        assert isinstance(self.parameters.outdir, Path)
        assert self.parameters.ndr is not None
        assert self.parameters.prefix is not None
        return File(self, self.parameters.outdir / f'{self.parameters.prefix}_{self.parameters.ndr}.save')

    @property
    def write_directory(self) -> File:
        """Return the directory where the calculation writes to."""
        assert isinstance(self.parameters.outdir, Path)
        assert self.parameters.ndw is not None
        assert self.parameters.prefix is not None
        return File(self, self.parameters.outdir / f'{self.parameters.prefix}_{self.parameters.ndw}.save')


def convert_flat_alphas_for_kcp(flat_alphas: List[float],
                                parameters: settings.KoopmansCPSettingsDict) -> List[List[float]]:
    """Convert a flat list of alpha values into a nested list indexed by [i_spin][i_orbital].

    Read alpha file returns a flat list ordered by filled spin up, filled spin down, empty spin up, empty spin down
    Here we reorder this into a nested list indexed by [i_spin][i_orbital]"
    """
    if parameters.nspin == 2:
        nbnd = len(flat_alphas) // 2
        nelec = parameters.nelec if parameters.nelec else parameters.nelup + parameters.neldw
        alphas = [flat_alphas[:parameters.nelup]
                  + flat_alphas[nelec:(nbnd + parameters.neldw)],
                  flat_alphas[parameters.nelup:nelec]
                  + flat_alphas[(nbnd + parameters.neldw):]]
    else:
        alphas = [flat_alphas]
    return alphas
