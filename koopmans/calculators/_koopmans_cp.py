"""

kcp calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os
import math
import numpy as np
from pathlib import Path
from scipy.linalg import block_diag
from typing import Optional, List
from pandas.core.series import Series
import xml.etree.ElementTree as ET
from ase import Atoms
from ase.calculators.espresso import Espresso_kcp
from koopmans import utils, settings, pseudopotentials
from koopmans.commands import ParallelCommand
from ._utils import CalculatorExt, CalculatorABC, kcp_bin_directory


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


class KoopmansCPCalculator(CalculatorExt, Espresso_kcp, CalculatorABC):
    # Subclass of CalculatorExt for performing calculations with kcp.x
    ext_in = '.cpi'
    ext_out = '.cpo'

    def __init__(self, atoms: Atoms, skip_qc: bool = False, alphas: Optional[List[float]] = None,
                 filling: Optional[List[bool]] = None, **kwargs):
        # Define the valid parameters
        self.parameters = settings.KoopmansCPSettingsDict()

        # Initialise first using the ASE parent and then CalculatorExt
        Espresso_kcp.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, skip_qc, **kwargs)

        # Add nelec if it is missing
        if 'nelec' not in self.parameters and 'pseudopotentials' in self.parameters:
            self.parameters.nelec = pseudopotentials.nelec_from_pseudos(
                self.atoms, self.pseudopotentials, self.parameters.pseudo_dir)

        if not isinstance(self.command, ParallelCommand):
            self.command = ParallelCommand(os.environ.get(
                'ASE_ESPRESSO_KCP_COMMAND', str(kcp_bin_directory) + os.path.sep + self.command))

        self.alphas = alphas
        self.filling = filling
        self.results_for_qc = ['energy', 'homo_energy', 'lumo_energy']

    def is_complete(self):
        return self.results.get('job_done', False)

    def is_converged(self):
        # Checks convergence of the calculation
        if 'conv_thr' not in self.parameters:
            raise ValueError('Cannot check convergence when "conv_thr" is not set')
        return self._ase_is_converged()

    def read_ham_files(self, bare=False) -> List[np.ndarray]:
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

            # Read the hamiltonian
            ham_filled = read_ham_file(ham_dir / filename)

            if self.parameters.empty_states_nbnd > 0:
                # Construct the filename
                filename = 'hamiltonian'
                if bare:
                    filename += '0'
                filename += '_emp'
                if self.parameters.nspin > 1:
                    filename += str(ispin)
                filename += '.xml'

                # Read the hamiltonian
                ham_empty = read_ham_file(ham_dir / filename)
                ham = block_diag(ham_filled, ham_empty)
            else:
                ham = ham_filled

            # Store the hamiltonian
            ham_matrix.append(ham)

        return ham_matrix

    def read_results(self):
        return_val = super().read_results()

        self.results['lambda'] = self.read_ham_files()
        if self.parameters.do_bare_eigs:
            self.results['bare lambda'] = self.read_ham_files(bare=True)

        return return_val

    def _ase_is_converged(self):
        if 'convergence' not in self.results:
            raise ValueError(
                'Could not locate calculation details to check convergence')

        # Check convergence for both filled and empty, allowing for the possibility
        # of do_outerloop(_empty) = False meaning the calculation is immediately
        # 'converged'
        do_outers = [self.parameters.do_outerloop, self.parameters.do_outerloop_empty]
        convergence_data = self.results['convergence'].values()
        converged = []
        for do_outer, convergence in zip(do_outers, convergence_data):
            if len(convergence) == 0:
                return False
            if not do_outer:
                converged.append(True)
            else:
                converged.append(
                    convergence[-1]['delta_E'] < self.parameters.conv_thr * utils.units.Hartree)
        return all(converged)

    @property
    def alphas(self):
        return self._alphas

    @alphas.setter
    def alphas(self, val):

        if val is None:
            return
        elif isinstance(val, Series):
            val = val.to_numpy()

        # alphas can be a 1D or 2D array. If it is 1D, convert it to 2D so that self.alphas
        # is indexed by [i_spin, i_orbital]
        if isinstance(val[0], float):
            val = [val for _ in range(self.parameters.nspin)]

        if len(val) == 1 and self.nspin == 2:
            val = [val[0] for _ in range(self.parameters.nspin)]

        self._alphas = val

    @property
    def filling(self):
        # Filling is indexed by [i_spin, i_orbital]
        # Filling is written in this way such that we can calculate it automatically,
        # but if it is explicitly set then we will always use that value instead
        if self._filling is None:
            n_filled_bands = self.parameters.nelec // 2
            n_empty_bands = self.parameters.empty_states_nbnd
            if n_empty_bands is None:
                n_empty_bands = 0
            filled_spin_channel = [True for _ in range(n_filled_bands)] + [False for _ in range(n_empty_bands)]
            return [filled_spin_channel for _ in range(self.parameters.nspin)]
        else:
            return self._filling

    @filling.setter
    def filling(self, val):

        if val is not None:
            # val can be a 1D or 2D array. If it is 1D, convert it to 2D so that filling
            # is indexed by [i_spin, i_orbital]
            if isinstance(val[0], bool):
                val = [val for _ in range(self.parameters.nspin)]

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

    def read_alphas(self):
        '''
        Reads in file_alpharef.txt and file_alpharef_empty.txt from this calculation's directory

        Output:
           alphas -- a list of alpha values (1 per orbital)
        '''

        if not self.parameters.do_orbdep or not self.parameters.odd_nkscalfact:
            return

        alphas = utils.read_alpha_file(self.directory)

        if self.parameters.nspin == 2:
            # Remove duplicates
            alphas = alphas[:self.parameters.nelec // 2] + \
                alphas[self.parameters.nelec:self.parameters.nelec + self.parameters.empty_states_nbnd]
            alphas = [alphas, alphas]
        return alphas

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
