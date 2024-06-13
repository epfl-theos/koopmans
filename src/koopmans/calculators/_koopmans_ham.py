"""

KCWHam calculator module for koopmans

Written by Edward Linscott Feb 2021

"""

import os
from typing import List, Optional

import numpy as np
from aiida_koopmans.helpers import from_kcwham_to_KcwCalculation
from ase import Atoms
from ase.calculators.espresso import KoopmansHam
from ase.dft.kpoints import BandPath

from koopmans import settings, utils
from koopmans.commands import ParallelCommand

from ._utils import CalculatorABC, KCWannCalculator, ReturnsBandStructure


class KoopmansHamCalculator(KCWannCalculator, KoopmansHam, ReturnsBandStructure, CalculatorABC):
    # Subclass of KCWannCalculator for calculating the Koopmans Hamiltonian with kcw.x
    ext_in = '.khi'
    ext_out = '.kho'

    def __init__(self, atoms: Atoms, alphas: Optional[List[int]] = None, *args, **kwargs):
        # Define the valid settings
        self.parameters = settings.KoopmansHamSettingsDict()

        # Initialize using the ASE parent, and then CalculatorExt
        KoopmansHam.__init__(self, atoms=atoms)
        KCWannCalculator.__init__(self, *args, **kwargs)

        self.command = ParallelCommand(f'kcw.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out} 2>&1')

        # Store the alphas
        self.alphas = alphas

    def write_alphas(self):
        # self.alphas is a list of alpha values indexed by spin index and then band index. Meanwhile, kcw.x takes a
        # single file for the alphas (rather than splitting between filled/empty) and does not have two columns for
        # spin up then spin down
        assert self.alphas is not None, 'You have not provided screening parameters to this calculator'
        if not len(self.alphas) == 1:
            raise NotImplementedError('KoopmansHamCalculator yet to be implemented for spin-polarized systems')
        [alphas] = self.alphas
        filling = [True for _ in range(len(alphas))]
        utils.write_alpha_file(self.directory, alphas, filling)

    def _calculate(self):
        self.write_alphas()
        super()._calculate()
        if isinstance(self.parameters.kpts, BandPath) and len(self.parameters.kpts.kpts) > 1:
            # Add the bandstructure to the results
            self.generate_band_structure()

    def get_k_point_weights(self):
        utils.warn('Need to properly define k-point weights')
        return np.ones(len(self.parameters['kpath'].kpts))

    def get_number_of_spins(self):
        return 1

    def get_eigenvalues(self, kpt=None, spin=0):
        if spin != 0:
            raise NotImplementedError(
                f'Koopmans Hamiltonian calculator is not implemented for spin-polarized systems')

        if 'band structure' not in self.results:
            raise ValueError('You must first calculate the band structure before you try to access the KS eigenvalues')

        if kpt is None:
            return self.results['band structure'].energies[spin, :]
        else:
            return self.results['band structure'].energies[spin, kpt]

    def get_fermi_level(self):
        return 0

    def vbm_energy(self) -> float:
        eigenvalues_np = self.eigenvalues_from_results()
        return np.max(eigenvalues_np[:, :, self.parameters.num_wann_occ - 1])

    def eigenvalues_from_results(self):
        assert 'eigenvalues' in self.results, 'Please call {0}.calculate() prior to calling {0}.band_structure'.format(
            self.__class__.__name__)

        return np.array([self.results['eigenvalues']])

    def is_converged(self):
        raise NotImplementedError('TODO')

    def check_convergence(self) -> None:
        # is_converged has not been implemented yet for this calculator
        return

    def read_input(self, **kwargs):
        # A .khi file doesn't have the requisite information to reconstruct the bandpath, so in the event that kpts
        # are already provided in self.parameters, don't overwrite them

        kpts = self.parameters.kpts

        super().read_input(**kwargs)

        if kpts is not None:
            self.parameters.kpts = kpts
        return
    
    # MB mod:
    def read_results(self, wchain=None):
        from ase import io
        if not wchain:
            output = io.read(self.label + '.kho')
        else:
            import pathlib
            import tempfile
            
            # Create temporary directory. However, see aiida-wannier90-workflows/src/aiida_wannier90_workflows/utils/workflows/pw.py for more advanced and smart ways.
            retrieved = wchain.outputs.retrieved
            with tempfile.TemporaryDirectory() as dirpath:
                # Open the output file from the AiiDA storage and copy content to the temporary file
                for filename in retrieved.base.repository.list_object_names():
                    if '.out' in filename:
                        # Create the file with the desired name
                        readable_filename = "kc.kho"
                        temp_file = pathlib.Path(dirpath) / readable_filename
                        with retrieved.open(filename, 'rb') as handle:
                            temp_file.write_bytes(handle.read())
                    
                        output = io.read(temp_file)
        
        self.calc = output.calc
        self.results = output.calc.results
        if hasattr(output.calc, 'kpts'):
            self.kpts = output.calc.kpts
            
    def calculate(self):

        if self.mode == "ase":
            return super().calculate()
        else: # MB mod
            builder = from_kcwham_to_KcwCalculation(kcw_calculator=self)
            from aiida.engine import run_get_node, submit
            running = run_get_node(builder)
            
            # once the running if completed
            self.calculation = running[-1]
            
            # with the correct filename (see the method above).
            self.read_results(wchain=self.calculation)
