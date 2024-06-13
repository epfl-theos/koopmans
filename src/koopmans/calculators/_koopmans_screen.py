"""

KCWScreen calculator module for koopmans

Written by Edward Linscott Feb 2021

"""

import os

import numpy as np
from aiida_koopmans.helpers import from_kcwscreen_to_KcwCalculation
from ase import Atoms
from ase.calculators.espresso import KoopmansScreen

from koopmans import settings, utils
from koopmans.commands import ParallelCommandWithPostfix

from ._utils import CalculatorABC, KCWannCalculator


class KoopmansScreenCalculator(KCWannCalculator, KoopmansScreen, CalculatorABC):
    # Subclass of KCWannCalculator for calculating screening parameters with kcw.x
    ext_in = '.ksi'
    ext_out = '.kso'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = settings.KoopmansScreenSettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        KoopmansScreen.__init__(self, atoms=atoms)
        KCWannCalculator.__init__(self, *args, **kwargs)
        super().__init__(*args, **kwargs)

        self.command = ParallelCommandWithPostfix(f'kcw.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out} 2>&1')

    # MB mod:
    def read_results(self, wchain=None):
        from ase import io
        if not wchain:
            output = io.read(self.label + '.kso')
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
                        readable_filename = "ks.kso"
                        temp_file = pathlib.Path(dirpath) / readable_filename
                        with retrieved.open(filename, 'rb') as handle:
                            temp_file.write_bytes(handle.read())
                    
                        output = io.read(temp_file)
        
        self.calc = output.calc
        self.results = output.calc.results
        if hasattr(output.calc, 'kpts'):
            self.kpts = output.calc.kpts
        
    def calculate(self):
        # Check eps infinity
        kpoints = [self.parameters.mp1, self.parameters.mp2, self.parameters.mp3]
        if np.max(kpoints) > 1 and self.parameters.eps_inf is None:
            utils.warn('You have not specified a value for eps_inf. This will mean that the screening parameters will '
                       'converge very slowly with respect to the k- and q-point grids')

        if self.mode == "ase":
            return super().calculate()
        else: # MB mod
            builder = from_kcwscreen_to_KcwCalculation(kcw_calculator=self)
            from aiida.engine import run_get_node, submit
            running = run_get_node(builder)
            
            # once the running if completed
            self.calculation = running[-1]
            
            #  with the correct filename (see the method above).
            self.read_results(wchain=self.calculation)

    def is_converged(self):
        raise NotImplementedError('TODO')

    def check_convergence(self) -> None:
        # is_converged has not been implemented yet for this calculator
        return
