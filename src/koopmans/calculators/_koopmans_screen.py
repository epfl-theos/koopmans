"""

KCWScreen calculator module for koopmans

Written by Edward Linscott Feb 2021

"""

import os

import numpy as np
from ase import Atoms
from ase.calculators.espresso import KoopmansScreen

from koopmans import settings, utils
from koopmans.commands import ParallelCommandWithPostfix

from ._utils import CalculatorABC, KCWannCalculator

def from_kcwscreen_to_KcwCalculation(kcw_calculator):
    
    """
    The input parent folder is meant to be set later, at least for now.
    """
    
    from aiida_koopmans.calculations.kcw import KcwCalculation
    from aiida import orm, load_profile
    load_profile()
    
    builder = KcwCalculation.get_builder()
    
    control_namelist = {
                    'kcw_iverbosity': kcw_calculator.parameters.kcw_iverbosity,
                    'kcw_at_ks'      :kcw_calculator.parameters.kcw_at_ks,
                    'calculation'    :kcw_calculator.parameters.calculation,
                    'lrpa'           :kcw_calculator.parameters.lrpa,
                    'mp1'            :kcw_calculator.parameters.mp1,
                    'mp2'            :kcw_calculator.parameters.mp2,
                    'mp3'            :kcw_calculator.parameters.mp3,
                    'homo_only'      :kcw_calculator.parameters.homo_only,
                    'read_unitary_matrix' : kcw_calculator.parameters.read_unitary_matrix,
                    'l_vcut'         :kcw_calculator.parameters.l_vcut,
                    'spin_component' :kcw_calculator.parameters.spin_component,
                    }

    if any(kcw_calculator.atoms.pbc): control_namelist['assume_isolated'] = "m-t"

    wannier_dict = {
                    "check_ks"       : kcw_calculator.parameters.check_ks,
                    "num_wann_occ"   : kcw_calculator.parameters.num_wann_occ,
                    "num_wann_emp"   : kcw_calculator.parameters.num_wann_emp,
                    "have_empty"     : kcw_calculator.parameters.have_empty,
                    "has_disentangle": kcw_calculator.parameters.has_disentangle,
                }
    
    screening_dict = {
                    'tr2'         : 1e-18,
                    'nmix'        : 4,
                    'niter'       : 33,
                    'check_spread': True,
                }


    kcw_screen_params = {
            "CONTROL":control_namelist,
            "WANNIER":wannier_dict,
            "SCREEN":screening_dict,
        }
    
    builder.parameters = orm.Dict(kcw_screen_params)
    builder.code = orm.load_code(kcw_calculator.mode["kcw_code"])
    builder.metadata = kcw_calculator.mode["metadata"]
    builder.parent_folder = kcw_calculator.parent_folder

    return builder

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

    def calculate(self):
        # Check eps infinity
        kpoints = [self.parameters.mp1, self.parameters.mp2, self.parameters.mp3]
        if np.max(kpoints) > 1 and self.parameters.eps_inf is None:
            utils.warn('You have not specified a value for eps_inf. This will mean that the screening parameters will '
                       'converge very slowly with respect to the k- and q-point grids')

        if self.mode == "ase":
            return super().calculate()
        else:
            builder = from_kcwscreen_to_KcwCalculation(self)
            from aiida.engine import run_get_node,submit
            running = run_get_node(builder)
            
            # once the running if completed
            self.calculation = running[-1]
            
            #self.read_results(wchain=self.calculation)

    def is_converged(self):
        raise NotImplementedError('TODO')

    def check_convergence(self) -> None:
        # is_converged has not been implemented yet for this calculator
        return
