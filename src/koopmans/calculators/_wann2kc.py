"""

Wann2KCW calculator module for koopmans

Written by Edward Linscott Feb 2021

"""

import os

from ase import Atoms
from ase.calculators.espresso import Wann2KC

from koopmans.commands import ParallelCommandWithPostfix
from koopmans.settings import Wann2KCSettingsDict

from ._utils import CalculatorABC, KCWannCalculator

# MB mod:
def from_wann2kc_to_KcwCalculation(wann2kc_calculator):
    
    """
    The input parent folder is meant to be set later, at least for now.
    """
    
    from aiida_koopmans.calculations.kcw import KcwCalculation
    from aiida import orm, load_profile
    load_profile()
    
    builder = KcwCalculation.get_builder()
    
    control_namelist = {
                    'kcw_iverbosity': wann2kc_calculator.parameters.kcw_iverbosity,
                    'kcw_at_ks'      :wann2kc_calculator.parameters.kcw_at_ks,
                    'calculation'    :wann2kc_calculator.parameters.calculation,
                    'lrpa'           :wann2kc_calculator.parameters.lrpa,
                    'mp1'            :wann2kc_calculator.parameters.mp1,
                    'mp2'            :wann2kc_calculator.parameters.mp2,
                    'mp3'            :wann2kc_calculator.parameters.mp3,
                    'homo_only'      :wann2kc_calculator.parameters.homo_only,
                    'read_unitary_matrix' : wann2kc_calculator.parameters.read_unitary_matrix,
                    'l_vcut'         :wann2kc_calculator.parameters.l_vcut,
                    'spin_component' :wann2kc_calculator.parameters.spin_component,
                    }

    if any(wann2kc_calculator.atoms.pbc): control_namelist['assume_isolated'] = "m-t"

    wannier_dict = {
                    "check_ks"       : wann2kc_calculator.parameters.check_ks,
                    "num_wann_occ"   : wann2kc_calculator.parameters.num_wann_occ,
                    "num_wann_emp"   : wann2kc_calculator.parameters.num_wann_emp,
                    "have_empty"     : wann2kc_calculator.parameters.have_empty,
                    "has_disentangle": wann2kc_calculator.parameters.has_disentangle,
                        }


    wann2kcw_params = {
            "CONTROL":control_namelist,
            "WANNIER":wannier_dict,
        }
    
    builder.parameters = orm.Dict(wann2kcw_params)
    builder.code = orm.load_code(wann2kc_calculator.mode["kcw_code"])
    builder.metadata = wann2kc_calculator.mode["metadata"]
    builder.parent_folder = wann2kc_calculator.parent_folder

    return builder

class Wann2KCCalculator(KCWannCalculator, Wann2KC, CalculatorABC):
    # Subclass of KCWannCalculator for converting Wannier functions to a KCW format with kcw.x
    ext_in = '.w2ki'
    ext_out = '.w2ko'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = Wann2KCSettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        Wann2KC.__init__(self, atoms=atoms)
        KCWannCalculator.__init__(self, *args, **kwargs)

        self.command = ParallelCommandWithPostfix(f'kcw.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out} 2>&1')

    def is_converged(self):
        return True

    def is_complete(self):
        return self.results['job_done']
       
    # MB mod:
    def calculate(self):
        
        if self.mode == "ase":
            return super().calculate()
        else:
            builder = from_wann2kc_to_KcwCalculation(self)
            from aiida.engine import run,submit
            running = run(builder)
            
            # once the running if completed
            self.calculation = running['remote_folder'].creator
            
            #self.read_results(wchain=self.calculation)
