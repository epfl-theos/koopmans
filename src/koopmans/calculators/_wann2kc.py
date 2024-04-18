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
    
    control_namelist = [
                    'kcw_iverbosity',
                    'kcw_at_ks'      ,
                    'lrpa'           ,
                    'mp1'            ,
                    'mp2'            ,
                    'mp3'            ,
                    'homo_only'      ,
                    'read_unitary_matrix' ,
                    'l_vcut'         ,
                    'spin_component' ,
                    ]
    
    control_dict = {k:v if k in control_namelist else None for k,v in wann2kc_calculator.parameters.items()}
    control_dict['calculation'] = "wann2kcw"

    if not any(wann2kc_calculator.atoms.pbc): control_dict['assume_isolated'] = "m-t"
    
    for k in list(control_dict):
        if control_dict[k] is None:
            control_dict.pop(k)

    wannier_namelist = [
                    "check_ks"       ,
                    "num_wann_occ"   ,
                    "num_wann_emp"   ,
                    "have_empty"     ,
                    "has_disentangle",
    ]

    wannier_dict = {k:v if k in wannier_namelist else None for k,v in wann2kc_calculator.parameters.items()}

    for k in list(wannier_dict):
        if wannier_dict[k] is None:
            wannier_dict.pop(k)
            
    wann2kcw_params = {
            "CONTROL":control_dict,
            "WANNIER":wannier_dict,
        }
    
    
    builder.parameters = orm.Dict(wann2kcw_params)
    builder.code = orm.load_code(wann2kc_calculator.mode["kcw_code"])
    builder.metadata = wann2kc_calculator.mode["metadata"]
    if "metadata_kcw" in wann2kc_calculator.mode: builder.metadata = wann2kc_calculator.mode["metadata_kcw"]
    builder.parent_folder = wann2kc_calculator.parent_folder
    
    if hasattr(wann2kc_calculator, "wannier90_files"):
        builder.wann_u_mat = wann2kc_calculator.wannier90_files["occ"]["u_mat"]
        builder.wann_emp_u_mat = wann2kc_calculator.wannier90_files["emp"]["u_mat"]
        builder.wann_emp_u_dis_mat = wann2kc_calculator.wannier90_files["emp"]["u_dis_mat"]
        builder.wann_centres_xyz = wann2kc_calculator.wannier90_files["occ"]["centres_xyz"]
        builder.wann_emp_centres_xyz = wann2kc_calculator.wannier90_files["emp"]["centres_xyz"]
        
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
            from aiida.engine import run_get_node,submit
            running = run_get_node(builder)
            
            # once the running if completed
            self.calculation = running[-1]
            
            #self.read_results(wchain=self.calculation)
