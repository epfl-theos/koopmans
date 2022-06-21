
from abc import ABC, abstractmethod
from xmlrpc.client import boolean

from koopmans import calculators
from ._workflow import Workflow
from koopmans import utils

from koopmans import ML_utils
import numpy as np

from typing import List, Dict, Any

import copy

# class MLModel(ABC):

#     @abstractmethod
#     def predict(self) -> None:
#         ...
    
#     @abstractmethod
#     def train(self) -> None:
#         ...

class MLModel():

    def predict(self):
        print("model.predict()")

    def train(self):
        print("model.train()")


class MLCapableWorkflow(Workflow):

    def __init__(self, *args, ml_model:MLModel=None, **kwargs):
        if ml_model is not None:
            self.ml_model = ml_model
        
        super().__init__(*args, **kwargs)
        
    
    # @abstractmethod
    # def convert_binary_to_xml(self) -> None:
    #     ...
    
    @property
    def wf_kwargs(self) -> Dict[str, Any]:
        h = super().wf_kwargs
        h.update({'ml_model': copy.deepcopy(self.ml_model)})
        return h


class MLFiitingWorkflow(MLCapableWorkflow):

    def __init__(self, calc_that_produced_orbital_densities, **kwargs):
        super().__init__(**kwargs)
        self.n_max                                = 4
        self.l_max                                = 4
        self.r_min                                = 0.5
        self.r_max                                = 4.0
        self.method_to_extract_from_binary        = 'from_ki'
        # TODO: don't hard-code the number of bands
        self.num_occ_bands                        = 4
        self.num_emp_bands                        = 2
        # end TODO
        self.calc_that_produced_orbital_densities = calc_that_produced_orbital_densities
        self.ML_dir                               = self.calc_that_produced_orbital_densities.directory / 'ML' / 'TMP'

        import ipdb 
        ipdb.set_trace()
        
        

    def _run(self):
        self.extract_input_vector_for_ML_model()


    def extract_input_vector_for_ML_model(self):
        self.convert_binary_to_xml()
        self.compute_decomposition()
        self.compute_power_spectrum()
            

    def convert_binary_to_xml(self):
        print("Convert binary to xml")
        orbital_densities_bin_dir            = self.calc_that_produced_orbital_densities.parameters.outdir/ f'kc_{self.calc_that_produced_orbital_densities.parameters.ndw}.save'
        

        if self.method_to_extract_from_binary == 'from_ki':
            utils.system_call(f'mkdir -p {self.ML_dir}')
            command  = str(calculators.bin_directory / 'bin2xml_real_space_density.x ') + ' '.join(str(x) for x in [orbital_densities_bin_dir, self.ML_dir, self.num_occ_bands, self.num_emp_bands])
            utils.system_call(command)

    def compute_decomposition(self):
        print("compute decomposition")
        self.r_cut = min(self.atoms.get_cell_lengths_and_angles()[:3])/2.5 #the maximum radius has to be smaller than half of the cell-size
        print(self.r_cut)
        if self.method_to_extract_from_binary == 'from_ki':
            centers_occ = np.array(self.calculations[4].results['centers'])
            centers_emp = np.array(self.calculations[7].results['centers'])
        ML_utils.precompute_radial_basis(self.n_max, self.l_max, self.r_min, self.r_max, self.ML_dir)
        ML_utils.func_compute_decomposition(self.n_max, self.l_max, self.r_min, self.r_max, self.r_cut, self.ML_dir, [self.num_emp_bands, self.num_occ_bands], self.atoms, centers_occ, centers_emp)
    
    def compute_power_spectrum(self):
        print("compute power spectrum")
        self.dir_power = self.ML_dir / f'power_spectra_{self.n_max}_{self.l_max}_{self.r_min}_{self.r_max}'
        ML_utils.main_compute_power(self.n_max, self.l_max, self.r_min, self.r_max, self.ML_dir, self.dir_power, [self.num_emp_bands, self.num_occ_bands])

    def predict(self):
        
        self.ml_model.predict()
    
    def train(self, new_orbitals:List[int]=None, new_orbitals_filling:List[bool]=None, new_alphas:List[float]=None):
        print("Now training")
        if new_orbitals != None:
            if not(len(new_orbitals) == len(new_orbitals_filling) == len(new_alphas)):
                raise ValueError('new_orbitals and new_alphas have to be of the same length')
            self.add_training_data(new_orbitals, new_orbitals_filling, new_alphas)
        self.ml_model.train()
    
    def add_training_data(self, new_orbitals:List[int], new_orbitals_filling:List[bool], new_alphas:List[float]):
        for orbital, orbital_filling, alpha in zip(new_orbitals, new_orbitals_filling, new_alphas):
            if orbital_filling:
                occ_string = 'occ'
            else:
                occ_string = 'emp'
            power_spectrum = np.loadtxt(self.dir_power / f'power_spectrum.orbital.{occ_string}.{orbital}.txt')
            print(power_spectrum)
            print("adding orbital ", orbital)
            print("filling orbital ", orbital_filling)
            print("alpha orbital ", alpha)

    
    def get_prediction_for_latest_alpha(self):
        return 0.65
    
    def use_prediction(self):
        print("Use prediction -> False")
        return False