
from abc import ABC, abstractmethod
from xmlrpc.client import boolean

from koopmans import calculators
from ._workflow import Workflow
from koopmans import utils

from koopmans import ML_utils
import numpy as np

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
        
    
    @abstractmethod
    def convert_binary_to_xml(self) -> None:
        ...
    
    


class MLFiitingWorkflow(MLCapableWorkflow):

    def __init__(self, ml_model:MLModel, orbital:int, calc_that_produced_orbital_densities, workflow:Workflow, filled:boolean=True, *args, **kwargs):
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
        self.ML_dir               = str(self.calc_that_produced_orbital_densities.directory) + '/ML/TMP'
        self.whole_workflow                        = workflow
        
        super().__init__(**kwargs, ml_model=ml_model)

    def _run(self):
        self.extract_input_vector_for_ML_model()


    def extract_input_vector_for_ML_model(self):
        self.convert_binary_to_xml()
        self.compute_decomposition()
        self.compute_power_spectrum()
            

    def convert_binary_to_xml(self):
        print("Convert binary to xml")
        orbital_densities_bin_dir            = str(self.calc_that_produced_orbital_densities.parameters.outdir) + '/kc_' + str(self.calc_that_produced_orbital_densities.parameters.ndw) + '.save'
        

        if self.method_to_extract_from_binary == 'from_ki':
            utils.system_call(f'mkdir -p {self.ML_dir}')
            command  = str(calculators.bin_directory) + '/bin2xml_real_space_density.x ' + orbital_densities_bin_dir + ' ' + self.ML_dir + ' ' + str(self.num_occ_bands) + ' ' + str(self.num_emp_bands)
            utils.system_call(command)

    def compute_decomposition(self):
        print("compute decomposition")
        self.r_cut = np.min(self.whole_workflow.atoms.get_cell())/2.0
        if self.method_to_extract_from_binary == 'from_ki':
            centers_occ = np.array(self.whole_workflow.calculations[4].results['centers'])
            centers_emp = np.array(self.whole_workflow.calculations[7].results['centers'])
        ML_utils.precompute_radial_basis(self.n_max, self.l_max, self.r_min, self.r_max, self.ML_dir)
        ML_utils.func_compute_decomposition(self.n_max, self.l_max, self.r_min, self.r_max, self.r_cut, self.ML_dir, [self.num_emp_bands, self.num_occ_bands], self.whole_workflow, centers_occ, centers_emp)
    
    def compute_power_spectrum(self):
        print("compute power spectrum")
        ML_utils.main_compute_power(self.n_max, self.l_max, self.ML_dir, [self.num_emp_bands, self.num_occ_bands])

    def predict(self):
        self.ml_model.predict()
    
    def train(self):
        self.ml_model.train()
    
    def get_prediction_for_latest_alpha(self):
        return 0.65
    
    def use_prediction(self):
        print("Use prediction -> False")
        return False