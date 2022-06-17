
from abc import ABC, abstractmethod
from xmlrpc.client import boolean

from koopmans import calculators
from ._workflow import Workflow
from koopmans import utils

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

    def __init__(self, ml_model:MLModel, orbital:int, calc_that_produced_orbital_densities, filled:boolean=True, *args, **kwargs):
        self.n_max                                = 4
        self.l_max                                = 4
        self.r_min                                = 0.5
        self.r_max                                = 4.0
        self.method_to_extract_from_binary        = 'from_ki'
        self.orbital                              = orbital #=-1 corresponds to the total density
        self.filled                               = filled
        self.calc_that_produced_orbital_densities = calc_that_produced_orbital_densities
        
        super().__init__(**kwargs, ml_model=ml_model)

    def _run(self):
        self.extract_input_vector_for_ML_model()


    def extract_input_vector_for_ML_model(self):
        self.convert_binary_to_xml()
        if self.orbital != -1:
            self.compute_decomposition()
            self.compute_power_spectrum()

    def convert_binary_to_xml(self):
        print("Convert binary to xml")
        orbital_densities_bin_dir            = str(self.calc_that_produced_orbital_densities.parameters.outdir) + '/kc_' + str(self.calc_that_produced_orbital_densities.parameters.ndw) + '.save'
        orbital_density_xml_dir              = str(self.calc_that_produced_orbital_densities.directory) + '/ML/TMP'

        if self.method_to_extract_from_binary == 'from_ki':
            utils.system_call(f'mkdir -p {orbital_density_xml_dir}')
            if self.orbital==-1:
                # extract total density
                filename = '/charge-density'
            else:
                # extract orbital density
                if self.filled:
                    occupation = 'occ'
                else:
                    occupation  = 'empty'
                filename     = '/sic_potential.' + occupation  + '.' + str(self.orbital) 
            filepath_bin = orbital_densities_bin_dir + filename + '.dat'
            filepath_xml = orbital_density_xml_dir   + filename + '.xml'
            command  = str(calculators.bin_directory) + '/bin2xml_real_space_density.x ' + filepath_bin + ' ' + filepath_xml
            utils.system_call(command)

    def compute_decomposition(self):
        print("compute decomposition")
    
    def compute_power_spectrum(self):
        print("compute power spectrum")

    def predict(self):
        self.ml_model.predict()
    
    def train(self):
        self.ml_model.train()
    
    def get_prediction_for_latest_alpha(self):
        return 0.65
    
    def use_prediction(self):
        print("Use prediction -> False")
        return False