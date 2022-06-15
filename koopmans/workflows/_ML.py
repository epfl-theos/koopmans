
from abc import ABC, abstractmethod
from ._workflow import Workflow

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

    def __init__(self, *args, ml_model:MLModel=None, **kwargs):
        self.n_max = 4
        self.l_max = 4
        self.r_min = 0.5
        self.r_max = 4.0
        super().__init__(**kwargs, ml_model=ml_model)

    def _run(self):
        self.extract_input_vector_for_ML_model()


    def extract_input_vector_for_ML_model(self):
        self.convert_binary_to_xml()
        self.compute_decomposition()
        self.compute_power_spectrum()

    def convert_binary_to_xml(self):
        print("Convert binary to xml")

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
        if(False):
            return False