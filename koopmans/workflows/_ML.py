
from abc import ABC, abstractmethod

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


class MLCapableWorkflow(ABC):

    def __init__(self, *args, **kwargs):
        return
    
    @abstractmethod
    def compute_real_space_density(self) -> None:
        ...
    


class MLFiitingWorkflow():
    def __init__(self, *args, ml_model:MLModel, **kwargs):
        super().__init__(*args, **kwargs)
        self.ml_model = ml_model

    def convert_binary_to_xml(self):
        print("Convert binary to xml")

    def compute_decomposition(self):
        print("compute decomposition")
    
    def compute_power_spectrum(self):
        print("compute power spectrum")
    
    def extract_input_vector_for_ML_model(self):
        self.convert_binary_to_xml()
        self.compute_decomposition()
        self.compute_power_spectrum()

    def run(self):
        self.extract_input_vector_for_ML_model()
    
    def predict(self):
        self.ml_model.predict()
    
    def train(self):
        self.ml_model.train()