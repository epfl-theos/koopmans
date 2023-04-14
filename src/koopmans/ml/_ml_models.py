import copy
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np
from deepdiff import DeepDiff
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

class AbstractPredictor(ABC):

    @abstractmethod
    def predict(self, x_test: np.ndarray) -> np.ndarray:
        ...

    @abstractmethod
    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        ...

    def todict(self):
        dct = dict(self.__dict__)
        return dct

    @classmethod
    def fromdict(cls, dct):
        return cls(**dct)

    def save_to_file(self, save_dir: Path) -> None:
        from koopmans.io import write
        write(self, save_dir / f'model.kwf')

    @classmethod
    def load_from_file(cls, save_dir: Path):
        from koopmans.io import read
        dct = read(save_dir / f'model.kwf')
        return cls.fromdict(dct)
    

class RidgeRegressionModel(AbstractPredictor):
    def __init__(self, scaler = StandardScaler(), model = Ridge(alpha=1.0), is_trained=False, name='ridge_regression') -> None:
        self.scaler = scaler
        self.model = model
        self.is_trained = is_trained
        self.name = name

    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        self.scaler = self.scaler.fit(X_train)
        X_train_scaled = self.scaler.transform(X_train)
        self.model.fit(X_train_scaled, Y_train)
        self.is_trained = True

    def predict(self, x_test: np.ndarray) -> np.ndarray:
        X_test = np.atleast_2d(x_test)
        X_test = self.scaler.transform(X_test)
        y_predict = self.model.predict(X_test)
        return y_predict
    

class LinearRegressionModel(AbstractPredictor):
    def __init__(self, model = Ridge(alpha=0.0), is_trained = False, name = 'linear_regression') -> None:
        self.model = model
        self.is_trained = is_trained
        self.name = name

    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        self.model.fit(X_train, Y_train)
        self.is_trained = True

    def predict(self, x_test: np.ndarray) -> np.ndarray:
        X_test = np.atleast_2d(x_test)
        y_predict = self.model.predict(X_test)
        return y_predict


class MeanModel(AbstractPredictor):
    def __init__(self, mean = 0.0, is_trained = True, name = 'mean') -> None:
        self.mean = mean
        self.is_trained = is_trained
        self.name = name

    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        self.mean = np.mean(Y_train)

    def predict(self, x_test: np.ndarray) -> np.ndarray:
        shape = np.shape(x_test)[0]
        return self.mean*np.ones(shape)


class MLModel():
    def __init__(self, type_of_ml_model: str ='ridge_regression',
                 X_train: Optional[np.ndarray] = None, Y_train: Optional[np.ndarray] = None, 
                 save_dir: Optional[Path]=None, sub_dir: Optional[str] = None):
        self.X_train = X_train
        self.Y_train = Y_train
        self.type_of_ml_model = type_of_ml_model
        self.model_class: AbstractPredictor = self.init_and_reset_model(save_dir, sub_dir)

    def init_and_reset_model(self, save_dir: Optional[Path]=None, sub_dir: Optional[str] = None) -> AbstractPredictor:
        model_class: AbstractPredictor
        model_classes = {'ridge_regression': RidgeRegressionModel, 'linear_regression': LinearRegressionModel, 'mean': MeanModel}
        cls = model_classes[self.type_of_ml_model]
        if save_dir is None:
            model_class = cls()
        else:
            if sub_dir is None:
                model_class = cls.load_from_file(save_dir)
            else:
                model_class = cls.load_from_file(save_dir / sub_dir)
        return model_class

    def __repr__(self):
        if self.X_train and self.Y_train:
            return f'{self.type_of_ml_model}(is_trained={self.model_class.is_trained}, ' \
                   f'number_of_training_vectors={self.X_train.shape[0]}, ' \
                   f'input_vector_dimension={self.Y_train.shape[0]})'
        else:
            return f'{self.type_of_ml_model}(is_trained={self.model_class.is_trained}, ' \
                   'no training data has been added so far)'

    def todict(self):
        # Shallow copy
        dct = dict(self.__dict__)
        return dct
    
    @classmethod
    def fromdict(cls, dct: Dict):
        predictor_class = dct.pop('model_class')
        ml_model = cls(**dct)
        ml_model.model_class = predictor_class
        return ml_model

    def predict(self, x_test: np.ndarray):
        """
        Make a prediction of using the trained model.
        """

        if self.model_class.is_trained:
            y_predict = self.model_class.predict(x_test)
            return y_predict
        else:
            return np.array([1.0])  # dummy value

    def train(self):
        """
        Reset the model and train the model (including the StandardScaler) with all training data added so far.
        """
        self.init_and_reset_model()
        self.model_class.fit(self.X_train, self.Y_train)

    def add_training_data(self, x_train: np.ndarray, y_train: Union[float, np.ndarray]):
        """
        Add training data to the model.
        """

        x_train = np.atleast_2d(x_train)
        y_train = np.atleast_1d(y_train)

        if self.X_train is None or self.Y_train is None:
            self.X_train = x_train
            self.Y_train = y_train
        else:
            self.X_train = np.concatenate([self.X_train, x_train])
            self.Y_train = np.concatenate([self.Y_train, y_train])

    def __eq__(self, other):
        items_to_pop = ['model_class']
        if isinstance(other, MLModel):
            self_dict = copy.deepcopy(self.__dict__)
            other_dict = copy.deepcopy(other.__dict__)
            for item in items_to_pop:
                if item in other_dict:
                    other_dict.pop(item)
                if item in self_dict:
                    self_dict.pop(item)
            return DeepDiff(self_dict, other_dict, significant_digits=8, number_format_notation='e') == {}
        else:
            return False
