import copy
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Optional, Union

import json
import numpy as np
from deepdiff import DeepDiff
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

class AbstractPredictor(ABC):

    def __init__(self) -> None:
        self.serialize_as_list: Dict = {}
        self.is_trained: bool = False

    @abstractmethod
    def predict(self, x_test: np.ndarray) -> np.ndarray:
        ...

    @abstractmethod
    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        ...

    def todict(self) -> Dict:
        dct = copy.deepcopy(dict(self.__dict__))
        serialize_as_list = dct.pop('serialize_as_list')
        if serialize_as_list != {}:
            ser_dct: Dict = {}
            for item in serialize_as_list.keys():
                to_serialize = dct.pop(item)
                ser_dct[item] = {}
                ser_dct[item]['init_params'] = to_serialize.get_params()
                ser_dct[item]['model_params'] = mp = {}
                for p in self.serialize_as_list[item]['list']:
                    if hasattr(to_serialize, p):
                        mp[p] = getattr(to_serialize,p).tolist()
                dct['serialize_as_list'] = ser_dct
        return dct
    
    @classmethod
    def fromdict(cls, dct: Dict):
        predictor = cls()
        if dct['serialize_as_list'] is not None:
            for item in dct['serialize_as_list'].keys():
                model = predictor.serialize_as_list[item]['class']
                model.set_params(**dct['serialize_as_list'][item]['init_params'])
                for name, p in dct['serialize_as_list'][item]['model_params'].items():
                    setattr(model, name, np.asarray(p))
        predictor.is_trained = dct['is_trained']
        return predictor

    def save_to_file(self, save_dir: Path) -> None:
        dct = self.todict()
        with open(save_dir / f'model.save', "w") as fp:
            json.dump(dct, fp)

    @classmethod
    def load_from_file(cls, save_dir: Path):
        with open(save_dir / f'model.save', "r") as fp:
            dct = json.load(fp)
        return cls.fromdict(dct)
    

class RidgeRegressionModel(AbstractPredictor):
    def __init__(self) -> None:
        self.scaler = StandardScaler()
        self.model = Ridge(alpha=1.0)
        self.is_trained = False
        self.name = 'ridge_regression'
        self.serialize_as_list = {'model': {'class': self.model, 'list': ['coef_', 'intercept_']},\
                                  'scaler': {'class': self.scaler, 'list': ['mean_', 'scale_']}}

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
    def __init__(self) -> None:
        self.model = Ridge(alpha=0.0)
        self.is_trained = False
        self.name = 'linear_regression'
        self.serialize_as_list = {'model': ['coef_', 'intercept_']}

    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        self.model.fit(X_train, Y_train)
        self.is_trained = True

    def predict(self, x_test: np.ndarray) -> np.ndarray:
        X_test = np.atleast_2d(x_test)
        y_predict = self.model.predict(X_test)
        return y_predict


class MeanModel(AbstractPredictor):
    def __init__(self) -> None:
        self.mean = 0.0
        self.is_trained = True
        self.name = 'mean'
        self.serialize_as_list = {}

    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        self.mean = np.mean(Y_train)

    def predict(self, x_test: np.ndarray) -> np.ndarray:
        shape = np.shape(x_test)[0]
        return self.mean*np.ones(shape)


class MLModel():
    def __init__(self, type_of_ml_model='ridge_regression',
                 X_train: Optional[np.ndarray] = None, Y_train: Optional[np.ndarray] = None, 
                 save_dir: Optional[Path]=None):
        self.X_train = X_train
        self.Y_train = Y_train
        self.type_of_ml_model = type_of_ml_model
        self.model: AbstractPredictor = self.init_and_reset_model(save_dir)

    def init_and_reset_model(self, save_dir: Optional[Path]=None) -> AbstractPredictor:
        model: AbstractPredictor
        if self.type_of_ml_model == 'ridge_regression':
            if save_dir is None:
                model = RidgeRegressionModel()
            else: 
                model = RidgeRegressionModel.load_from_file(save_dir)
        elif self.type_of_ml_model == 'linear_regression':
            if save_dir is None:
                model = LinearRegressionModel()
            else:
                model = LinearRegressionModel.load_from_file(save_dir)
        elif self.type_of_ml_model == 'mean':
            if save_dir is None:
                model = MeanModel()
            else:
                model = MeanModel.load_from_file(save_dir)
        else:
            raise ValueError(f"{self.type_of_ml_model} is not implemented as a valid ML model.")
        return model

    def __repr__(self):
        if self.X_train and self.Y_train:
            return f'{self.type_of_ml_model}(is_trained={self.model.is_trained}, ' \
                   f'number_of_training_vectors={self.X_train.shape[0]}, ' \
                   f'input_vector_dimension={self.Y_train.shape[0]})'
        else:
            return f'{self.type_of_ml_model}(is_trained={self.model.is_trained}, ' \
                   'no training data has been added so far)'

    def todict(self):
        dct = copy.deepcopy(dict(self.__dict__))
        dct['model'] = self.model.todict()
        return dct
    
    @classmethod
    def fromdict(cls, dct: Dict):
        predictor_dct = dct.pop('model')
        ml_model = cls(**dct)
        ml_model.model = ml_model.model.fromdict(predictor_dct)
        return ml_model

    def predict(self, x_test: np.ndarray):
        """
        Make a prediction of using the trained model.
        """

        if self.model.is_trained:
            y_predict = self.model.predict(x_test)
            return y_predict
        else:
            return np.array([1.0])  # dummy value

    def train(self):
        """
        Reset the model and train the model (including the StandardScaler) with all training data added so far.
        """
        self.init_and_reset_model()
        self.model.fit(self.X_train, self.Y_train)

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
        items_to_pop = ['model']
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
