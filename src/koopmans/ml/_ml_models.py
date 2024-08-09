import copy
from abc import ABC, abstractmethod
from typing import Callable, Dict, List, Optional, Union

import numpy as np
from deepdiff import DeepDiff
from sklearn.base import BaseEstimator
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from koopmans.bands import Band
from koopmans.utils import get_binary_content


class AbstractPredictor(ABC):

    def __init__(self, model, is_trained=False) -> None:
        self.model = model
        self.is_trained = is_trained

    @abstractmethod
    def predict(self, x_test: np.ndarray) -> np.ndarray:
        ...

    @abstractmethod
    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        ...

    def todict(self):
        dct = dict(self.__dict__)
        dct['__koopmans_module__'] = self.__module__
        dct['__koopmans_name__'] = self.__class__.__name__
        return dct

    @classmethod
    def fromdict(cls, dct: Dict):
        return cls(**dct)


class RidgeRegressionModel(AbstractPredictor):
    def __init__(self, model=None, scaler=None, is_trained=False) -> None:
        model = Ridge(alpha=1.0) if model is None else model
        super().__init__(model, is_trained)
        self.scaler = StandardScaler() if scaler is None else scaler

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
    def __init__(self, model=None, is_trained=False) -> None:
        model = Ridge(alpha=0.0) if model is None else model
        super().__init__(model, is_trained)

    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        self.model.fit(X_train, Y_train)
        self.is_trained = True

    def predict(self, x_test: np.ndarray) -> np.ndarray:
        X_test = np.atleast_2d(x_test)
        y_predict = self.model.predict(X_test)
        return y_predict


class MeanModel(AbstractPredictor):
    def __init__(self, mean=0.0, is_trained=True) -> None:
        self.mean = mean
        self.is_trained = is_trained

    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        self.mean = np.mean(Y_train)

    def predict(self, x_test: np.ndarray) -> np.ndarray:
        shape = np.shape(x_test)[0]
        return self.mean*np.ones(shape)


def model_factory(type_of_ml_model: str) -> AbstractPredictor:
    if type_of_ml_model == 'ridge_regression':
        return RidgeRegressionModel()
    elif type_of_ml_model == 'linear_regression':
        return LinearRegressionModel()
    elif type_of_ml_model == 'mean':
        return MeanModel()
    else:
        raise ValueError(f"{type_of_ml_model} is not implemented as a valid ML model.")


class AbstractMLModel(ABC):

    def __init__(self, type_of_ml_model: str = 'ridge_regression', descriptor='orbital_density', is_trained: bool = False,
                 model: AbstractPredictor | None = None, descriptor_from_band: Callable[[Band], np.ndarray] | None = None):
        self.type_of_ml_model = type_of_ml_model
        self.model = model_factory(type_of_ml_model) if model is None else model
        self.model.is_trained = is_trained
        self.descriptor_from_band = descriptor_from_band_factory(
            descriptor) if descriptor_from_band is None else descriptor_from_band

    @abstractmethod
    def add_training_data(self, bands: List[Band]) -> None:
        ...

    @abstractmethod
    def train(self) -> None:
        ...

    @property
    @abstractmethod
    def is_trained(self) -> bool:
        ...

    @abstractmethod
    def predict(self, band: Band) -> float:
        ...

    def __eq__(self, other):
        items_to_pop = ['model', 'model_occ', 'model_emp']
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

    _fields_to_exclude: List[str] = []

    def todict(self):
        dct = copy.deepcopy(dict(self.__dict__))
        dct['__koopmans_module__'] = self.__module__
        dct['__koopmans_name__'] = self.__class__.__name__
        return dct

    @classmethod
    def fromdict(cls, dct: Dict):
        return cls(**dct)

    @abstractmethod
    def _repr_fields(self) -> str:
        ...

    def __repr__(self):
        return f'{self.__class__.__name__}({self._repr_fields()})'


def power_spectrum_from_band(band: Band) -> np.ndarray:
    assert band.power_spectrum is not None
    return np.frombuffer(get_binary_content(*band.power_spectrum))


def self_hartree_from_band(band: Band) -> np.ndarray:
    return np.array([band.self_hartree])


def descriptor_from_band_factory(descriptor: str) -> Callable[[Band], np.ndarray]:
    if descriptor == 'orbital_density':
        return power_spectrum_from_band
    elif descriptor == 'self_hartree':
        return self_hartree_from_band
    else:
        raise ValueError(f"{descriptor} is not implemented as a valid descriptor.")


class MLModel(AbstractMLModel):

    _fields_to_exclude = ['model']

    def __init__(self, type_of_ml_model='ridge_regression', descriptor: str = 'orbital_density', is_trained: bool = False,
                 X_train: Optional[np.ndarray] = None, Y_train: Optional[np.ndarray] = None, **kwargs):
        super().__init__(type_of_ml_model=type_of_ml_model, descriptor=descriptor, is_trained=is_trained, **kwargs)
        self.X_train = X_train
        self.Y_train = Y_train

    @property
    def is_trained(self):
        return self.model.is_trained

    def _repr_fields(self):
        out = f"model={self.model.__class__.__name__}(...), is_trained={self.is_trained}"
        if self.X_train is not None:
            out += f", len(X_train)={len(self.X_train)}"
        if self.Y_train is not None:
            out += f", shape(Y_train)={self.Y_train.shape}"
        return out

    def predict(self, band: Band) -> float:
        """
        Make a prediction of using the trained model.
        """

        if not self.is_trained:
            raise ValueError(f'{self.__class__.__name__} must be trained before calling predict()')

        descriptor = self.descriptor_from_band(band)

        return self.model.predict(descriptor)[0]

    def train(self):
        """
        Reset the model and train the model (including the StandardScaler) with all training data added so far.
        """
        self.model = model_factory(self.type_of_ml_model)
        self.model.fit(self.X_train, self.Y_train)

    def add_training_data(self, bands: List[Band]):  # x_train: np.ndarray, y_train: Union[float, np.ndarray]):
        """
        Add training data to the model.
        """

        x_train = np.atleast_2d([self.descriptor_from_band(b) for b in bands])
        y_train = np.atleast_1d([b.alpha for b in bands])

        if self.X_train is None or self.Y_train is None:
            self.X_train = x_train
            self.Y_train = y_train
        else:
            self.X_train = np.concatenate([self.X_train, x_train])
            self.Y_train = np.concatenate([self.Y_train, y_train])


class OccEmpMLModels(AbstractMLModel):
    """A model that contains two sub-models, one for occupied and one for empty."""

    _fields_to_exclude = ['model_occ', 'model_emp']

    def __init__(self, type_of_ml_model='ridge_regression', descriptor='orbital_from_band', is_trained: bool = False,
                 X_train_occ: Optional[np.ndarray] = None, Y_train_occ: Optional[np.ndarray] = None,
                 X_train_emp: Optional[np.ndarray] = None, Y_train_emp: Optional[np.ndarray] = None,
                 model_occ=None, model_emp=None, descriptor_from_band: Callable[[Band], np.ndarray] | None = None):

        self.model_occ = MLModel(type_of_ml_model=type_of_ml_model, descriptor=descriptor, is_trained=is_trained,
                                 X_train=X_train_occ, Y_train=Y_train_occ, model=model_occ, descriptor_from_band=descriptor_from_band)
        self.model_emp = MLModel(type_of_ml_model=type_of_ml_model, descriptor=descriptor, is_trained=is_trained,
                                 X_train=X_train_emp, Y_train=Y_train_emp, model=model_emp, descriptor_from_band=descriptor_from_band)
        self.type_of_ml_model = type_of_ml_model

    def add_training_data(self, bands: List[Band]) -> None:
        for band in bands:
            if band.filled:
                self.model_occ.add_training_data([band])
            else:
                self.model_emp.add_training_data([band])

    def train(self) -> None:
        self.model_occ.train()
        self.model_emp.train()

    @property
    def is_trained(self):
        return self.model_occ.is_trained and self.model_emp.is_trained

    def predict(self, band: Band) -> np.ndarray:
        if band.filled:
            return self.model_occ.predict(band)
        else:
            return self.model_emp.predict(band)

    def _repr_fields(self):
        out = f"model='{self.model_occ.__class__.__name__}(...)', is_trained={self.is_trained}"
        if self.model_occ.X_train is not None:
            out += f", len(X_train, occ)={len(self.model_occ.X_train)}"
        if self.model_occ.Y_train is not None:
            out += f", shape(Y_train, occ)={self.model_occ.Y_train.shape}"
        if self.model_emp.X_train is not None:
            out += f", len(X_train, emp)={len(self.model_emp.X_train)}"
        if self.model_emp.Y_train is not None:
            out += f", shape(Y_train, emp)={self.model_emp.Y_train.shape}"
        return out
