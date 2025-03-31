import copy
from abc import ABC, abstractmethod
from functools import partial
from typing import Callable, Dict, List, Optional

import numpy as np
from deepdiff import DeepDiff
from sklearn.base import BaseEstimator
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from koopmans.bands import Band
from koopmans.engines import Engine


class WrappedEstimator(ABC):
    """Abstract class for wrapping a scikit-learn estimator."""

    def __init__(self, estimator: BaseEstimator, is_trained: bool = False) -> None:
        self.estimator = estimator
        self.is_trained = is_trained

    @abstractmethod
    def predict(self, x_test: np.ndarray) -> np.ndarray:
        """Make a prediction for the provided input."""
        ...

    @abstractmethod
    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        """Fit the estimator to the provided training data."""
        ...

    def todict(self):
        """Convert the class to a dictionary."""
        dct = dict(self.__dict__)
        dct['__koopmans_module__'] = self.__module__
        dct['__koopmans_name__'] = self.__class__.__name__
        return dct

    @classmethod
    def fromdict(cls, dct: Dict):
        """Construct an instance of this class from a dictionary."""
        return cls(**dct)


class RidgeRegressionEstimator(WrappedEstimator):
    """An estimator that uses ridge regression."""

    def __init__(self, estimator=None, scaler=None, is_trained=False) -> None:
        estimator = Ridge(alpha=1.0) if estimator is None else estimator
        super().__init__(estimator, is_trained)
        self.scaler = StandardScaler() if scaler is None else scaler

    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        """Fit the estimator with appropriate scaling."""
        self.scaler = self.scaler.fit(X_train)
        X_train_scaled = self.scaler.transform(X_train)
        self.estimator.fit(X_train_scaled, Y_train)
        self.is_trained = True

    def predict(self, x_test: np.ndarray) -> np.ndarray:
        """Rescale and fit the data."""
        X_test = np.atleast_2d(x_test)
        X_test = self.scaler.transform(X_test)
        y_predict = self.estimator.predict(X_test)
        return y_predict

    def todict(self):
        """Convert the class to a dictionary."""
        dct = {k: v for k, v in self.__dict__.items() if not isinstance(v, BaseEstimator)}
        dct['__koopmans_module__'] = self.__module__
        dct['__koopmans_name__'] = self.__class__.__name__
        return dct


class LinearRegressionEstimator(WrappedEstimator):
    """An estimator that uses linear regression."""

    def __init__(self, estimator=None, is_trained=False) -> None:
        estimator = Ridge(alpha=0.0) if estimator is None else estimator
        super().__init__(estimator, is_trained)

    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        """Fit the estimator to the training data."""
        self.estimator.fit(X_train, Y_train)
        self.is_trained = True

    def predict(self, x_test: np.ndarray) -> np.ndarray:
        """Make a prediction."""
        X_test = np.atleast_2d(x_test)
        y_predict = self.estimator.predict(X_test)
        return y_predict


class MeanEstimator(WrappedEstimator):
    """An estimator that simply uses the mean of the training data."""

    def __init__(self, estimator=None, is_trained=False) -> None:
        estimator = DummyRegressor(strategy='mean') if estimator is None else estimator
        super().__init__(estimator, is_trained)

    def fit(self, X_train: np.ndarray, Y_train: np.ndarray):
        """Fit the estimator to the mean of the training data."""
        self.estimator.fit(X_train, Y_train)
        self.is_trained = True

    def predict(self, x_test: np.ndarray) -> np.ndarray:
        """Return the mean of the training data."""
        X_test = np.atleast_2d(x_test)
        y_predict = self.estimator.predict(X_test)
        return y_predict


def estimator_factory(estimator: str) -> WrappedEstimator:
    """Create an estimator from a string."""
    if estimator == 'ridge_regression':
        return RidgeRegressionEstimator()
    elif estimator == 'linear_regression':
        return LinearRegressionEstimator()
    elif estimator == 'mean':
        return MeanEstimator()
    else:
        raise ValueError(f"`{estimator}` is not implemented as a valid ML estimator.")


class AbstractMLModel(ABC):
    """An abstract class for ML models that predict the screening parameters of Bands."""

    def __init__(self, estimator_type: str = 'ridge_regression', descriptor='orbital_density', is_trained: bool = False,
                 estimator: WrappedEstimator | None = None,
                 descriptor_from_band: Callable[[Band], np.ndarray] | None = None, engine: Engine | None = None):
        self.estimator_type = estimator_type
        self.estimator = estimator_factory(estimator_type) if estimator is None else estimator
        self.estimator.is_trained = is_trained

        if isinstance(self.estimator, MeanEstimator):
            # MeanEstimator does not require a descriptor
            self.descriptor_from_band = dummy_from_band if descriptor_from_band is None else descriptor_from_band
        else:
            if descriptor_from_band is None:
                assert engine is not None
                self.descriptor_from_band = descriptor_from_band_factory(descriptor, engine)
            else:
                self.descriptor_from_band = descriptor_from_band
        self.descriptor_type = descriptor

    @abstractmethod
    def add_training_data(self, bands: List[Band]) -> None:
        """Add training data to the estimator."""
        ...

    @abstractmethod
    def train(self) -> None:
        """Train the estimator."""
        ...

    @property
    @abstractmethod
    def is_trained(self) -> bool:
        """Return whether the estimator is trained."""
        ...

    @abstractmethod
    def predict(self, band: Band) -> float:
        """Predict the screening parameter of the provided Band."""
        ...

    def __eq__(self, other):
        items_to_pop = ['estimator']
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

    def todict(self):
        """Convert the class to a dictionary."""
        dct = {k: v for k, v in self.__dict__.items()}
        dct['__koopmans_module__'] = self.__module__
        dct['__koopmans_name__'] = self.__class__.__name__
        return dct

    @classmethod
    def fromdict(cls, dct: Dict):
        """Create an instance of the class from a dictionary."""
        return cls(**dct)

    @abstractmethod
    def _repr_fields(self) -> str:
        ...

    def __repr__(self):
        return f'{self.__class__.__name__}({self._repr_fields()})'


def power_spectrum_from_band(band: Band, engine: Engine) -> np.ndarray:
    """Return the power spectrum of a band."""
    assert band.power_spectrum is not None
    assert band.power_spectrum.parent_process is not None
    binary_content = engine.read_file(band.power_spectrum, binary=True)
    assert isinstance(binary_content, bytes)
    return np.frombuffer(binary_content)


def self_hartree_from_band(band: Band) -> np.ndarray:
    """Return the self-Hartree energy of a band."""
    return np.array([band.self_hartree])


def dummy_from_band(band: Band) -> np.ndarray:
    """Return a dummy descriptor (for use in the MeanEstimator)."""
    return np.array([np.nan])


def descriptor_from_band_factory(descriptor: str, engine: Engine | None = None) -> Callable[[Band], np.ndarray]:
    """Create a descriptor function from a string."""
    if descriptor == 'orbital_density':
        assert engine is not None
        return partial(power_spectrum_from_band, engine=engine)
    elif descriptor == 'self_hartree':
        return self_hartree_from_band
    else:
        raise ValueError(f"`{descriptor}` is not implemented as a valid descriptor.")


class MLModel(AbstractMLModel):
    """A ML model that contains a single estimator."""

    def __init__(self, estimator_type='ridge_regression', descriptor: str = 'orbital_density', is_trained: bool = False,
                 X_train: Optional[np.ndarray] = None, Y_train: Optional[np.ndarray] = None, **kwargs):
        super().__init__(estimator_type=estimator_type, descriptor=descriptor, is_trained=is_trained, **kwargs)
        self.X_train = X_train
        self.Y_train = Y_train

    @property
    def is_trained(self):
        """Return whether the estimator is trained."""
        return self.estimator.is_trained

    def _repr_fields(self):
        out = f"estimator={self.estimator.__class__.__name__}(...), is_trained={self.is_trained}"
        if self.X_train is not None:
            out += f", len(X_train)={len(self.X_train)}"
        if self.Y_train is not None:
            out += f", shape(Y_train)={self.Y_train.shape}"
        return out

    def predict(self, band: Band) -> float:
        """Make a prediction using the trained estimator."""
        if not self.is_trained:
            raise ValueError(f'`{self.__class__.__name__}` must be trained before calling `predict()`')

        descriptor = self.descriptor_from_band(band)

        return self.estimator.predict(descriptor)[0]

    def train(self):
        """Reset and train the estimator (including the StandardScaler)."""
        self.estimator = estimator_factory(self.estimator_type)
        self.estimator.fit(self.X_train, self.Y_train)

    def add_training_data(self, bands: List[Band]):  # x_train: np.ndarray, y_train: Union[float, np.ndarray]):
        """Add training data to the estimator."""
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

    def __init__(self, estimator_type='ridge_regression', descriptor='orbital_density', is_trained: bool = False,
                 X_train_occ: Optional[np.ndarray] = None, Y_train_occ: Optional[np.ndarray] = None,
                 X_train_emp: Optional[np.ndarray] = None, Y_train_emp: Optional[np.ndarray] = None,
                 estimator_occ=None, estimator_emp=None,
                 descriptor_from_band: Callable[[Band], np.ndarray] | None = None, engine: Engine | None = None):

        self.model_occ = MLModel(estimator_type=estimator_type, descriptor=descriptor, is_trained=is_trained,
                                 X_train=X_train_occ, Y_train=Y_train_occ, estimator=estimator_occ,
                                 descriptor_from_band=descriptor_from_band, engine=engine)
        self.model_emp = MLModel(estimator_type=estimator_type, descriptor=descriptor, is_trained=is_trained,
                                 X_train=X_train_emp, Y_train=Y_train_emp, estimator=estimator_emp,
                                 descriptor_from_band=descriptor_from_band, engine=engine)
        self.estimator_type = estimator_type
        self.descriptor_type = descriptor

    def add_training_data(self, bands: List[Band]) -> None:
        """Add training data to the sub-models."""
        for band in bands:
            if band.filled:
                self.model_occ.add_training_data([band])
            else:
                self.model_emp.add_training_data([band])

    def train(self) -> None:
        """Train both sub-models."""
        self.model_occ.train()
        self.model_emp.train()

    @property
    def is_trained(self) -> bool:
        """Return whether both sub-models are trained."""
        return self.model_occ.is_trained and self.model_emp.is_trained

    def predict(self, band: Band) -> float:
        """Predict the screening parameter of a given band."""
        if band.filled:
            return self.model_occ.predict(band)
        else:
            return self.model_emp.predict(band)

    def _repr_fields(self):
        out = f"estimator='{self.model_occ.__class__.__name__}(...)', is_trained={self.is_trained}"
        if self.model_occ.X_train is not None:
            out += f", len(X_train, occ)={len(self.model_occ.X_train)}"
        if self.model_occ.Y_train is not None:
            out += f", shape(Y_train, occ)={self.model_occ.Y_train.shape}"
        if self.model_emp.X_train is not None:
            out += f", len(X_train, emp)={len(self.model_emp.X_train)}"
        if self.model_emp.Y_train is not None:
            out += f", shape(Y_train, emp)={self.model_emp.Y_train.shape}"
        return out
