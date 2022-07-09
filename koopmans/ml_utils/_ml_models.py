import copy
from curses import has_key
from multiprocessing.sharedctypes import Value
from typing import List, Union
import numpy as np
from deepdiff import DeepDiff
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler


class RidgeRegression():
    def __init__(self, is_trained: bool = False, X_train: np.ndarray = None, Y_train: np.ndarray = None):
        self.is_trained = is_trained
        self.X_train = X_train
        self.Y_train = Y_train
        self.scaler = StandardScaler()
        self.model = Ridge(alpha=1.0)

    def __repr__(self):
        if self.X_train is None:
            return f'RidgeRegression(is_trained={self.is_trained}, no training data has been added so far)'
        else:
            return f'RidgeRegression(is_trained={self.is_trained},number_of_training_vectors={self.X_train.shape[0]},input_vector_dimension={self.Y_train.shape[0]})'

    def todict(self):
        # TODO Yannick: implement the todict-function
        # Shallow copy
        dct = dict(self.__dict__)
        items_to_pop = ['model', 'scaler']
        for item in items_to_pop:
            dct.pop(item)

        # # Adding information required by the json decoder
        # dct['__koopmans_name__'] = self.__class__.__name__
        # dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @classmethod
    def fromdict(cls, dct):
        return cls(**dct)

    def predict(self, x_test: np.ndarray):
        """
        Make a prediction of using the trained model. 
        """

        if self.is_trained:
            X_test = np.atleast_2d(x_test)
            X_test = self.scaler.transform(X_test)
            y_predict = self.model.predict(X_test)
            return y_predict
        else:
            return np.array([1.0])  # dummy value

    def train(self):
        """
        Reset the model and train the model (including the StandardScaler) with all training data added so far. 
        """

        self.model = Ridge(alpha=1.0)
        self.scaler = self.scaler.fit(self.X_train)
        X_train_scaled = self.scaler.transform(self.X_train)
        self.model.fit(X_train_scaled, self.Y_train)
        self.is_trained = True

    def add_training_data(self, x_train: np.ndarray, y_train: Union[float, np.ndarray]):
        """
        Add training data to the model.
        """

        x_train = np.atleast_2d(x_train)
        y_train = np.atleast_1d(y_train)
        if self.X_train is None:
            self.X_train = x_train
            self.Y_train = y_train
        else:
            self.X_train = np.concatenate([self.X_train, x_train])
            self.Y_train = np.concatenate([self.Y_train, y_train])

    def __eq__(self, other):
        items_to_pop = ['model', 'scaler']
        if isinstance(other, RidgeRegression):
            self_dict = copy.deepcopy(self.__dict__)
            other_dict = copy.deepcopy(other.__dict__)
            for item in items_to_pop:
                if item in other_dict:
                    other_dict.pop(item)
                if item in self_dict:
                    self_dict.pop(item)
            # raise ValueError()
            return DeepDiff(self_dict, other_dict, significant_digits=8, number_format_notation='e') == {}
        else:
            return False
