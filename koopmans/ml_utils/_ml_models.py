from typing import Union
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler


class RidgeRegression():
    def __init__(self):
        self.is_trained = False
        self.X_train = []
        self.Y_train = []
        self.scaler = StandardScaler()
        self.model = Ridge(alpha=1.0)

    def __repr__(self):
        return f'RidgeRegression(is_trained={self.is_trained},number_of_training_vectors={self.X_train.shape[0]},input_vector_dimension={self.Y_train.shape[0]})'

    def todict(self):
        # TODO: implement
        return

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
        if self.X_train == []:
            self.X_train = x_train
            self.Y_train = y_train
        else:
            self.X_train = np.concatenate([self.X_train, x_train])
            self.Y_train = np.concatenate([self.Y_train, y_train])
