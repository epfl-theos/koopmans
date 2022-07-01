import numpy as np
from sklearn.linear_model import Ridge 
from sklearn.preprocessing import StandardScaler
# from sklearn.metrics import r2_score
# from sklearn.metrics import mean_absolute_error as mae

class RidgeRegression():
    def __init__(self):
        self.is_trained = False
        self.X_train    = []
        self.Y_train    = []
        self.scaler     = StandardScaler()
        self.model      = Ridge(alpha=1.0)
    
    def __repr__(self):
        return f'RidgeRegression(is_trained={self.is_trained})'

    def todict(self):
        # TODO: implement
        return 

    def predict(self, x_test):
        if self.is_trained:
            X_test = np.atleast_2d(x_test)
            X_test = self.scaler.transform(X_test)
            print("Yannick Debug: shape(X_test) = ", np.shape(X_test))
            y_predict = self.model.predict(X_test)
            return y_predict
        else:
            return np.array([1.0]) # dummy value 
            # return np.array([np.nan])

    def train(self):
        self.model      = Ridge(alpha=1.0)
        self.scaler     = self.scaler.fit(self.X_train)
        X_train_scaled  = self.scaler.transform(self.X_train)

        print("Yannick Debug: shape(X_train) = ", np.shape(self.X_train))
        print("Yannick Debug: shape(Y_train) = ", np.shape(self.Y_train))
        self.model.fit(X_train_scaled,self.Y_train)
        self.is_trained = True
    
    def add_training_data(self, x_train, y_train):
        x_train = np.atleast_2d(x_train)
        y_train = np.atleast_1d(y_train)
        if self.X_train == []:
            self.X_train = x_train
            self.Y_train = y_train
        else:    
            self.X_train = np.concatenate([self.X_train, x_train])
            self.Y_train = np.concatenate([self.Y_train, y_train])
        
        
        
