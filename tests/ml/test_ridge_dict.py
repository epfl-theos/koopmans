"""
Test for to- and fromdict functions of the RidgeRegression class
"""
import numpy as np
from pathlib import Path

from koopmans import utils
from koopmans.ml import MLModel


def test_ridge_write_then_read_empty_from_dict(tmpdir):
    with utils.chdir(tmpdir):
        ridge_out = MLModel('ridge_regression')
        ridge_dict = ridge_out.todict()
        ridge_in = MLModel.fromdict(ridge_dict)

        assert ridge_out == ridge_in


def test_ridge_write_then_read_non_empty_from_dict(tmpdir):
    with utils.chdir(tmpdir):
        ridge_out = MLModel('ridge_regression')
        X = np.array([[2.7, 5.9, 4.9], [4.2, 4.9, 6.9]])
        y = np.array([2.9, 7.0])
        ridge_out.add_training_data(X, y)
        ridge_out.train()
        ridge_dict = ridge_out.todict()

        ridge_in = MLModel.fromdict(ridge_dict)

        assert ridge_out == ridge_in


def test_ridge_equality_method(tmpdir):
    with utils.chdir(tmpdir):
        ridge_out = MLModel('ridge_regression')
        X = np.array([[2.7, 5.9, 4.9], [4.2, 4.9, 6.9]])
        y = np.array([2.9, 7.0])
        ridge_out.add_training_data(X, y)
        ridge_out.train()

        ridge_in = MLModel('ridge_regression')
        X = np.array([[2.7, 5.8, 4.9], [4.2, 4.9, 6.9]])
        y = np.array([2.9, 7.0])
        ridge_in.add_training_data(X, y)
        ridge_in.train()

        assert (ridge_out == ridge_in) == False


def test_ridge_save_to_and_load_non_empty_model_from_file(tmpdir):
    with utils.chdir(tmpdir):
        ridge_out = MLModel('ridge_regression')
        X = np.array([[2.7, 5.9, 4.9], [4.2, 4.9, 6.9]])
        X_test = np.array([[2.4, 5.7, 4.6], [4.4, 4.6, 6.2]])
        y = np.array([2.9, 7.0])
        ridge_out.add_training_data(X, y)
        ridge_out.train()
        y_pred_out = ridge_out.predict(X_test)
        save_dir = Path('.') / 'trained_model'
        save_dir.mkdir(parents=True, exist_ok=True)
        ridge_out.model_class.save_to_file(save_dir)

        ridge_in = MLModel('ridge_regression', save_dir=save_dir)
        y_pred_in = ridge_in.predict(X_test)

        assert np.all(y_pred_out == y_pred_in)


def test_ridge_save_to_and_load_empty_model_from_file(tmpdir):
    with utils.chdir(tmpdir):
        ridge_out = MLModel('ridge_regression')
        save_dir = Path('.') / 'trained_model'
        save_dir.mkdir(parents=True, exist_ok=True)
        ridge_out.model_class.save_to_file(save_dir)

        ridge_in = MLModel('ridge_regression', save_dir=save_dir)

        X = np.array([[2.7, 5.9, 4.9], [4.2, 4.9, 6.9]])
        X_test = np.array([[2.4, 5.7, 4.6], [4.4, 4.6, 6.2]])
        y = np.array([2.9, 7.0])

        ridge_out.add_training_data(X, y)
        ridge_out.train()
        y_pred_out = ridge_out.predict(X_test)

        ridge_in.add_training_data(X, y)
        ridge_in.train()
        y_pred_in = ridge_in.predict(X_test)


        assert np.all(y_pred_out == y_pred_in)


