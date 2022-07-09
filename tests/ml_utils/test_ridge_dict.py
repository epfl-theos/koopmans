"""
Test for to- and fromdict functions of the RidgeRegression class
"""
from koopmans import utils
from koopmans.ml_utils import RidgeRegression
import numpy as np


def test_ridge_write_then_read_empty_from_dict(tmpdir):
    with utils.chdir(tmpdir):
        ridge_out = RidgeRegression()
        ridge_dict = ridge_out.todict()
        ridge_in = RidgeRegression.fromdict(ridge_dict)

        assert ridge_out == ridge_in


def test_ridge_write_then_read_non_empty_from_dict(tmpdir):
    with utils.chdir(tmpdir):
        ridge_out = RidgeRegression()
        X = np.array([[2.7, 5.9, 4.9], [4.2, 4.9, 6.9]])
        y = np.array([2.9, 7.0])
        ridge_out.add_training_data(X, y)
        ridge_out.train()
        ridge_dict = ridge_out.todict()

        ridge_in = RidgeRegression.fromdict(ridge_dict)

        assert ridge_out == ridge_in


def test_ridge_equality_method(tmpdir):
    with utils.chdir(tmpdir):
        ridge_out = RidgeRegression()
        X = np.array([[2.7, 5.9, 4.9], [4.2, 4.9, 6.9]])
        y = np.array([2.9, 7.0])
        ridge_out.add_training_data(X, y)
        ridge_out.train()

        ridge_in = RidgeRegression()
        X = np.array([[2.7, 5.8, 4.9], [4.2, 4.9, 6.9]])
        y = np.array([2.9, 7.0])
        ridge_in.add_training_data(X, y)
        ridge_in.train()

        assert (ridge_out == ridge_in) == False
