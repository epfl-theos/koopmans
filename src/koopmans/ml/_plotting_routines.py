from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np

from koopmans import utils


def plot_calculated_vs_predicted(y: np.ndarray, y_pred: np.ndarray, qoi: str, filename: Path,
                                 metric: Tuple[str, float]):
    """
    Plot the the calculated vs predicted alpha parameters (or other quantities of interest) for a given number of
    training samples.
    """
    _, ax = plt.subplots(1, 1, figsize=(7.5, 7.5))
    lb_x = 0.98 * min(y)
    ub_x = 1.02 * max(y)
    x = np.linspace(lb_x, ub_x, 1000)
    ax.set_xlim((lb_x, ub_x))
    ax.set_ylim((lb_x, ub_x))
    ax.plot(y, y_pred, 'o', color='green', label="{} = {:1.6f}".format(metric[0], metric[1]))
    ax.plot(x, x, color='grey', linewidth=0.75)
    plt.ylabel(f"predicted {qoi}")
    plt.xlabel(f"calculated {qoi}")
    plt.legend()
    utils.savefig(fname=str(filename))
    plt.close()


def plot_error_histogram(y: np.ndarray, y_pred: np.ndarray, qoi: str, filename: Path, metric: Tuple[str, float]):
    """
    Plot the error histogram of the alpha parameters (or other quantities of interest) for a given number of
    training samples.
    """
    error = np.abs(y-y_pred)
    _, ax = plt.subplots(1, 1, figsize=(7.5, 7.5))
    lb_x = 0.0
    ub_x = 0.5
    plt.xlabel(f"error")
    bins = np.linspace(lb_x, ub_x, 80)
    plt.xlim((lb_x, ub_x))
    ax.hist(error, bins,  label="{} = {:1.6f}".format(metric[0], metric[1]))
    plt.ylabel(f"number of {qoi}")
    plt.legend()
    utils.savefig(fname=str(filename))
    plt.close()


def plot_convergence(convergence_points: np.ndarray, means: np.ndarray, stddevs: np.ndarray, qoi: str,
                     metric_name: str, spin: int, filename: Path):
    """
    Plot the convergence of the MAE (or other error-metrics) of the alpha parameters (or other quantities of interest)
    as a function of the training samples.
    """

    _, ax = plt.subplots(1, 1, figsize=(15.0, 7.5))
    ax.xaxis.get_major_locator().set_params(integer=True)
    ax.set_xlabel("Number of snapshots for training")
    ax.set_ylabel(f"{metric_name} of predicted {qoi}")
    ax.errorbar(np.array(convergence_points)+1, means, stddevs, marker='o', linestyle="--")
    utils.savefig(fname=str(filename.with_suffix(".png")))
    plt.close()
