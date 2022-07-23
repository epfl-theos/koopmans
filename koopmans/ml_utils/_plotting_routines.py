from pathlib import Path
from typing import Tuple
import matplotlib.pyplot as plt
from koopmans import utils
import numpy as np


def plot_calculated_vs_predicted(y: np.ndarray, y_pred: np.ndarray, qoi: str, filename: Path, metric: Tuple[str, float]):
    fig, ax = plt.subplots(1, 1, figsize=(7.5, 7.5))
    lb_x = 0.995*min(y)
    ub_x = 1.005*max(y)
    # if qoi == 'alphas':
    #     lb_x = 0.995*min(y)
    #     ub_x = 1.005*max(y)
    # elif qoi == 'evs':
    #     lb_x = 0.995*min(y)
    #     ub_x = 1.005*max(y)
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
    error = np.abs(y-y_pred)
    fig, ax = plt.subplots(1, 1, figsize=(7.5, 7.5))
    lb_x = 0.0
    ub_x = 10*max(error)
    plt.xlabel(f"error")
    if qoi == 'alphas':
        lb_x = 0.0
        ub_x = 10*max(error)
    elif qoi == 'evs':
        lb_x = 0.0
        ub_x = 0.7
        plt.xlabel(f"error in eV")
    bins = np.linspace(lb_x, ub_x, 100)
    # ax.set_xlim((lb_x, ub_x))
    # ax.set_ylim((0, 1000))
    ax.hist(error, bins,  label="{} = {:1.6f}".format(metric[0], metric[1]))
    plt.ylabel(f"number of {qoi}")
    plt.legend()
    utils.savefig(fname=str(filename))
    plt.close()


def plot_convergence(convergence_points: np.ndarray, means: np.ndarray, stddevs: np.ndarray, qoi: str, metric_name: str, spin: int, filename: Path):
    np.savetxt(filename.with_suffix(".txt"), np.array([convergence_points, means, stddevs]).T)

    # def snap2orb(x):
    #     return x * self.bands.n_bands[spin]

    # def orb2snap(x):
    #     return x/self.bands.n_bands[spin]

    fig, ax = plt.subplots(1, 1, figsize=(15.0, 7.5))
    lb_y = 0.0
    ub_y = 4*max(means)
    # if qoi == 'alphas':
    #     lb_y = 0.0
    #     ub_y = 0.0175
    # elif qoi == 'evs':
    #     lb_y = 0.0
    #     ub_y = 0.0175

    ax.set_ylim(lb_y, ub_y)
    ax.xaxis.get_major_locator().set_params(integer=True)
    ax.set_xlabel("Number of snapshots for training")
    ax.set_ylabel(f"{metric_name} of predicted {qoi}")
    ax.errorbar(convergence_points, means, stddevs, marker='o', linestyle="--")
    # secax = ax.secondary_xaxis('top', functions=(snap2orb, orb2snap))
    # secax.set_xlabel('Number of orbitals for training')
    utils.savefig(fname=str(filename.with_suffix(".png")))
    plt.close()
