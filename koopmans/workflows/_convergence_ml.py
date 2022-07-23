from ._trajectory import TrajectoryWorkflow
from abc import ABC, abstractmethod

import copy
import json as json_ext
import os
import statistics
from pathlib import Path
from tracemalloc import Snapshot
from typing import Any, Dict, List, Optional, Tuple, Union
from unittest import result

import matplotlib.pyplot as plt
import numpy as np
from sklearn import metrics
from sklearn.metrics import mean_absolute_error, r2_score

from ase import Atoms, io
from ase.calculators.calculator import Calculator
from koopmans import calculators, utils
from koopmans.ml_utils import RidgeRegression

from ._workflow import Workflow


class ConvergenceMLWorkflow(Workflow):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _run(self) -> None:
        if self.parameters.number_of_training_snapshots >= len(self.snapshots):
            raise ValueError(
                "There are not enough test-snapshots available for a convergence_analysis. Please increase the number of snapshots in the xyz-file or decrease 'number_of_training_snapshots'")

        workflow = TrajectoryWorkflow(self.snapshots, **self.wf_kwargs)

        # Yannick TODO: find more elegant solution
        tmp_number_of_training_snapshots = self.parameters.number_of_training_snapshots
        # set the number of training snapshots to a very high value to not use the prediction for the first run
        self.parameters.number_of_training_snapshots = 10000
        self.print(
            f'Obtaining ab-initio results for the last {len(self.test_indices)} snapshot(s)', style='heading')
        # get the ab-initio result for the test_indices
        wf_kwargs = {'indices': self.test_indices, 'save_dir': self.dirs['convergence_true']}
        self.run_subworkflow(workflow, **wf_kwargs)
        # set the number of training snapshots back to its original value
        self.parameters.number_of_training_snapshots = tmp_number_of_training_snapshots
        self.ml_model = RidgeRegression()
        for convergence_point in self.convergence_points:
            train_indices = [convergence_point]
            self.print(
                f'Training on {len(train_indices)} snapshot(s) and then testing on the last {len(self.test_indices)} snapshot(s)', style='heading')
            wf_kwargs = {'indices': train_indices}
            self.run_subworkflow(workflow, **wf_kwargs)  # train the model
            # test the model (without retraining the model)
            delete_final_dir = False
            if 'evs' in self.quantities_of_interest:
                delete_final_dir = True
            wf_kwargs = {'indices': self.test_indices,
                         'save_dir': self.dirs[f'convergence_{convergence_point}'],  'delete_final_dir': delete_final_dir}
            self.run_subworkflow(workflow, **wf_kwargs)
        self.get_result_dict()
        # self.make_plots_from_result_dict()

    def make_plots_from_result_dict(self):

        for spin in range(self.bands.n_spin):
            spin_id = str("spin_"+str(spin))
            for qoi in self.quantities_of_interest:
                for i, convergence_point in enumerate(self.convergence_points):
                    res = self.result_dict[spin_id][qoi]
                    plot_calculated_vs_predicted(res['true_array'].flatten(), res['pred_array'][i, :, :].flatten(
                    ), qoi, self.dirs[f'convergence_figures_{i}'] / f"{spin_id}_{qoi}_calculated_vs_predicted.png", ('MAE', res['MAE']['mean'][i]))
                    plot_error_histogram(res['true_array'].flatten(), res['pred_array'][i, :, :].flatten(
                    ), qoi, self.dirs[f'convergence_figures_{i}'] / f"{spin_id}_{qoi}_error_histogram.png", ('MAE', res['MAE']['mean'][i]))
                for metric in self.metrics:
                    res = self.result_dict[spin_id][qoi][metric]
                    plot_convergence(self.convergence_points, res['mean'], res['stdd'], qoi, metric,
                                     spin, self.dirs['convergence_figures'] / f"{spin_id}_{qoi}_{metric}_convergence")

    def get_result_dict(self):
        result_dict = {}

        for spin in range(self.bands.n_spin):
            spin_id = str("spin_"+str(spin))
            result_dict[spin_id] = {}
            for qoi in self.quantities_of_interest:
                result_dict[spin_id][qoi] = {}
                result_dict[spin_id][qoi]['pred_array'] = np.zeros(
                    (len(self.convergence_points), len(self.test_indices), self.bands.n_bands[spin]))
                result_dict[spin_id][qoi]['true_array'] = np.zeros(
                    (len(self.test_indices), self.bands.n_bands[spin]))
                for metric_name in self.metrics.keys():
                    result_dict[spin_id][qoi][metric_name] = {}
                    for statistic_name in self.statistics.keys():
                        result_dict[spin_id][qoi][metric_name][statistic_name] = np.zeros(
                            len(self.convergence_points))

        self.result_dict = result_dict

        for spin in range(self.bands.n_spin):
            spin_id = str("spin_"+str(spin))
            for qoi in self.quantities_of_interest:
                for j, test_index in enumerate(self.test_indices):
                    self.result_dict[spin_id][qoi]['true_array'][j, :] = np.loadtxt(
                        self.dirs['convergence_true'] / f"{qoi}_snapshot_{test_index+1}.txt")[spin, :]
                for i, convergence_point in enumerate(self.convergence_points):
                    for metric in self.metrics.keys():
                        tmp_array = np.zeros(len(self.test_indices))
                        for j, test_index in enumerate(self.test_indices):
                            self.result_dict[spin_id][qoi]['pred_array'][i, j, :] = np.loadtxt(
                                self.dirs[f'convergence_{convergence_point}'] / f"{qoi}_snapshot_{test_index+1}.txt")[spin, :]

                            tmp_array[j] = self.metrics[metric](
                                self.result_dict[spin_id][qoi]['pred_array'][i, j, :], self.result_dict[spin_id][qoi]['true_array'][j, :])

                        for statistic in self.statistics:
                            self.result_dict[spin_id][qoi][metric][statistic][i] = self.statistics[statistic](tmp_array)

    @ classmethod
    def _fromjsondct(cls, bigdct: Dict[str, Any]):
        """
        Reads the atomic positions for each snapshot from the xyz-file specified by the user in the snapshots-file.
        """

        try:
            snapshots_file = bigdct['setup'].pop('snapshots')
        except:
            raise ValueError(
                f'To calculate a trajectory, please provide a xyz-file containing the atomic positions of the snapshots in the setup-block of the json-input file.')

        snapshots = io.read(snapshots_file, index=':')
        if isinstance(snapshots, Atoms):
            snapshots = [snapshots]
        bigdct['setup']['atomic_positions'] = utils.construct_atomic_positions_block(snapshots[0])
        wf = super(ConvergenceMLWorkflow, cls)._fromjsondct(bigdct)
        wf.snapshots = snapshots
        wf.number_of_snapshots = len(snapshots)

        wf.convergence_points = list(range(0, wf.parameters.number_of_training_snapshots))
        wf.test_indices = list(range(wf.parameters.number_of_training_snapshots,
                                     wf.number_of_snapshots))

        convergence_dir = Path.cwd() / 'convergence_analysis'
        utils.system_call(f'rm -rf {convergence_dir}')
        wf.dirs = {
            'convergence': convergence_dir,
            'convergence_true': convergence_dir / 'true',
            'convergence_pred': convergence_dir / 'pred',
            'convergence_figures': convergence_dir / 'figures'
        }
        for convergence_point in wf.convergence_points:
            wf.dirs[f'convergence_{convergence_point}'] = wf.dirs[f'convergence_pred'] / \
                f"predicted_after_{convergence_point+1}"
            wf.dirs[f'convergence_figures_{convergence_point}'] = wf.dirs[f'convergence_figures'] / \
                f"predicted_after_{convergence_point+1}"
        for dir in wf.dirs.values():
            utils.system_call(f'mkdir -p {dir}')
        wf.quantities_of_interest = ['alphas', 'evs']
        wf.metrics = {'MAE': mean_absolute_error, 'R2S': r2_score}
        wf.statistics = {'mean': np.mean, 'stdd': np.std}

        return wf

    def todict(self):

        dct = dict(self.__dict__)

        # items_to_pop = ['atoms']
        # for item in items_to_pop:
        #     dct.pop(item)

        items_to_pop = ['metrics', 'statistics']
        for item in items_to_pop:
            dct.pop(item)
        return dct


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
