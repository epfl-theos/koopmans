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
from koopmans.ml_utils import MLModel
from koopmans.ml_utils._plotting_routines import plot_calculated_vs_predicted, plot_error_histogram, plot_convergence

from ._workflow import Workflow


class ConvergenceMLWorkflow(Workflow):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def convert_to_list(self, param, type):
        if isinstance(param, type):
            return [param]
        else:
            return param

    def _run(self) -> None:

        n_maxs = self.convert_to_list(self.parameters.n_max, int)
        l_maxs = self.convert_to_list(self.parameters.l_max, int)
        r_mins = self.convert_to_list(self.parameters.r_min, float)
        r_maxs_added = self.convert_to_list(self.parameters.r_max, float)

        # Yannick TODO: find more elegant solution to prevent ML model in the first run
        self.parameters.n_max = 2
        self.parameters.l_max = 2
        self.parameters.r_min = 4.0
        self.parameters.r_max = 6.0
        if self.parameters.number_of_training_snapshots >= len(self.snapshots):
            raise ValueError(
                "There are not enough test-snapshots available for a convergence_analysis. Please increase the number of snapshots in the xyz-file or decrease 'number_of_training_snapshots'")
        delete_final_dir = False
        if 'evs' in self.quantities_of_interest:
            delete_final_dir = True

        tmp_number_of_training_snapshots = self.parameters.number_of_training_snapshots
        # set the number of training snapshots to a very high value to not use the prediction for the first run
        self.parameters.number_of_training_snapshots = 10000
        self.print(
            f'Obtaining ab-initio results for the last {len(self.test_indices)} snapshot(s)', style='heading')
        # get the ab-initio result for the test_indices
        twf = TrajectoryWorkflow(snapshots=self.snapshots, **self.wf_kwargs)
        self.run_subworkflow(twf, indices=self.test_indices,
                             save_dir=self.dirs['convergence_true'], delete_final_dir=delete_final_dir)
        # set the number of training snapshots back to its original value
        self.parameters.number_of_training_snapshots = tmp_number_of_training_snapshots

        from_scratch = False
        # TODO Yannick: rewrite this workflow s.t. it is also sensible to do grid search with more than one snapshot to train on
        f = open("result_grid_search.txt", "w")
        f.write('n_max\tl_max\tr_min\tr_cut\tMAE\n')
        f.close()
        for self.parameters.n_max in n_maxs:
            for self.parameters.l_max in l_maxs:
                for self.parameters.r_min in r_mins:
                    for r_max_added in r_maxs_added:
                        self.parameters.r_max = self.parameters.r_min + r_max_added
                        if self.parameters.occ_and_emp_together:
                            self.ml_model = MLModel(self.parameters.type_ml_model)
                        else:
                            self.ml_model_occ = MLModel(self.parameters.type_ml_model)
                            self.ml_model_emp = MLModel(self.parameters.type_ml_model)
                        for convergence_point in self.convergence_points:
                            train_indices = [convergence_point]
                            self.print(
                                f'Training on {len(train_indices)} snapshot(s) and then testing on the last {len(self.test_indices)} snapshot(s)', style='heading')
                            twf = TrajectoryWorkflow(snapshots=self.snapshots, **self.wf_kwargs)
                            self.run_subworkflow(twf, indices=train_indices)  # train the model
                            # test the model (without retraining the model)
                            twf = TrajectoryWorkflow(snapshots=self.snapshots, **self.wf_kwargs)
                            self.run_subworkflow(twf, from_scratch=from_scratch, indices=self.test_indices,
                                                 save_dir=self.dirs[f'convergence_{convergence_point}'], delete_final_dir=delete_final_dir)
                        self.get_result_dict()
                        self.make_plots_from_result_dict()
                        f = open("result_grid_search.txt", "a")
                        f.write("{:5.4f}\t{:5.4f}\t{:5.4f}\t{:5.4f}\t{:5.8f}\n".format(self.parameters.n_max,
                                self.parameters.l_max, self.parameters.r_min, self.parameters.r_max, (self.result_dict[str("spin_"+str(0))]['alphas']['MAE']['mean'][-1])))
                        f.close()

    def make_plots_from_result_dict(self):

        for spin in range(self.bands.n_spin):
            spin_id = str("spin_"+str(spin))
            for qoi in self.quantities_of_interest:
                # for i, convergence_point in enumerate(self.convergence_points):
                #     res = self.result_dict[spin_id][qoi]
                #     plot_calculated_vs_predicted(res['true_array'].flatten(), res['pred_array'][i, :, :].flatten(
                #     ), qoi, self.dirs[f'convergence_figures_{i}'] / f"{spin_id}_{qoi}_calculated_vs_predicted.png", ('MAE', res['MAE']['mean'][i]))
                #     plot_error_histogram(res['true_array'].flatten(), res['pred_array'][i, :, :].flatten(
                #     ), qoi, self.dirs[f'convergence_figures_{i}'] / f"{spin_id}_{qoi}_error_histogram.png", ('MAE', res['MAE']['mean'][i]))
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
        wf.quantities_of_interest = ['alphas']  # , 'evs']
        wf.metrics = {'MAE': mean_absolute_error, 'R2S': r2_score}
        wf.statistics = {'mean': np.mean, 'stdd': np.std}

        return wf

    def todict(self):

        dct = dict(self.__dict__)

        # items_to_pop = ['atoms']
        # for item in items_to_pop:
        #     dct.pop(item)

        items_to_replace_with_strings = ['metrics', 'statistics']
        for item in items_to_replace_with_strings:
            for fct_name in dct[item].keys():
                # replace the acutal function with just the corresponding string to make it JSON-writeable
                dct[item][fct_name] = fct_name
        return dct
