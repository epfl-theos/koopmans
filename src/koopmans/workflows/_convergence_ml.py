import itertools
import pickle
import shutil
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from ase import Atoms, io
from sklearn.metrics import mean_absolute_error, r2_score

from koopmans import utils
from koopmans.ml import (MLModel, plot_calculated_vs_predicted,
                         plot_convergence, plot_error_histogram,
                         precompute_parameters_of_radial_basis)

from ._trajectory import TrajectoryWorkflow
from ._workflow import Workflow


class ConvergenceMLWorkflow(Workflow):

    def __init__(self, snapshots: List[Atoms], * args, **kwargs):

        if isinstance(snapshots, Atoms):
            kwargs['atoms'] = snapshots
        else:
            kwargs['atoms'] = snapshots[0]

        super().__init__(*args, **kwargs)

        self.snapshots = snapshots
        self.number_of_snapshots = len(self.snapshots)
        self.convergence_points = list(range(0, self.ml.number_of_training_snapshots))
        self.test_indices = list(range(self.ml.number_of_training_snapshots,
                                       self.number_of_snapshots))

    @ classmethod
    def _fromjsondct(cls, bigdct: Dict[str, Any], override: Dict[str, Any] = {}):
        """
        Reads the atomic positions for each snapshot from the xyz-file specified by the user in the snapshots-file
        and initialize the snapshot-indices for training and testing
        """
        try:
            snapshots_file = bigdct['atoms']['atomic_positions'].pop('snapshots')
        except:
            raise ValueError(f'To calculate a trajectory, please provide a xyz-file containing the atomic positions '
                             'of the snapshots in the setup-block of the json-input file.')

        snapshots = io.read(snapshots_file, index=':')

        if isinstance(snapshots, Atoms):
            snapshots = [snapshots]
        bigdct['atoms']['atomic_positions'] = utils.construct_atomic_positions_block(snapshots[0])

        wf = super(ConvergenceMLWorkflow, cls)._fromjsondct(bigdct, override)

        # Initialize the set of snapshots for training and testing
        wf.snapshots = snapshots
        wf.number_of_snapshots = len(snapshots)
        wf.convergence_points = list(range(0, wf.ml.number_of_training_snapshots))
        wf.test_indices = list(range(wf.ml.number_of_training_snapshots,
                                     wf.number_of_snapshots))

        return wf

    def toinputjson(self) -> Dict[str, Dict[str, Any]]:
        bigdct = super().toinputjson()
        snapshots_file = "snapshots.json"
        io.write(snapshots_file, self.snapshots)
        bigdct['atoms']['atomic_positions'] = {"snapshots": snapshots_file}
        return bigdct

    def todict(self):

        dct = dict(self.__dict__)

        items_to_pop = ['atoms']
        for item in items_to_pop:
            dct.pop(item)

        return dct

    def inititalize_directories(self, grid_search_mode: bool):
        '''
        Delete the old folder of the convergence analysis and initialize
        all new folders needed for the convergence analysis.
        '''
        convergence_dir = Path.cwd() / 'convergence_analysis'
        shutil.rmtree(convergence_dir, ignore_errors=True)
        self.dirs = {
            'convergence': convergence_dir,  # parent convergence_directory
            'convergence_true': convergence_dir / 'true',  # directory where we store the ab-initio results
            'convergence_pred': convergence_dir / 'pred',  # folder where we store the predicted results
            'convergence_final_results': convergence_dir / 'final_results'  # folder where we store final resilts
        }
        for convergence_point in self.convergence_points:
            # for the grid search we are only interested in the result after all training snapshots have been added
            if not (grid_search_mode and convergence_point != self.convergence_points[-1]):
                self.dirs[f'convergence_{convergence_point}'] = self.dirs[f'convergence_pred'] / \
                    f"predicted_after_{convergence_point+1}"
            # in case of the grid search we won't produce all the plots, hence we also don't need the corresponding folders
            if not grid_search_mode:
                self.dirs[f'convergence_final_results_{convergence_point}'] = self.dirs[f'convergence_final_results'] / \
                    f"predicted_after_{convergence_point+1}"

        # create all directories
        for dir in self.dirs.values():
            dir.mkdir(parents=True, exist_ok=True)

    def _run(self) -> None:
        """
        Depending on the provided input parameters this workflow performs either of the following:

        i)  If for at least one the hyper parameters (n_max, l_max, r_min, r_max) a list with more than one number is
            provided, we perform a grid search w.r.t. these parameters. That means that the MAE (or other error metrics)
            of the alpha parameters (or other quantities of interest) is computed for every valid (r_min<r_max)
            combination of (n_max, l_max, r_min, r_max) when trained on (1, ... , number_of_training_samples) and
            evaluated on the remaining snapshots (number_of_training_samples, ... , number of snapshots provided in the
            xyz file).

        ii) Else, a convergence analyis with respect to the number of training samples is performed.  That means that
            the MAE (or other error-metrics) of the alpha parameters (or other quantities of interest) are calculated
            when trained on (1, ... , number_of_training_samples) and evaluated on the remaining snapshots
            (number_of_training_samples, ... , number of snapshots provided in the xyz file).
        """
        if self.ml.number_of_training_snapshots >= len(self.snapshots):
            raise ValueError(
                "There are not enough test-snapshots available for a convergence_analysis. Please increase the number "
                "of snapshots in the xyz-file or decrease 'number_of_training_snapshots'")

        def convert_to_list(param, type):
            if isinstance(param, type):
                # If param is an int or a float convert it for the checks to a list
                return [param]
            else:
                # If param is not an int or a float check that it is a list of ints / floats
                assert(isinstance(param, list))
                for value in param:
                    assert(isinstance(value, type))
                return param

        n_maxs = convert_to_list(self.ml.n_max, int)
        l_maxs = convert_to_list(self.ml.l_max, int)
        r_mins = convert_to_list(self.ml.r_min, float)
        r_maxs = convert_to_list(self.ml.r_max, float)

        # Perform a grid search iff there is at least one parameter that contains more than one number
        if len(n_maxs)+len(l_maxs)+len(r_mins)+len(r_maxs) > 4:
            grid_search_mode = True
        else:
            grid_search_mode = False

        self.inititalize_directories(grid_search_mode)

        # If the eigenvalues are quantities of interest we have to perform the final calcualtion afresh for every
        # snapshot since its result (and thereby the eigenvalues) depend on the prediction of the alpha parameters
        get_evs = False
        if 'evs' in self.ml.quantities_of_interest:
            get_evs = True

        # Get the ab-initio result for the test_indices
        self.print(f'Obtaining ab-initio results for the last {len(self.test_indices)} snapshot(s)', style='heading')
        # Make sure that we compute these snapshots ab-initio and don't use the ml-predicted alpha values
        use_ml = self.ml.use_ml
        self.ml.use_ml = False
        twf = TrajectoryWorkflow.fromparent(self, snapshots=self.snapshots, indices=self.test_indices,
                                            save_dir=self.dirs['convergence_true'], get_evs=get_evs,
                                            overwrite_atoms=False)
        twf.run()
        self.ml.use_ml = use_ml

        from_scratch = self.parameters.from_scratch

        with open(self.dirs['convergence'] / "result_grid_search.txt", "w") as fd:
            fd.write('n_max\tl_max\tr_min\tr_max\tMAE\n')

        # Loop over all possible combinations
        for (self.ml.n_max, self.ml.l_max, self.ml.r_min, self.ml.r_max) in \
                itertools.product(n_maxs, l_maxs, r_mins, r_maxs):

            if self.ml.r_max <= self.ml.r_min:
                # Skip this combination if r_max<r_min
                continue
            try:
                precompute_parameters_of_radial_basis(
                    self.ml.n_max, self.ml.l_max, self.ml.r_min, self.ml.r_max)
            except:
                # skip this set of parameters if it is not possible to find the coefficients of the radial basis
                # function
                utils.warn(
                    f"Failed to precompute the radial basis. You might want to try a larger r_min, e.g. r_min=1.0.")
                continue

            self.print(f"Running the calculation with (n_max, l_max, r_min, r_max) = ({self.ml.n_max}, "
                       f"{self.ml.l_max}, {self.ml.r_min}, {self.ml.r_max})")

            # Reset the MLModel for every new combination of (n_max, l_max, r_min, r_max)
            if self.ml.occ_and_emp_together:
                self.ml.ml_model = MLModel(self.ml.type_of_ml_model)
            else:
                self.ml.ml_model_occ = MLModel(self.ml.type_of_ml_model)
                self.ml.ml_model_emp = MLModel(self.ml.type_of_ml_model)

            # train the model on (1, ... , number_of_training_snapshots) samples
            for convergence_point in self.convergence_points:
                self.print(f'Adding snapshot {convergence_point+1} to the training data', style='heading')

                # this is the snapshot we want to add to our training data
                train_indices = [convergence_point]

                # Train on this snapshot
                twf = TrajectoryWorkflow.fromparent(
                    self, snapshots=self.snapshots, indices=train_indices, overwrite_atoms=False)
                # we want to make sure that no snapshot is calculated afresh just because the final directory of the
                # preceding calculation was deleted
                twf.run(from_scratch=from_scratch)

                # for the grid search we are only interested in the result after all training snapshots have been added
                if not (grid_search_mode and convergence_point != self.convergence_points[-1]):

                    # Test the model (and recalculate the final calculation in case we are interested in eigenvalues)
                    # and save the results
                    self.print(f'Testing on the last {len(self.test_indices)} snapshot(s)', style='heading')
                    twf = TrajectoryWorkflow.fromparent(self, snapshots=self.snapshots, indices=self.test_indices,
                                                        save_dir=self.dirs[f'convergence_{convergence_point}'],
                                                        get_evs=get_evs, overwrite_atoms=False)
                    twf.run(from_scratch=False)

            # Gather all the important results
            self.get_result_dict()

            if not grid_search_mode:
                # Create the plots for the convergence analysis
                self.make_convergence_analysis_plots()
            else:
                # Add the result to the result-file of the grid search
                val = self.result_dict[str("spin_"+str(0))]['alphas']['MAE']['mean'][-1]
                with open(self.dirs['convergence'] / "result_grid_search.txt", "a") as fd:
                    fd.write(f"{self.ml.n_max:5.4f}\t{self.ml.l_max:5.4f}\t{self.ml.r_min:5.4f}\t"
                             f"{self.ml.r_max:5.4f}\t{val:5.8f}\n")

    def get_result_dict(self):
        '''
        Create a directory that contains all relevant results. This intermediate step is supposed to simplify
        debugging and postprocessing.
        '''

        # Define the quantities we are interested in (calculating the eigenvalues requires performing the final
        # calculation afresh for every snapshot and is hence a little bit more computationally expensive)
        metrics = {'MAE': mean_absolute_error, 'R2S': r2_score}
        statistics = {'mean': np.mean, 'stdd': np.std}

        result_dict = {}

        for spin in range(self.bands.n_spin):
            spin_id = str("spin_"+str(spin))
            result_dict[spin_id] = {}
            for qoi in self.ml.quantities_of_interest:
                result_dict[spin_id][qoi] = {}
                result_dict[spin_id][qoi]['pred_array'] = np.zeros(
                    (len(self.convergence_points), len(self.test_indices), self.bands.n_bands[spin]))
                result_dict[spin_id][qoi]['true_array'] = np.zeros(
                    (len(self.test_indices), self.bands.n_bands[spin]))
                for metric_name in metrics.keys():
                    result_dict[spin_id][qoi][metric_name] = {}
                    for statistic_name in statistics.keys():
                        result_dict[spin_id][qoi][metric_name][statistic_name] = np.zeros(
                            len(self.convergence_points))

        self.result_dict = result_dict

        for spin in range(self.bands.n_spin):
            spin_id = str("spin_"+str(spin))
            for qoi in self.ml.quantities_of_interest:
                for j, test_index in enumerate(self.test_indices):
                    self.result_dict[spin_id][qoi]['true_array'][j, :] = np.loadtxt(
                        self.dirs['convergence_true'] / f"{qoi}_snapshot_{test_index+1}.txt")[spin, :]
                for i, convergence_point in enumerate(self.convergence_points):
                    for metric in metrics.keys():
                        tmp_array = np.zeros(len(self.test_indices))
                        for j, test_index in enumerate(self.test_indices):
                            self.result_dict[spin_id][qoi]['pred_array'][i, j, :] = np.loadtxt(
                                self.dirs[f'convergence_{convergence_point}']
                                / f"{qoi}_snapshot_{test_index+1}.txt")[spin, :]

                            tmp_array[j] = metrics[metric](self.result_dict[spin_id][qoi]['pred_array'][i, j, :],
                                                           self.result_dict[spin_id][qoi]['true_array'][j, :])

                        for statistic in statistics:
                            self.result_dict[spin_id][qoi][metric][statistic][i] = statistics[statistic](tmp_array)

    def make_convergence_analysis_plots(self):
        """
        Create all relevant plots for the convergence analysis
        """
        metrics = {'MAE': mean_absolute_error}

        if self.parameters.spin_polarized:
            max_spin = self.bands.n_spin
        else:
            max_spin = 1

        for spin in range(max_spin):
            spin_id = str("spin_"+str(spin))
            for qoi in self.ml.quantities_of_interest:
                for i, convergence_point in enumerate(self.convergence_points):
                    res = self.result_dict[spin_id][qoi]
                    if qoi == 'alphas':
                        plot_calculated_vs_predicted(res['true_array'].flatten(),
                                                     res['pred_array'][i, :, :].flatten(), qoi,
                                                     self.dirs[f'convergence_final_results_{i}'] /
                                                     f"{spin_id}_{qoi}_calculated_vs_predicted.png",
                                                     ('MAE', res['MAE']['mean'][i]))
                    elif qoi == 'evs':
                        plot_error_histogram(res['true_array'].flatten(), res['pred_array'][i, :, :].flatten(),
                                             qoi, self.dirs[f'convergence_final_results_{i}'] /
                                             f"{spin_id}_{qoi}_error_histogram.png",
                                             ('MAE', res['MAE']['mean'][i]))
                for metric in metrics:
                    res = self.result_dict[spin_id][qoi][metric]
                    plot_convergence(self.convergence_points, res['mean'], res['stdd'], qoi, metric, spin,
                                     self.dirs['convergence_final_results'] / f"{spin_id}_{qoi}_{metric}_convergence")
