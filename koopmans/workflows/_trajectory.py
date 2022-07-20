

from pathlib import Path
from ase.calculators.calculator import Calculator
from sklearn.metrics import mean_absolute_error as mae
from ._workflow import Workflow
from typing import List, Dict, Any
from ase import Atoms, io
import json as json_ext
from koopmans import utils
from koopmans.ml_utils import RidgeRegression
import numpy as np
from koopmans import calculators
import os
import matplotlib.pyplot as plt

load_results_from_output = True


class TrajectoryWorkflow(Workflow):

    def __init__(self, *args, **kwargs):
        self.snapshots: List[Atoms] = kwargs.pop('snapshots', [])
        super().__init__(*args, **kwargs)

    def _run(self) -> None:
        """
        Runs the trajectory workflow.
        """

        # Initialize the RidgeRegression() model
        if self.parameters.use_ml:
            self.ml_model = RidgeRegression()

        if self.parameters.use_ml and self.parameters.mode == 'convergence':
            if self.parameters.number_of_training_snapshots >= len(self.snapshots):
                raise ValueError(
                    "There are not enough test-snapshots available. Please increase the number of snapshots in the xyz-file or decrease 'number_of_training_snapshots'")
            convergence_points = list(range(1, self.parameters.number_of_training_snapshots+1))
            convergence_dir = Path.cwd() / 'convergence_analysis'
            utils.system_call(f'rm -rf {convergence_dir}')
            convergence_dir_true = convergence_dir / "true"
            convergence_dir_pred = convergence_dir / "predicted"
            utils.system_call(f'mkdir -p {convergence_dir_true}')
            utils.system_call(f'mkdir -p {convergence_dir_pred}')

            test_indices = list(range(self.parameters.number_of_training_snapshots, len(self.snapshots)))
            tmp_number_of_training_snapshots = self.parameters.number_of_training_snapshots
            # set the number of training snapshots to a very high value to not use the prediction for the first run
            self.parameters.number_of_training_snapshots = 10000
            self.print(f'Obtaining ab-initio results for the last {len(test_indices)} snapshot(s)', style='heading')
            # get the ab-initio result for the test_indices
            self.run_trajectory(test_indices, get_evs=True, folder=convergence_dir_true)
            # set the number of training snapshots back to its original value
            self.parameters.number_of_training_snapshots = tmp_number_of_training_snapshots
            for convergence_point in convergence_points:
                convergence_dir_pred_curr = convergence_dir_pred / ("predicted_after_" + str(convergence_point))
                utils.system_call(f'mkdir -p {convergence_dir_pred_curr}')
                self.ml_model = RidgeRegression()  # reset the ML-model
                train_indices = list(range(0, convergence_point))  # get the the indices on which we train on
                self.print(
                    f'Training on {len(train_indices)} snapshot(s) and then testing on the last {len(test_indices)} snapshot(s)', style='heading')
                self.run_trajectory(train_indices)  # train the model
                # test the model (without retraining the model)
                self.run_trajectory(test_indices, get_evs=True, folder=convergence_dir_pred_curr)
            self.analyze_convergence_results(convergence_dir, convergence_dir_true,
                                             convergence_dir_pred, convergence_points)
        else:
            snapshot_indices = list(range(0, len(self.snapshots)))
            self.run_trajectory(snapshot_indices)

    def run_trajectory(self, indices: List[int], get_evs: bool = False, folder: Path() = Path.cwd()):
        # Import it like this so if they have been monkey-patched, we will get the monkey-patched version
        from koopmans.workflows import KoopmansDFPTWorkflow, KoopmansDSCFWorkflow

        for i in indices:
            snapshot = self.snapshots[i]
            self.parameters.current_snapshot = i
            self.atoms.set_positions(snapshot.positions)

            self.print(f'Performing Koopmans calculation on snapshot {i+1} / {len(self.snapshots)}', style='heading')

            if self.parameters.method == 'dfpt':
                workflow = KoopmansDFPTWorkflow(**self.wf_kwargs)
                raise NotImplementedError("TODO Yannick: What is the final_calculator of dftp")
            elif self.parameters.method == 'dscf':
                workflow = KoopmansDSCFWorkflow(**self.wf_kwargs)
                final_calculator = calculators.KoopmansCPCalculator
            else:
                raise NotImplementedError("The trajectory workflow is currently only implemented with dfpt and dscf.")
            # reset the bands to the initial guesses (i.e. either from file or to 0.6 but not from the previous calculation)
            self.bands = workflow.bands
            subdirectory = 'snapshot_' + str(i+1)
            self.run_subworkflow(workflow, subdirectory=subdirectory)

            if get_evs:
                final_calc = [c for c in workflow.calculations if isinstance(c, final_calculator)][-1]
                evs = final_calc.results['eigenvalues']
                alphas = self.bands.alphas
                file_evs = folder / ("evs" + "_" + subdirectory + ".txt")
                file_alphas = folder / ("alphas" + "_" + subdirectory + ".txt")
                np.savetxt(file_evs, evs)
                np.savetxt(file_alphas, alphas)

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
        wf = super(TrajectoryWorkflow, cls)._fromjsondct(bigdct)
        wf.snapshots = snapshots
        return wf

    def analyze_convergence_results(self, convergence_dir: Path, convergence_dir_true: Path, convergence_dir_pred: Path, convergence_points: List[int]):

        for spin in range(self.bands.n_spin):
            MAEs_alpha = np.zeros(len(convergence_points))
            stddevs = np.zeros(len(convergence_points))
            test_indices = list(range(self.parameters.number_of_training_snapshots, len(self.snapshots)))
            evs_true = np.zeros((len(test_indices), self.bands.n_bands[spin]))
            alphas_true = np.zeros((len(test_indices), self.bands.n_bands[spin]))
            for j, test_index in enumerate(test_indices):
                evs_true[j, :] = np.loadtxt(convergence_dir_true / f"evs_snapshot_{test_index+1}.txt")[spin, :]
                alphas_true[j, :] = np.loadtxt(convergence_dir_true / f"alphas_snapshot_{test_index+1}.txt")[spin, :]
            for i, convergence_point in enumerate(convergence_points):
                evs_pred = np.zeros((len(test_indices), self.bands.n_bands[spin]))
                alphas_pred = np.zeros((len(test_indices), self.bands.n_bands[spin]))
                MAE_alpha = []
                for j, test_index in enumerate(test_indices):
                    evs_pred[j, :] = np.loadtxt(
                        convergence_dir_pred / f"predicted_after_{convergence_point}" / f"evs_snapshot_{test_index+1}.txt")[spin, :]
                    alphas_pred[j, :] = np.loadtxt(
                        convergence_dir_pred / f"predicted_after_{convergence_point}" / f"alphas_snapshot_{test_index+1}.txt")[spin, :]
                    MAE_alpha.append(mae(alphas_true[j, :], alphas_pred[j, :]))
                MAEs_alpha[i] = np.mean(MAE_alpha)
                stddevs[i] = np.std(MAE_alpha)
                self.plot_evs_error_histogram(evs_true.reshape(-1), evs_pred.reshape(-j), MAEs_alpha[i],
                                              convergence_dir_pred / f"predicted_after_{convergence_point}")
                self.plot_alphas_calculated_vs_predicted(alphas_true.reshape(-1), alphas_pred.reshape(-1), MAEs_alpha[i],
                                                         convergence_dir_pred / f"predicted_after_{convergence_point}")
            self.plot_convergence(convergence_points, MAEs_alpha, stddevs, convergence_dir, spin)

    def plot_alphas_calculated_vs_predicted(self, y, y_pred, MAE, result_folder):
        lb_x = 0.995*min(y)
        ub_x = 1.005*max(y)
        x = np.linspace(lb_x, ub_x, 1000)
        fig, ax = plt.subplots(1, 1, figsize=(7.5, 7.5))
        ax.set_xlim((lb_x, ub_x))
        ax.set_ylim((lb_x, ub_x))
        ax.plot(y, y_pred, 'o', color='green', label="MAE = " + str(round(MAE, 15)))
        ax.plot(x, x, color='grey', linewidth=0.75)
        plt.ylabel(r"predicted $\alpha$")
        plt.xlabel(r"calculated $\alpha$")
        plt.legend()
        plt.savefig(result_folder / 'alphas_calculated_vs_predicted.png')
        plt.close()

    def plot_evs_calculated_vs_predicted(self, y, y_pred, MAE, result_folder):
        lb_x = 0.995*min(y)  # -17.5#
        ub_x = 1.005*max(y)  # -7.5#
        x = np.linspace(lb_x, ub_x, 1000)
        fig, ax = plt.subplots(1, 1, figsize=(7.5, 7.5))
        ax.set_xlim((lb_x, ub_x))
        ax.set_ylim((lb_x, ub_x))
        ax.plot(y, y_pred, 'o', color='green', label="estimated MAE = " + str(round(MAE, 5)))
        ax.plot(x, x, color='grey', linewidth=0.75)
        plt.ylabel(r"predicted ev in eV")
        plt.xlabel(r"calculated ev in eV")
        plt.legend()
        plt.savefig(result_folder / 'evs_calculated_vs_predicted.png')
        plt.close()

    def plot_evs_error_histogram(self, y, y_pred, MAE, result_folder):
        error = np.abs(y-y_pred)
        lb_x = 0.0  # -17.5#
        ub_x = 0.7  # -7.5#
        bins = np.linspace(lb_x, ub_x, 100)
        fig, ax = plt.subplots(1, 1, figsize=(7.5, 7.5))
        # ax.set_xlim((lb_x, ub_x))
        # ax.set_ylim((0, 1000))

        ax.hist(error, bins, label="MAE = " + str(round(MAE, 15)))
        plt.ylabel(r"number of eigenvalues")  # $\alpha$")
        plt.xlabel(r"error in eV")  # $\alpha$")
        plt.legend()
        plt.savefig(result_folder / 'error_histogram_evs.png')
        plt.close()

    def plot_convergence(self, convergence_points, MAEs, stddevs, convergence_folder, spin):
        def snap2orb(x):
            return x * self.bands.n_bands[spin]

        def orb2snap(x):
            return x/self.bands.n_bands[spin]

        fig, ax = plt.subplots(1, 1, figsize=(15.0, 7.5))
        ax.xaxis.get_major_locator().set_params(integer=True)
        ax.set_xlabel("Number of snapshots for training")
        ax.set_ylabel(r"MAE of predicted $\alpha$")
        ax.set_ylim(0.0, 0.0175)
        ax.errorbar(convergence_points, MAEs, stddevs, marker='o', linestyle="--")
        secax = ax.secondary_xaxis('top', functions=(snap2orb, orb2snap))
        secax.set_xlabel('Number of orbitals for training')
        plt.savefig(convergence_folder / 'convergence_plot.png')
        plt.close()
