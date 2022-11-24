from pathlib import Path

import matplotlib
import numpy as np

from koopmans import io

matplotlib.use('Agg')  # nopep8
import matplotlib.pyplot as plt  # nopep8


def make_bar_diagrams(HOMOs_pred, HOMOs_calc):

    # Creating the figure
    _, ax = plt.subplots(figsize=(15.0, 7.5))
    y_1 = -HOMOs_pred
    y_2 = -HOMOs_calc
    x_axis = np.arange(len(y_1))
    x_labels = x_axis+11
    ax.bar(x_axis-0.2, y_1, 0.4, label='predicted')
    ax.bar(x_axis+0.2, y_2, 0.4, label='calculated')
    ax.set_ylabel('IP in eV')
    ax.set_xlabel('snapshot number')
    ax.set_xticks(x_axis, x_labels)
    ax.set_ylim(0, 15)
    ax.legend()
    plt.savefig('bar_diagram_predicted_and_calculated.png', facecolor=(1, 1, 1, 0))


if __name__ == '__main__':
    tutorial_folder = Path('tutorial_5b')
    num_train_snapshots = 10
    num_test_snapshots = 10
    num_occ_orbitals = 4

    # get all the HOMO energies
    true_results_folder = tutorial_folder / 'convergence_analysis' / 'true'
    pred_results_folder = tutorial_folder / 'convergence_analysis' / 'pred' / 'predicted_after_5'

    true_results = [np.loadtxt(true_results_folder /
                               f'evs_snapshot_{num_train_snapshots + i + 1}.txt') for i in range(num_test_snapshots)]
    pred_results = [np.loadtxt(pred_results_folder /
                               f'evs_snapshot_{num_train_snapshots + i + 1}.txt') for i in range(num_test_snapshots)]

    HOMOs_true = np.array([true_result[0, num_occ_orbitals-1] for true_result in true_results])

    HOMOs_pred = np.array([pred_result[0, num_occ_orbitals-1] for pred_result in pred_results])

    # plot the bar diagram of the HOMO energies (=-IPs)
    make_bar_diagrams(HOMOs_pred, HOMOs_true)
