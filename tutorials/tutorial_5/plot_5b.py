from pathlib import Path

import matplotlib
import numpy as np

from koopmans import io

matplotlib.use('Agg')  # nopep8
import matplotlib.pyplot as plt  # nopep8


def make_bar_diagrams(gaps_predicted, gaps_calculated):

    # Creating the figure
    _, ax = plt.subplots(figsize=(7.5, 7.5))

    y_1 = gaps_predicted
    y_2 = gaps_calculated
    x_axis = np.arange(len(y_1))
    x_labels = x_axis+1

    ax.bar(x_axis-0.2, y_1, 0.4, label='predicted')
    ax.bar(x_axis+0.2, y_2, 0.4, label='calculated')
    ax.set_ylabel('LUMO-HOMO in eV')
    ax.set_xlabel('snapshot number')
    ax.set_xticks(x_axis, x_labels)
    ax.set_ylim(0, 15)
    ax.legend()

    plt.savefig('bar_diagram_predicted_and_calculated.png', facecolor=(1, 1, 1, 0))


# if __name__ == '__main__':
#     # Read in the ml workflow
#     wf = io.read(Path('tutorial_5a') / 'h2o_trajectory_ml.kwf')

#     # get all the final calculations (one final calculation per snapshot)
#     final_calculations = [c for c in wf['calculations'] if c.prefix == 'ki_final']

#     # from the final calculations we can extract the HOMO and LUMO energies
#     homos_predicted = np.array([final_calculation.results['homo_energy'] for final_calculation in final_calculations])
#     lumos_predicted = np.array([final_calculation.results['lumo_energy'] for final_calculation in final_calculations])

#     # make_bar_diagrams(lumos_predicted-homos_predicted)

#     # Read in the ml-convergence workflow
#     wf = io.read(Path('tutorial_5a') / 'h2o_trajectory_ml.kwf')

#     # get all the final calculations (one final calculation per snapshot)
#     final_calculations = [c for c in wf['calculations'] if c.prefix == 'ki_final']

#     # from the final calculations we can extract the HOMO and LUMO energies
#     homos_calculated = np.array([final_calculation.results['homo_energy'] for final_calculation in final_calculations])
#     lumos_calculated = np.array([final_calculation.results['lumo_energy'] for final_calculation in final_calculations])

#     make_bar_diagrams(lumos_predicted-homos_predicted, lumos_calculated-homos_calculated)

if __name__ == '__main__':
    tutorial_folder = Path('tutorial_5b')
    num_train_snapshots = 10
    num_test_snapshots = 10
    num_occ_orbitals = 4

    true_results_folder = tutorial_folder / 'convergence_analysis' / 'true'
    pred_results_folder = tutorial_folder / 'convergence_analysis' / 'pred' / 'predicted_after_5'

    true_results = [np.loadtxt(true_results_folder /
                               f'evs_snapshot_{num_train_snapshots + i + 1}.txt') for i in range(num_test_snapshots)]
    pred_results = [np.loadtxt(pred_results_folder /
                               f'evs_snapshot_{num_train_snapshots + i + 1}.txt') for i in range(num_test_snapshots)]

    HOMOs_true = np.array([true_result[0, num_occ_orbitals-1] for true_result in true_results])
    LUMOs_true = np.array([true_result[0, num_occ_orbitals] for true_result in true_results])

    HOMOs_pred = np.array([pred_result[0, num_occ_orbitals-1] for pred_result in pred_results])
    LUMOs_pred = np.array([pred_result[0, num_occ_orbitals] for pred_result in pred_results])

    make_bar_diagrams(LUMOs_pred-HOMOs_pred, LUMOs_true-HOMOs_true)
