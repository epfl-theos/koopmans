from pathlib import Path

import matplotlib
import numpy as np

from koopmans import io

matplotlib.use('Agg')  # nopep8
import matplotlib.pyplot as plt  # nopep8


def make_bar_diagrams(gaps_predicted):

    # Creating the figure
    _, ax = plt.subplots(figsize=(7.5, 7.5))

    y_1 = gaps_predicted
    x_axis = np.arange(len(y_1))
    x_labels = x_axis+1

    ax.bar(x_axis-0.2, y_1, 0.4, label='predicted')
    # ax.bar(x_axis+0.2, y_2, 0.4, label='calculated')
    ax.set_ylabel('LUMO-HOMO in eV')
    ax.set_xlabel('snapshot number')
    ax.set_xticks(x_axis, x_labels)
    ax.set_ylim(0, 15)
    ax.legend()

    plt.savefig('bar_diagram_predicted.png', facecolor=(1, 1, 1, 0))


if __name__ == '__main__':
    num_train_snapshots = 10
    num_test_snapshots = 10

    # Read in the ml workflow
    wf = io.read(Path('tutorial_5a') / 'h2o_trajectory_ml.kwf')

    # get all the final calculations (one final calculation per snapshot)
    final_calculations = [c for c in wf['calculations'] if c.prefix == 'ki_final']

    # from the final calculations we can extract the HOMO and LUMO energies
    homos_predicted = np.array([final_calculation.results['homo_energy']
                               for final_calculation in final_calculations[num_train_snapshots:num_train_snapshots+num_test_snapshots]])
    lumos_predicted = np.array([final_calculation.results['lumo_energy']
                               for final_calculation in final_calculations[num_train_snapshots:num_train_snapshots+num_test_snapshots]])

    make_bar_diagrams(lumos_predicted-homos_predicted)
