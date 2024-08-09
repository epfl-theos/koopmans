from pathlib import Path

import matplotlib

matplotlib.use('Agg')  # nopep8
import matplotlib.pyplot as plt  # nopep8
import numpy as np

from koopmans import io

if __name__ == '__main__':
    # Extract the eigenvalue information from the .kwf file
    wf = io.read('tutorial_5a/predict/h2o_predict.kwf')
    calcs = [c for c in wf.calculations if c.prefix == 'ki_final']
    eigenvalues = np.array([c.results['eigenvalues'][0] for c in calcs])

    # Create a violin plot figure
    fig, ax = plt.subplots(figsize=(5, 3))
    plt.violinplot(eigenvalues, vert=False, showmeans=True)

    # Figure aesthetics
    n_orbs_occ = 4
    n_orbs = 6
    labels = [f'HOMO - {n_orbs_occ - i}' for i in range(1, n_orbs_occ)] + ['HOMO',
                                                                           'LUMO'] + [f'LUMO + {i - n_orbs_occ}' for i in range(n_orbs_occ + 1, n_orbs)]
    ax.set_yticks(range(1, n_orbs + 1), labels)
    ax.set_xlabel('energy (eV)')
    plt.tight_layout()

    # Save the figure
    plt.savefig('predicting.png', facecolor=(1, 1, 1, 0), dpi=300)
