

'''
Script for running the ozone example
'''

import numpy as np

# isort: off
import koopmans.mpl_config
import matplotlib.pyplot as plt
# isort: on

from ase.build import molecule
from koopmans import io, utils, workflows


def run(from_scratch=False):
    names = ['dscf', 'dfpt', 'dfpt_with_dscf_screening']
    ozone = molecule('O3', vacuum=6)
    parameters = {'functional': 'ki', 'init_orbitals': 'kohn-sham', 'frozen_orbitals': True}

    wfs = {}
    for name in names:
        with utils.chdir(name):
            wf = workflows.SinglepointWorkflow(atoms=ozone, parameters=parameters,
                                               nbnd=10, method=name[:4], ecutwfc=65.0)
            wf.parameters.from_scratch = from_scratch
            if 'screening' in name:
                wf.parameters.calculate_alpha = False
                wf.parameters.alpha_guess = [wfs['dscf'].bands.alphas[0]]

            # io.write(wf, wf.name + '.json')
            wf.run()
        wfs[name] = wf
    return wfs


def plot(wfs):
    # Setting up figure
    _, axes = plt.subplots(ncols=1, nrows=3, sharex=True)
    x = np.array(wfs['dscf'].bands.indices[0])

    # Plotting self-Hartrees
    ax = axes[0]
    dscf_sh = wfs['dscf'].calculations[-1].results['orbital_data']['self-Hartree'][0]
    dfpt_sh = wfs['dfpt'].calculations[-2].results['orbital_data']['self-Hartree']
    bar1 = ax.bar(x - 0.125, dfpt_sh, width=0.25, label='DFPT')
    bar2 = ax.bar(x + 0.125, dscf_sh, width=0.25, label=r'$\Delta$SCF')
    ax.set_ylabel('self-Hartree (eV)')
    ax.legend(ncol=2, bbox_to_anchor=(1, 1), loc='upper right', frameon=False)

    # Plotting alphas
    ax = axes[1]
    alphas = np.array([wf.bands.alphas[0] for wf in wfs.values()])
    ax.bar(x, alphas[1] - alphas[0], width=0.25, label='DFPT')
    ax.set_ylabel(r'$\alpha_i^\mathrm{DFPT} - \alpha_i^\mathrm{\Delta SCF}$')
    ax.axhline(y=0, color='k', linestyle='-', lw=0.5)

    # Plotting eigenvalues
    ax = axes[2]
    ax.axhline(y=0, color='k', linestyle='-', lw=0.5)
    eigs = np.array([wf.calculations[-1].results['eigenvalues'][0] for wf in wfs.values()])
    bar1 = ax.bar(x - 0.125, eigs[1] - eigs[0], width=0.25, label='DFPT')
    bar2 = ax.bar(x + 0.125, eigs[2] - eigs[0], width=0.25, label=r'DFPT@$\alpha_{\Delta SCF}$')
    label_bars(bar1 + bar2)
    ax.set_ylabel(r'$\varepsilon_i - \varepsilon_i^\mathrm{\Delta SCF}$')
    ax.set_xlabel('orbital index')
    ax.set_xticks(x)
    ax.legend(ncol=2, bbox_to_anchor=(0, 1), loc='upper left', frameon=False)

    # Saving figure and closing
    plt.tight_layout()
    plt.savefig('comparison.pdf', format='pdf')
    plt.close('all')


def label_bars(bars):
    # For adding text labels to a bar plot
    for rect in bars:
        height = rect.get_height()
        if height < 0:
            va = 'top'
        else:
            va = 'bottom'
        plt.text(rect.get_x() + rect.get_width() / 2.0, height,
                 '{:.3f}'.format(height), ha='center', va=va, fontsize=3, color=rect.get_facecolor())


if __name__ == '__main__':
    # Run calculations
    wfs = run(from_scratch=True)

    # Plot the results
    plot(wfs)
