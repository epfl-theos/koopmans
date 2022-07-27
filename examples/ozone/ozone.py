'''
Script for running the ozone example

Written by Edward Linscott, Jun 2021
'''

import matplotlib.pyplot as plt
import numpy as np

import koopmans.mpl_config
from koopmans import io, utils


def run(from_scratch=False):
    names = ['dscf', 'dfpt', 'dfpt_with_dscf_screening']

    wfs = {}
    for name in names:
        with utils.chdir(name):
            wf = io.read('ozone.json')
            wf.parameters.from_scratch = from_scratch
            if 'screening' in name:
                wf.parameters.calculate_alpha = False
                wf.parameters.alpha_guess = wfs['dscf'].bands.alphas
            wf.run()

            # Save workflow to file
            io.write(wf, 'ozone.kwf')

        wfs[name] = wf
    return wfs


def plot(wfs):
    # Setting up figure
    _, axes = plt.subplots(ncols=1, nrows=3, sharex=True)
    x = np.array(wfs['dscf'].bands.indices)

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
    alphas = np.array([wf.bands.alphas for wf in wfs.values()])
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
