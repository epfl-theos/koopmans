'''
Script for running environ PBE DSCF tutorial

Written by Edward Linscott, Jan 2021
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from koopmans import io

sns.set_style('darkgrid')
sns.set_style("darkgrid", {"axes.facecolor": "d2d2d9"})

if __name__ == '__main__':
    # Read in json file
    workflow = io.read_json('o2_environ_dscf.json')

    # Run workflow
    workflow.run()

    # Load results into a pandas dataframe
    epsilons = sorted(list(set([c.environ_settings['ENVIRON']['env_static_permittivity'] for c in
                                workflow.calculations])))
    results = set([k for c in workflow.calculations for k in c.results.keys()])

    columns = [f'{label} {result}' for label in ['charged', 'neutral'] for result in results]
    columns.append('EA')
    df = pd.DataFrame(index=epsilons, columns=columns)
    for calc in workflow.calculations:
        epsilon = float(calc.environ_settings['ENVIRON']['env_static_permittivity'])
        if calc.parameters.tot_charge == 0:
            label = 'neutral'
        else:
            label = 'charged'
        df.loc[epsilon][f'{label} calc'] = calc
        for key, val in calc.results.items():
            if val is not None:
                df.loc[epsilon][f'{label} {key}'] = val

    # Calculate the electron affinities via Nattino et al equation 11 (http://dx.doi.org/10.1021/acs.jctc.9b00552)
    df['EA'] = df.apply(lambda row: row['neutral energy'] - row['neutral electrostatic embedding']
                        - row['charged energy'] + row['charged electrostatic embedding'], axis=1)

    # Remove the rows without EA values
    df = df.dropna()

    # Plot the results
    ax = df.plot(y='EA', style='o', label='data')

    # Interpolate with a quartic
    coeffs = np.polyfit(df.index, df['EA'], deg=4)
    x_interp = np.linspace(1, max(epsilons))
    y_interp = np.polyval(coeffs, x_interp)
    ax.plot(x_interp, y_interp, label='quartic fit')

    # Label the extrapolated value
    ax.annotate(f'EA = {y_interp[0]:.3f} eV', xy=(1, y_interp[0]), xycoords='data', xytext=(1.5, y_interp[0]),
                ha='left', va='center', arrowprops={'arrowstyle': '->', 'color': 'k'})

    # Figure aesthetics
    ax.set_xlabel(r'$\varepsilon$')
    ax.set_ylabel(r"$\Delta E'$ (eV)")
    ax.set_xlim([1, max(epsilons)])
    ax.legend()
    plt.tight_layout()

    # Save the figure to file
    plt.savefig('o2_dscf_ea_result.png', facecolor=(1, 1, 1, 0))
