import matplotlib
import pandas as pd

from koopmans import io
from koopmans.cell import cell_to_parameters

matplotlib.use('Agg')  # nopep8
import matplotlib.pyplot as plt  # nopep8

# Read in the workflow
wf = io.read('h2o_conv.kwf')
calcs = wf.calculations

# Create a simple function for extracting celldm1 from a calculator


def get_celldm1(calc):
    return cell_to_parameters(calc.atoms.cell)['celldms'][1]


# Extract the list of ecutwfcs and cell sizes tested
ecutwfcs = sorted(list(set([c.parameters.ecutwfc for c in calcs])))
cell_sizes = sorted(list(set([get_celldm1(c) for c in calcs])))

# Extract the homo energies and store them in a pandas dataframe
df = pd.DataFrame(columns=cell_sizes, index=ecutwfcs)
for c in calcs:
    df.loc[c.parameters.ecutwfc, get_celldm1(c)] = c.results['homo_energy']

# Convert the homo energies to the absolute difference from the most accurate value
df.loc[:, :] -= df.loc[max(ecutwfcs), max(cell_sizes)]
df = df.abs()

# Plot the data (with a small tweak to make the labels more readable)
ax = df.rename(columns={c: f'{c:.1f}' for c in df.columns}).plot(logy=True, marker='o')

# Plot a horizontal line for the convergence threshold
ax.axhline(0.01, c='grey', ls='--')

# Figure aesthetics
ax.set_ylabel(r'$|\Delta\varepsilon_{HOMO}|$ (eV)')
ax.set_xlabel('energy cutoff (Ha)')
ax.legend(title='celldm1 ($\AA$)', ncol=2)
plt.tight_layout()

# Save the figure
plt.savefig('convergence.png', facecolor=(1, 1, 1, 0))
