import copy

import matplotlib.pyplot as plt
import numpy as np

from koopmans import io, utils
from koopmans.ml import MLModel
from koopmans.workflows import InitializationWorkflow, KoopmansDSCFWorkflow

# Load the workflow python object from the .kwf file
wf = io.read('tutorial_5a/test/h2o_test.kwf')

# Select the calculations that correspond to the final KI calculations using predicted/ab initio screening parameters
calculations_with_predicted_alpha = [c for c in wf.calculations if c.prefix == 'ki_final_ml']
calculations_with_true_alpha = [c for c in wf.calculations if c.prefix == 'ki_final']

# Extract the orbital energies for these calculations and compute the difference
predicted_orbital_energies = [e for c in calculations_with_predicted_alpha for e in c.results['eigenvalues'][0]]
true_orbital_energies = [e for c in calculations_with_true_alpha for e in c.results['eigenvalues'][0]]
errors = np.array(predicted_orbital_energies) - np.array(true_orbital_energies)

# Plotting
fig, [ax0, ax1] = plt.subplots(1, 2, figsize=(6, 3))

# Scatterplot
ax0.plot(true_orbital_energies, predicted_orbital_energies, 'o')
xylimits = ax0.get_xlim() + ax0.get_ylim()
limits = [min(xylimits), max(xylimits)]
ax0.plot(limits, limits, color='k', linestyle='--', zorder=-1, lw=1)
ax0.set_xlim(limits)
ax0.set_ylim(limits)
ax0.set_aspect('equal', adjustable='box')
ax0.set_xlabel(r'$\varepsilon^\mathrm{true}_i$ (eV)')
ax0.set_ylabel(r'$\varepsilon^\mathrm{predicted}_i$ (eV)')

# Histogram
ax1.hist(errors, bins=20, orientation='vertical')
ax1.axvline(0, color='black', linestyle='--', lw=1)
ax1.set_yticks([])
ax1.set_xlabel(r'$\varepsilon^\mathrm{predicted}_i - \varepsilon^\mathrm{true}_i$ (eV)')

# Save the figure
plt.subplots_adjust(left=0.12, right=0.99, top=0.98, bottom=0.18, wspace=0.07)
plt.savefig('testing.png', dpi=300)
