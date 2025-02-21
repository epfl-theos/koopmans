import matplotlib.pyplot as plt
import numpy as np

from koopmans import io

# Load the workflow python object from the .pkl file
wf = io.read('h2o_test.pkl')

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
ax0.plot(true_orbital_energies, predicted_orbital_energies, 'o', alpha=0.5, markeredgewidth=0)
xylimits = ax0.get_xlim() + ax0.get_ylim()
limits = [min(xylimits), max(xylimits)]
ax0.plot(limits, limits, color='k', linestyle='--', zorder=-1, lw=1)
ax0.set_xlim(limits)
ax0.set_ylim(limits)
ax0.set_aspect('equal', adjustable='box')
ax0.set_xlabel(r'$\varepsilon^\mathrm{true}_i$ (eV)')
ax0.set_ylabel(r'$\varepsilon^\mathrm{predicted}_i$ (eV)')

# Histogram
errors_meV = errors * 1000
ax1.hist(errors_meV, bins=20, orientation='vertical')
ax1.axvline(0, color='black', linestyle='--', lw=1)
ax1.set_yticks([])
ax1.set_xlabel(r'$\varepsilon^\mathrm{predicted}_i - \varepsilon^\mathrm{true}_i$ (meV)')
mean = np.mean(errors_meV)
ax1.axvline(mean, color='grey', linestyle='-', lw=1)
std = np.std(errors_meV)
ymax = ax1.get_ylim()[1]
ax1.fill_betweenx([0, ymax], mean - std, mean + std, color='gray', alpha=0.5, lw=0)
ax1.set_ylim(0, ymax)
ax1.text(0.98, 0.98, f'mean = {mean:.0f} meV\n$\sigma = {std:.0f}$ meV', transform=ax1.transAxes, ha='right', va='top')

# Save the figure
plt.subplots_adjust(left=0.12, right=0.99, top=0.98, bottom=0.18, wspace=0.07)
plt.savefig('testing.png', dpi=300)
