import matplotlib

from koopmans import io

matplotlib.use('Agg')  # nopep8
import matplotlib.pyplot as plt  # nopep8

# Read in the workflow
wf = io.read('h2o_conv.kwf')
calcs = wf.calculations

# Adding cell_size attribute to make things easier
for calc in calcs:
    calc.cell_size = calc.atoms.cell[0, 0]

# Extracting the list of ecutwfcs and cell sizes tested
ecutwfcs = sorted(list(set([c.parameters.ecutwfc for c in calcs])))
cell_sizes = sorted(list(set([c.cell_size for c in calcs])))

# Creating the figure
_, ax = plt.subplots()

# Extract the reference HOMO level from the most accurate calculation
[reference_calc] = [c for c in calcs if c.parameters.ecutwfc == max(ecutwfcs) and c.cell_size == max(cell_sizes)]
ref_homo = reference_calc.results['homo_energy']

for cell_size in cell_sizes:
    # Plot a line of delta e_HOMO for this particular cell size
    selected_calcs = sorted([c for c in calcs if c.cell_size == cell_size], key=lambda x: x.parameters.ecutwfc)
    assert [c.parameters.ecutwfc for c in selected_calcs] == ecutwfcs
    x = ecutwfcs
    y = [abs(c.results['homo_energy'] - ref_homo) for c in selected_calcs]
    if cell_size == max(cell_sizes):
        x = x[:-1]
        y = y[:-1]
    cell_size /= min(cell_sizes)
    plt.semilogy(x, y, 'o-', label=f'{cell_size:.1f}')

# Plot the convergence threshold
ax.axhline(0.01, c='grey', ls='--')

# Figure aesthetics
ax.set_ylabel(r'$|\Delta\varepsilon_{HOMO}|$ (eV)')
ax.set_xlabel('energy cutoff (Ha)')
ax.legend(title='cell size', ncol=2)
plt.tight_layout()

# Save the figure
plt.savefig('convergence.png', facecolor=(1, 1, 1, 0))
