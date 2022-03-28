from koopmans.io import read

# Load the workflow object
wf = read('si.kwf')

# Access the band structure from the very last calculation
results = wf.calculations[-1].results
bs = results['band structure']

# Print the band strucutre to file
bs.plot(filename='ki_bandstructure.png')

# Extract the band gap
n_occ = wf.projections.num_wann(occ=True)
gap = bs.energies[:, :, n_occ:].min() - bs.energies[:, :, :n_occ].max()
print(f'KI band gap = {gap:.2f} eV')
