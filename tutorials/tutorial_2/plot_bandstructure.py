"""Plot the bandstructure of bulk silicon and extract the band gap."""

from koopmans.io import read

# Load the workflow object
wf = read('si.pkl')

# Access the band structure from the outputs of the workflow
bs = wf.outputs.band_structures['ki']

# Shift the band structure so that the VBM is at 0 eV
bs = bs.subtract_reference()

# Print the band structure to file
bs.plot(filename='si_bandstructure.png')

# Extract the band gap
n_occ = wf.number_of_electrons() // 2
gap = bs.energies[:, :, n_occ:].min() - bs.energies[:, :, :n_occ].max()
print(f'KI band gap = {gap:.2f} eV')
