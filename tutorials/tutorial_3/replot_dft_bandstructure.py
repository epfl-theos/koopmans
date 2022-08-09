import pickle

import matplotlib.pyplot as plt

# Use the python library "pickle" to load the *.fig.pkl file
fig = pickle.load(open('dft_bands/zno_dftbands_bandstructure.fig.pkl', 'rb'))

# Rescale the y axes
fig.axes[0].set_ylim([-5, 15])

# Show the figure
plt.show()
