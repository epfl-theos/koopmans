"""Plot the band structore of ZnO, modifying the pre-saved plot file."""

import pickle

import matplotlib.pyplot as plt

# Use the python library "pickle" to load the *.fig.pkl file
fig = pickle.load(open('zno_bandstructure.fig.pkl', 'rb'))

# Rescale the y axes
fig.axes[0].set_ylim([-5, 15])

# Show/save the figure (uncomment as desired)
plt.savefig('zno_bandstructure_rescaled.png')
# plt.show()
