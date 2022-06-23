import pickle
import matplotlib.pyplot as plt

with open('2x2x2/si_bandstructure.fig.pkl', 'rb') as fd:
    pickle.load(fd)

ax = plt.gca()
ax.set_ylim([-5, 5])
ax.set_ylabel(r'$\omega$ (eV)')
# plt.show() uncomment to view the figure interactively
plt.savefig('2x2x2/si_bandstructure_rescaled.png')
