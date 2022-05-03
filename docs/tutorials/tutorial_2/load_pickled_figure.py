import pickle
import matplotlib.pyplot as plt

with open('interpolated_bandstructure_2x2x2.fig.pkl', 'rb') as fd:
    pickle.load(fd)

ax = plt.gca()
ax.set_ylim([-5, 5])
ax.set_ylabel(r'$\omega$ (eV)')
# plt.show() uncomment to view the figure interactively
plt.savefig('interpolated_bandstructure_2x2x2_modified.png')
