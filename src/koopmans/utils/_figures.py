import pickle

# isort: off
from koopmans import mpl_config
import matplotlib.pyplot as plt
# isort: on


def savefig(fname, **kwargs):
    # Saves a figure twice
    #  1) as a png (or other) using matplotlib's plt.savefig()
    #  2) in editable form using pickle's dump()
    splitfname = fname.rsplit('.', 1)
    fname = splitfname[0]
    if len(splitfname) == 1:
        fmt = 'png'
    else:
        fmt = splitfname[1]

    plt.savefig(fname=f'{fname}.{fmt}', dpi=1000, **kwargs)

    with open(fname + '.fig.pkl', 'wb') as fd:
        pickle.dump(plt.gcf(), fd)
