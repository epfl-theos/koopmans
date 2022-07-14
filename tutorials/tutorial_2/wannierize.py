from koopmans.io import read
from koopmans.utils import chdir

for grid_size in [2, 4, 8]:

    # Read in the input file
    wf = read('si.json')

    # Modify the kgrid
    wf.kgrid = [grid_size, grid_size, grid_size]

    # Run the workflow in a subdirectory
    with chdir('{0}x{0}x{0}'.format(grid_size)):
        wf.run()
