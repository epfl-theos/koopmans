import shutil

from koopmans.io import read
from koopmans.utils import chdir

for grid_size in [2, 4, 8]:
    # Read the input file
    wf = read('si.json')

    # Modify the kgrid
    wf.kpoints.grid = [grid_size, grid_size, grid_size]

    # Run the workflow
    wf.run(subdirectory='{0}x{0}x{0}'.format(grid_size))
