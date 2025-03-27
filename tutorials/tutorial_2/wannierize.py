"""Wannierize the band structure of bulk silicon for three different k-point grids."""

from koopmans.io import read

for grid_size in [2, 4, 8]:
    # Read the input file
    wf = read('si.json')

    # Modify the kgrid
    wf.kpoints.grid = [grid_size, grid_size, grid_size]

    # Run the workflow
    wf.directory = '{0}x{0}x{0}'.format(grid_size)
    wf.run()
