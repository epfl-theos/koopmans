from koopmans.io import read

for grid_size in [2, 4, 8]:
    print('grid_size' in locals())
    print('grid_size' in globals())
    wf = read('si.json')
    wf.kgrid = [grid_size, grid_size, grid_size]
    wf.name = 'si_{0}x{0}x{0}'.format(grid_size)
    # wf.run()
