from koopmans.io import read

for grid_size in [2, 4, 8]:

    print('\n {0}x{0}x{0} grid'.format(grid_size))

    wf = read('si.json', override={'setup': {'k_points': {'kgrid': [grid_size for _ in range(3)]}}})
    wf.run()
