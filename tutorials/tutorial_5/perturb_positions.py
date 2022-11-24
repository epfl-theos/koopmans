import copy

import numpy as np

from ase import io

np.random.seed(0)

cell = np.diag([6.892900, 6.892900, 6.892900])
unperturbed_h2o = io.read("h2o_unpertubed_positions.xyz", index='0')
unperturbed_h2o.set_cell(cell)

trajectory = []
for i in range(20):
    random_pertubation = np.random.normal(0, 1, (3, 3))
    new_perturbed_h2o = copy.deepcopy(unperturbed_h2o)
    new_perturbed_h2o.positions += 0.1*random_pertubation
    trajectory.append(new_perturbed_h2o)

io.write('snapshots.gif', trajectory, interval=500, rotation='-100y,-10x')
io.write('tutorial_5a/snapshots.xyz', trajectory)
io.write('tutorial_5b/snapshots.xyz', trajectory)
