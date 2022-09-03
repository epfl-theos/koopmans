import numpy as np
import copy
from ase import io

np.random.seed(0)
cell = np.diag([6.892900, 6.892900, 6.892900])

unperturbed_h2o = io.read("h2o_unpertubed_positions.xyz", index='0')
unperturbed_h2o.set_cell(cell)


trajectory = []

for i in range(10):
    random_pertubation = np.random.rand(3, 3)-0.5
    new_perturbed_h2o = copy.deepcopy(unperturbed_h2o)
    new_perturbed_h2o.positions += 0.5*random_pertubation
    trajectory.append(new_perturbed_h2o)

io.write('snapshots.gif', trajectory, interval=500, rotation='-100y,-10x')
io.write('snapshots.xyz', trajectory)
