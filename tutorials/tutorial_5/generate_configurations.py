"""Generate input configurations of water molecules."""

import copy

import numpy as np
from ase_koopmans import io

np.random.seed(0)

cell = np.diag([6.892900, 6.892900, 6.892900])
unperturbed_h2o = io.read("h2o.xyz")
unperturbed_h2o.set_cell(cell)
unperturbed_h2o.pbc = True

trajectory = []
for i in range(20):
    random_pertubation = np.random.normal(0, 1, (3, 3))
    new_perturbed_h2o = copy.deepcopy(unperturbed_h2o)
    new_perturbed_h2o.positions += 0.1 * random_pertubation
    trajectory.append(new_perturbed_h2o)

io.write('snapshots.gif', trajectory, interval=500, rotation='-100y,-10x')
io.write('01-train/training_snapshots.xyz', trajectory[:5])
io.write('02-predict/predicting_snapshots.xyz', trajectory[5:])
io.write('03-test/testing_snapshots.xyz', trajectory[5:])
io.write('04-advanced-testing/snapshots.xyz', trajectory)
