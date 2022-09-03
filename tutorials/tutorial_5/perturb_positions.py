import numpy as np

from ase import Atoms, io

snapshots = io.read("h2o_unpertubed_positions.xyz", index='0')
print(snapshots.positions)
