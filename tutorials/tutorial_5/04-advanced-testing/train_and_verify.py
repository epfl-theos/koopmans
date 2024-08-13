import copy

import matplotlib.pyplot as plt
import numpy as np

from koopmans import io
from koopmans.ml import MLModel
from koopmans.workflows import InitializationWorkflow, KoopmansDSCFWorkflow

# Optional: reduce n_total and n_train to make the script run faster
n_total = 20
n_train = 10
n_test = n_total - n_train

# # Step 1: do ab initio calculations on all the snapshots
ab_initio_wfs = []
for i in range(n_total):
    # Create the workflow from the template JSON input file
    ab_initio_wf = io.read('h2o_singlepoint.json', override={'ml': {'train': True}})

    # Set the atomic positions for the current snapshot
    ab_initio_wf.atoms = ab_initio_wf.snapshots[i]
    ab_initio_wf.parameters.from_scratch = False

    # Run the workflow
    ab_initio_wf.run(subdirectory=f'ab_initio/snapshot_{i + 1}')

    # Store the result
    ab_initio_wfs.append(ab_initio_wf)

# Step 2: train an ML model on the first x configurations, for x in (1, ..., n_train)
models = []
training_data = []
for i, ab_initio_wf in enumerate(ab_initio_wfs[:n_train]):
    model = MLModel('ridge_regression', descriptor='orbital_density')
    training_data += ab_initio_wf.bands.to_solve
    model.add_training_data(training_data)
    model.train()
    models.append(model)

# Step 3: test each ML model on the last n_test configurations
n_eigs = len(ab_initio_wfs[0].calculations[-1].results['eigenvalues'][0])
errors = [[] for _ in range(n_train)]
for i_snapshot in range(n_train, n_total):
    ab_initio_wf = ab_initio_wfs[i_snapshot]
    eigenvalue_errors = []

    # Run the initialization workflow (which we don't want to repeat for each model)
    init_wf = InitializationWorkflow.fromparent(ab_initio_wf)
    init_wf.bands = copy.deepcopy(ab_initio_wf.bands)
    init_wf.parent = None
    init_wf.ml.train = False
    init_wf.run(subdirectory=f'ml/snapshot_{i_snapshot + 1}/init')
    descriptors = [b.power_spectrum for b in init_wf.bands.to_solve]

    for i_train, model in enumerate(models):
        # Create a Koopmans DSCF workflow, linking the outputs of the initialization workflow
        ml_wf = KoopmansDSCFWorkflow.fromparent(ab_initio_wf,
                                                initial_variational_orbital_files=init_wf.outputs.variational_orbital_files,
                                                previous_cp_calc=init_wf.outputs.final_calc,
                                                precomputed_descriptors=descriptors)
        ml_wf.bands = copy.deepcopy(ab_initio_wf.bands)
        ml_wf.parent = None

        # Attach the model that we want to test to the workflow
        ml_wf.ml_model = model
        ml_wf.ml.train = False
        ml_wf.ml.predict = True

        # Run the workflow
        ml_wf.run(subdirectory=f'ml/snapshot_{i_snapshot + 1}/predict_with_ntrain_{i_train + 1}')

        # Extract and store the error in the orbital energies when using the ML model
        predicted_orbital_energies = ml_wf.calculations[-1].results['eigenvalues'][0]
        actual_orbital_energies = ab_initio_wf.calculations[-1].results['eigenvalues'][0]
        errors[i_train] += [p - a for p, a in zip(predicted_orbital_energies, actual_orbital_energies)]

# Plotting
fig, axarr = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]})

# violin plot
ax = axarr[0]
ax.axhline(0, color='black', linestyle='--', lw=1)
ax.violinplot(errors, showmeans=True)
ax.set_ylabel('errors in orbital \n energies (eV)')

# MAE plot
ax = axarr[1]
ax.plot(range(1, n_train + 1), [np.mean(np.abs(err)) for err in errors], 'o--')
ax.set_yscale('log')
ax.set_ylabel('MAE in orbital \n energies (eV)')
ax.set_xlabel('number of configurations used for training')
ax.set_ylim(1e-2, 1e0)

# Save the figure
plt.tight_layout()
plt.savefig('ntrain_convergence.png', dpi=300)
