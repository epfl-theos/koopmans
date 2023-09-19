'''
A simple script that performs a convergence test on the number of snapshots for training
'''

from koopmans import io
from koopmans.workflows import (ConvergenceVariable,
                                ConvergenceWorkflowFactory, Workflow)

# Create a subworkflow
subworkflow = io.read('h2o_trajectory_ml.json')

# Define a function that extracts the number of snapshots from a workflow


def get_nsnap(workflow: Workflow):
    return workflow.ml.number_of_training_snapshots

# Define a function that sets the number of snapshots for a workflow


def set_nsnap(workflow: Workflow, value: int):
    workflow.ml.number_of_training_snapshots = value


nsnap_variable = ConvergenceVariable(name='number_of_training_snapshots',
                                     increment=1,
                                     get_value=get_nsnap,
                                     set_value=set_nsnap)

# Define an observable


def get_alphas(workflow: Workflow) -> float:
    import ipdb
    ipdb.set_trace()


# Create the convergence workflow using the convergence factory. Because koopmans knows how to
# converge with respect to ecutwfc, we don't need to implement a custom ConvergenceVariable for it
# and instead can just tell it the variable name
workflow = ConvergenceWorkflowFactory(subworkflow,
                                      observable=get_alphas,
                                      threshold=1e-3,
                                      variables=[nsnap_variable])

# Run the workflow
workflow.run()
