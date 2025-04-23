"""A simple script that converges the HOMO energy of a water molecule with respect to nr1b, nr2b, and nr3b."""

from ase_koopmans.build import molecule

from koopmans import engines, workflows

# Initialize the engine to run the workflow
engine = engines.LocalhostEngine()

# Use ASE to construct a water molecule
atoms = molecule('H2O', vacuum=5.0)

# Create a subworkflow which calculates (among other things) the PBE HOMO energy of water
subworkflow = workflows.DFTCPWorkflow(atoms=atoms, ecutwfc=30.0, base_functional='pbe',
                                      pseudo_library='PseudoDojo/0.4/PBE/SR/standard/upf',
                                      calculator_parameters={'kcp': {'nr1b': 6, 'nr2b': 6, 'nr3b': 6}},
                                      engine=engine)

# koopmans doesn't implement convergence with respect to nrb, so we need to define a custom
# ConvergenceVariable. To do so, we must first...
# ... define a function that extracts nr1-3b from a workflow


def get_nrb(workflow):
    """Get the values for nr1-3b from a workflow."""
    return [workflows.get_calculator_parameter(workflow, f'nr{i}b') for i in [1, 2, 3]]

# ... define a function that sets nr1-3b


def set_nrb(workflow, value):
    """Set nr1-3b equal to `value` for a workflow."""
    workflows.set_calculator_parameter(workflow, 'nr1b', value[0])
    workflows.set_calculator_parameter(workflow, 'nr2b', value[1])
    workflows.set_calculator_parameter(workflow, 'nr3b', value[2])


nrb_variable = workflows.ConvergenceVariable(name='nrb',
                                             increment=[2, 2, 2],
                                             get_value=get_nrb,
                                             set_value=set_nrb)

# Create the convergence workflow using the convergence factory
workflow = workflows.ConvergenceWorkflowFactory(subworkflow,
                                                observable='total energy',
                                                threshold=1e-3,
                                                variables=[nrb_variable])

# Run the workflow
workflow.run()
