import pytest

from ase.io.espresso import cell_to_ibrav
from koopmans import workflows

# A dummy implementation of the Workflow ABC


class WorkflowExample(workflows.Workflow):

    def _run(self) -> None:
        raise NotImplementedError()


def test_creation_without_cell(silicon):
    atoms = silicon['atoms']

    raise ValueError()
    wf = WorkflowExample()
    pass
