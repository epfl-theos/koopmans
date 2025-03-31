"""Testing `koopmans.engines.localhost`."""

import pytest  # noqa: F401

from koopmans.engines.localhost import LocalhostEngine
from koopmans.process_io import IOModel
from koopmans.status import Status
from koopmans.workflows import Workflow


class DummyOutput(IOModel):
    """Dummy output model for the dummy workflow."""

    message: 'str'


def test_localhost_running_dummy_workflow(silicon):
    """Test the LocalhostEngine by running a dummy workflow."""
    dummy_message = 'Dummy workflow completed'

    class dummy_workflow(Workflow):

        output_model = DummyOutput

        def _run(self):
            self.status = Status.COMPLETED
            self.outputs = DummyOutput(message=dummy_message)

    workflow = dummy_workflow(**silicon)
    workflow.run()
    assert workflow.status == Status.COMPLETED
    assert workflow.outputs.message == dummy_message


def test_localhost_installing_custom_pseudopotentials(datadir):
    """Test installing a custom pseudopotential using a LocalhostEngine."""
    # Install a pseudopotential, fetch it, and check that it is the same file

    engine = LocalhostEngine()

    # Install the pseudopotentials contained in datadir/pseudos
    pseudos = []
    for p in (datadir / 'pseudos').glob('*.upf'):
        pseudos.append(p)
        engine.install_pseudopotential(p, library='testing')

    # Fetch the pseudopotentials and check that they point to the original files
    for p in pseudos:
        upfdict = engine.get_pseudopotential('testing', filename=p.name)
        assert upfdict.filename == p

    # Uninstall the testing library
    engine.uninstall_pseudopotential_library('testing')
