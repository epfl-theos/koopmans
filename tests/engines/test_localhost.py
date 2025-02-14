import pytest

from koopmans import utils
from koopmans.engines.localhost import LocalhostEngine
from koopmans.outputs import OutputModel
from koopmans.workflows import Workflow
from koopmans.status import Status


class DummyOutput(OutputModel):
    message: 'str'


def test_localhost(silicon, tmp_path):

    dummy_message = 'Dummy workflow completed'

    class dummy_workflow(Workflow):

        output_model = DummyOutput

        def _run(self):
            self.status = Status.COMPLETED
            self.outputs = DummyOutput(message=dummy_message)

    with utils.chdir(tmp_path):
        workflow = dummy_workflow(**silicon)
        workflow.run()
        assert workflow.status == Status.COMPLETED
        assert workflow.outputs.message == dummy_message
