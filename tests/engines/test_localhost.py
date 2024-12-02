import pytest

from koopmans import utils
from koopmans.engines.localhost import LocalhostEngine
from koopmans.outputs import OutputModel
from koopmans.workflows import Workflow


class DummyOutput(OutputModel):
    message: 'str'


def test_localhost(silicon):

    class dummy_workflow(Workflow):

        output_model = DummyOutput

        def _steps_generator(self):
            for i in range(10):
                utils.warn(f'This is step {i} of a dummy workflow')
                yield tuple()

            self.outputs = DummyOutput(message='Dummy workflow completed')

    engine = LocalhostEngine()
    workflow = dummy_workflow(**silicon)
    engine.run(workflow)
