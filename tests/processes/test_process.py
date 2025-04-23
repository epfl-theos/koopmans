"""Testing the `Process` class."""

from pathlib import Path

from pydantic import BaseModel

from koopmans.engines import LocalhostEngine
from koopmans.processes import Process


class DummyInput(BaseModel):
    """A dummy input model."""

    question: str


class DummyOutput(BaseModel):
    """A dummy output model."""

    answer: int


class DummyProcess(Process):
    """A dummy process."""

    input_model = DummyInput
    output_model = DummyOutput

    def _run(self):
        """Run the process."""
        self.outputs = self.output_model(answer=42)

    def dump_inputs(self):
        """Dump the inputs."""
        # Can't pickle DummyInput
        pass

    def dump_outputs(self):
        """Dump the outputs."""
        # Can't pickle DummyOutput
        pass


def test_process():
    """Test the `Process` class."""
    process = DummyProcess(question="What is the meaning of life, the universe, and everything?")
    process.directory = Path()
    process.engine = LocalhostEngine()
    process.run()
    assert process.outputs.answer == 42
