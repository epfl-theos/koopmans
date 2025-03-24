from pathlib import Path

from pydantic import BaseModel

from koopmans.engines import LocalhostEngine
from koopmans.processes import Process


class DummyInput(BaseModel):
    question: str


class DummyOutput(BaseModel):
    answer: int


class DummyProcess(Process):

    input_model = DummyInput
    output_model = DummyOutput

    def _run(self):
        self.outputs = self.output_model(answer=42)

    def dump_inputs(self):
        # Can't pickle DummyInput
        pass

    def dump_outputs(self):
        # Can't pickle DummyOutput
        pass


def test_process():
    process = DummyProcess(question="What is the meaning of life, the universe, and everything?")
    process.directory = Path()
    process.engine = LocalhostEngine()
    process.run()
    assert process.outputs.answer == 42
