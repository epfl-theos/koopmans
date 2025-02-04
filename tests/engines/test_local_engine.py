import pytest

from koopmans.engines.localhost import LocalEngine
from koopmans.workflows import DFTBandsWorkflow


def test_local_engine(silicon):
    engine = LocalEngine()
    workflow = DFTBandsWorkflow(**silicon)
    engine.run(workflow)
