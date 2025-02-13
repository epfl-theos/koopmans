import pytest

from koopmans.engines.localhost import LocalhostEngine
from koopmans.workflows import DFTBandsWorkflow


def test_local_engine(silicon):
    engine = LocalhostEngine()
    workflow = DFTBandsWorkflow(**silicon)
    engine.run(workflow)
