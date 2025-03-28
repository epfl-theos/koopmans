"""Class that defines the status of a step in a workflow."""

from enum import Enum


class Status(Enum):
    """Status of a step in a workflow."""

    NOT_STARTED = 'not started'
    RUNNING = 'running'
    COMPLETED = 'completed'
    FAILED = 'failed'
