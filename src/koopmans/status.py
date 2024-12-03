from enum import Enum


class Status(Enum):
    NOT_STARTED = 'not started'
    RUNNING = 'running'
    COMPLETED = 'completed'
    FAILED = 'failed'
