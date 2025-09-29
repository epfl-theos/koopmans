from enum import Enum


class Spin(Enum):
    """Different possibilities for the spin channel of a system."""

    UP = "up"
    DOWN = "down"
    SPINOR = "spinor"
    NONE = None

    def __str__(self):
        return str(self.value)
