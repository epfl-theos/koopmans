from enum import Enum


class Spin(Enum):
    UP = "up"
    DOWN = "down"
    SPINOR = "spinor"
    NONE = None

    def __str__(self):
        return str(self.value)
