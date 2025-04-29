"""Pydantic model for inputs/outputs of Processes."""

import numpy as np
from pydantic import BaseModel


class IOModel(BaseModel):
    """A base class for input/output models of Processes."""

    model_config = {'arbitrary_types_allowed': True, 'extra': 'forbid'}

    def __eq__(self, other):
        """Check if two instances of the class have matching contents (need not be the same instance)."""
        if not isinstance(other, self.__class__):
            return False
        for key in self.dict().keys():
            cond = (getattr(self, key) == getattr(other, key))
            if isinstance(cond, np.ndarray):
                cond = cond.all()
            if not cond:
                return False
        return True

    def todict(self):
        """Convert the class to a dictionary."""
        dct = self.model_dump()
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__module__
        return dct

    @classmethod
    def fromdict(cls, dct):
        """Create an instance of the class from a dictionary."""
        return cls(**dct)
