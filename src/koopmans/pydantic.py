"""Pydantic base model to use throughout `koopmans`."""

import numpy as np
from pydantic import BaseModel as _BaseModel
from pydantic import ConfigDict

default_config = ConfigDict(extra="forbid",
                            arbitrary_types_allowed=True,
                            strict=False,
                            validate_assignment=True,
                            json_encoders={np.ndarray: lambda v: v.tolist()}
                            )


class BaseModel(_BaseModel):
    """Base model with a modified default configuration."""

    model_config = default_config

    def __eq__(self, other):
        """Check if two instances of the class have matching contents (need not be the same instance)."""
        if not isinstance(other, self.__class__):
            return False
        for key in self.model_dump().keys():
            cond = (getattr(self, key) == getattr(other, key))
            if isinstance(cond, np.ndarray):
                cond = cond.all()
            if not cond:
                return False
        return True
