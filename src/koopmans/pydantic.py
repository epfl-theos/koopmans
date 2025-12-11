"""Pydantic base model to use throughout `koopmans`."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
from pydantic import BaseModel as _BaseModel
from pydantic import ConfigDict

if TYPE_CHECKING:
    from koopmans.files import LinkedFile

default_config = ConfigDict(extra="forbid",
                            arbitrary_types_allowed=True,
                            strict=False,
                            validate_assignment=True,
                            json_encoders={np.ndarray: lambda v: v.tolist()}
                            )


class BaseModel(_BaseModel):
    """Base model with a modified default configuration."""

    model_config = default_config

    def __eq__(self, other: Any):
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

    def fields_that_are_files(self) -> list[LinkedFile]:
        """Return the fields of this model that are Files."""
        from koopmans.files import LinkedFile
        file_inputs: list[LinkedFile] = []
        for field in self.model_fields:
            value = getattr(self, field)
            if isinstance(value, LinkedFile):
                file_inputs.append(value)
        return file_inputs
