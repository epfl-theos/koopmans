"""Pydantic base model to use throughout `koopmans`."""

from pydantic import BaseModel as _BaseModel
from pydantic import ConfigDict


class BaseModel(_BaseModel):
    """Base model with a modified default configuration."""

    model_config = ConfigDict(extra="forbid",
                              arbitrary_types_allowed=True,
                              validate_assignment=True)
