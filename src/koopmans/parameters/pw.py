
from pathlib import Path
from typing import Any

from pydantic import Field, field_validator
from pydantic_espresso.models.qe_7_4.pw import \
    ControlNamelist as _ControlNamelist
from pydantic_espresso.models.qe_7_4.pw import (ElectronsNamelist,
                                                SystemNamelist)

from koopmans.pydantic import BaseModel


class ControlNamelist(_ControlNamelist):

    @field_validator('outdir', mode='before')
    @classmethod
    def enforce_outdir_is_tmp(cls, v: Any) -> Path:
        return Path('tmp')


class PWInputParameters(BaseModel):
    control: ControlNamelist = Field(default_factory=lambda: ControlNamelist())
    system: SystemNamelist = Field(default_factory=lambda: SystemNamelist())
    electrons: ElectronsNamelist = Field(default_factory=lambda: ElectronsNamelist())
