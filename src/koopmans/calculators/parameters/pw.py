from pydantic import Field
from pydantic_espresso.models.qe_7_4.pw import (ControlNamelist,
                                                ElectronsNamelist,
                                                SystemNamelist)

from koopmans.pydantic import BaseModel


class PWInputParameters(BaseModel):
    control: ControlNamelist = Field(default_factory=lambda: ControlNamelist())
    system: SystemNamelist = Field(default_factory=lambda: SystemNamelist())
    electrons: ElectronsNamelist = Field(default_factory=lambda: ElectronsNamelist())
