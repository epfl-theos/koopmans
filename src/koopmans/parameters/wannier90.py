from typing import Any, ClassVar

from pydantic import create_model
from wannier90_input.models.latest import Wannier90Input
from wannier90_input.models.parameters import (Coordinate,
                                               FractionalCoordinate,
                                               Projection)

from koopmans.pydantic import BaseModel

fields = {field_name: (field_info.annotation, field_info)
          for field_name, field_info in Wannier90Input.model_fields.items()}

# A modified input class which removes input parameters (such as k-points) that are stored centrally in a workflow rather than as calculator parameters
Wannier90InputParameters = create_model(
    "Wannier90InputParameters",
    __base__=BaseModel,
    **{k: v for k, v in fields.items() if k not in ['unit_cell_cart', 'mp_grid', 'kpoints']}
)

# As above, but also excludes auto-generated keywords
RestrictedWannier90InputParameters = create_model(
    "RestrictedWannier90InputParameters",
    __base__=BaseModel,
    **{k: v for k, v in fields.items() if k not in ['num_wann', 'num_bands', 'exclude_bands', 'unit_cell_cart', 'mp_grid', 'kpoints', 'projections']}
)
