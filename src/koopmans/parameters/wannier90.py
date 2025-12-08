from typing import Any, ClassVar

from pydantic import model_validator
from wannier90_input.models.latest import Wannier90Input
from wannier90_input.models.parameters import (Coordinate,
                                               FractionalCoordinate,
                                               Projection)


class Wannier90InputParameters(Wannier90Input):
    """A class that represents input parameters unique to Wannier90

    We obscure Wannier90 input parameters (such as k-points) that are stored centrally"""
    num_wann: ClassVar[int] = -1
    num_bands: ClassVar[int | None] = -1
    exclude_bands: ClassVar[list[int]] = []
    unit_cell_cart: ClassVar[list[Coordinate]] = []
    mp_grid: ClassVar[tuple[int, int, int]] = (1, 1, 1)
    kpoints: ClassVar[list[FractionalCoordinate]] = []
    projections: ClassVar[list[Projection]] = []

    @model_validator(mode="before")
    @classmethod
    def set_default_num_bands(cls, values: dict[str, Any]) -> dict[str, Any]:
        """Skip default num_bands setting."""
        return values

    @model_validator(mode="before")
    @classmethod
    def atoms_frac_xor_cart(cls, values: dict[str, Any]) -> dict[str, Any]:
        """Skip atoms_frac/atoms_cart check."""
        return values
