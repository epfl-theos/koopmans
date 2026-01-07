"""projwfc calculator module for koopmans."""

from __future__ import annotations

from typing import Annotated, Any

from ase.spectrum.doscollection import GridDOSCollection
from ase.spectrum.dosdata import GridDOSData
from pint import Quantity
from pydantic import BeforeValidator, Field, model_validator
from pydantic_pint import PydanticPintQuantity

from koopmans.files import File, FileViaSymlink
from koopmans.parameters.projwfc import ProjwfcInputParameters
from koopmans.serialization import deserialize_ase_griddoscollection

from .calculator import CalculatorInputs, CalculatorOutputs, CalculatorProcess


class ProjwfcInputs(CalculatorInputs):
    """Inputs for calculations with projwfc.x."""

    parameters: ProjwfcInputParameters
    pw_outdir: FileViaSymlink
    spin_polarized: bool
    expected_subshells: dict[str, list[str]] = Field(
        description="a list of orbitals (e.g. 1s, 2s, 2p, ...) expected for each element in the system.")

    @model_validator(mode='before')
    @classmethod
    def match_outdir_with_pw_outdir(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Ensure that pw_outdir matches outdir in parameters."""
        if 'pw_outdir' in data and 'parameters' in data:
            pw_outdir = data['pw_outdir']
            data['parameters'].outdir = pw_outdir.path
        return data


class ProjwfcOutputs(CalculatorOutputs):
    """Outputs of calculations with projwfc.x."""

    dos: Annotated[GridDOSCollection, BeforeValidator(deserialize_ase_griddoscollection)]


class ProjwfcCalculator(CalculatorProcess[ProjwfcInputs, ProjwfcOutputs]):
    """Calculator process for performing calculations with projwfc.x."""

    ext_in = '.pri'
    ext_out = '.pro'
    code = "projwfc"
    id = "projwfc"
    input_model = ProjwfcInputs
    output_model = ProjwfcOutputs
