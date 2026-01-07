"""pw2wannier90 calculator module for koopmans."""

from __future__ import annotations

from typing import Annotated, Any

from pint import Quantity
from pydantic import Field, model_validator
from pydantic_pint import PydanticPintQuantity

from koopmans.files import File, FileViaSymlink
from koopmans.parameters.pw2wannier90 import PW2Wannier90InputParameters

from .calculator import CalculatorInputs, CalculatorOutputs, CalculatorProcess


class PW2Wannier90Inputs(CalculatorInputs):
    """Inputs for calculations with pw2wannier90.x."""

    parameters: PW2Wannier90InputParameters
    nnkp_file: FileViaSymlink
    pw_outdir: FileViaSymlink

    @model_validator(mode='before')
    @classmethod
    def match_seedname_and_prefix(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Ensure that the seedname and prefix match the corresponding filenames."""
        data['parameters'].seedname = data['nnkp_file'].stem
        return data

    @model_validator(mode='before')
    @classmethod
    def set_outdir(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Set the outdir parameter based on the pw_outdir file."""
        data['parameters'].outdir = data['pw_outdir'].path
        return data


class PW2Wannier90Outputs(CalculatorOutputs):
    """Outputs of calculations with pw2wannier90.x."""

    eig_file: File
    mmn_file: File | None = None
    amn_file: File | None = None
    unk_files: list[File] = Field(default_factory=list)


class PW2Wannier90Calculator(CalculatorProcess[PW2Wannier90Inputs, PW2Wannier90Outputs]):
    """Calculator process for performing calculations with wannier90.x."""

    ext_in = '.p2wi'
    ext_out = '.p2wo'
    code = "pw2wannier90"
    id = "pw2wannier90"
    input_model = PW2Wannier90Inputs
    output_model = PW2Wannier90Outputs
