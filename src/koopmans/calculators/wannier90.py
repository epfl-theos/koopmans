"""Wannier90 calculator module for koopmans."""

from typing import Annotated

from ase.spectrum.band_structure import BandStructure
from pydantic import BeforeValidator, Field

from koopmans.files import File, FileViaSymlink
from koopmans.kpoints import Kpoints
from koopmans.parameters.wannier90 import Wannier90InputParameters
from koopmans.serialization import deserialize_ase_bandstructure

from .calculator import CalculatorInputs, CalculatorOutputs, CalculatorProcess


class Wannier90PreprocInputs(CalculatorInputs):
    """Inputs for preprocessing calculations with wannier90.x."""

    parameters: Wannier90InputParameters
    kpoints: Kpoints


class Wannier90PreprocOutputs(CalculatorOutputs):
    """Outputs of preprocessing calculations with wannier90.x."""

    nnkp_file: File
    unk_files: list[File] = Field(default_factory=list)


class Wannier90PreprocCalculator(CalculatorProcess[Wannier90PreprocInputs, Wannier90PreprocOutputs]):
    """Calculator process for performing calculations with wannier90.x."""

    ext_in = '.win'
    ext_out = '.wout'
    code = "wannier90"
    id = "wannier90"
    input_model = Wannier90PreprocInputs
    output_model = Wannier90PreprocOutputs


class Wannier90Inputs(CalculatorInputs):
    """Inputs for calculations with wannier90.x."""

    parameters: Wannier90InputParameters
    kpoints: Kpoints
    amn_file: FileViaSymlink
    mmn_file: FileViaSymlink
    eig_file: FileViaSymlink
    unk_files: list[FileViaSymlink] = Field(default_factory=list)


class Wannier90Outputs(CalculatorOutputs):
    """Outputs of calculations with wannier90.x."""

    centers: list[list[float]]
    spreads: list[float]
    band_structure: Annotated[BandStructure, BeforeValidator(deserialize_ase_bandstructure)] | None = None
    imre_ratio: list[float] = Field(default_factory=list)
    hr_file: File | None = None
    u_file: File | None = None
    centers_file: File | None = None
    u_dis_file: File | None = None
    chk_file: File


class Wannier90Calculator(CalculatorProcess[Wannier90Inputs, Wannier90Outputs]):
    """Calculator process for performing calculations with wannier90.x."""

    ext_in = '.win'
    ext_out = '.wout'
    code = "wannier90"
    id = "wannier90"
    input_model = Wannier90Inputs
    output_model = Wannier90Outputs
