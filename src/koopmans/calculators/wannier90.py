"""Wannier90 calculator module for koopmans."""

from ase.spectrum.band_structure import BandStructure
from pydantic import Field

from koopmans.files import File
from koopmans.kpoints import GammaOnlyKpoints, Kpoints
from koopmans.parameters.wannier90 import Wannier90InputParameters
from koopmans.projections import (ExplicitProjectionsBlock,
                                  ImplicitProjectionsBlock)

from .calculator import CalculatorInputs, CalculatorOutputs, CalculatorProcess


class Wannier90Inputs(CalculatorInputs):
    """Inputs for calculations with wannier90.x."""

    parameters: Wannier90InputParameters
    projections: ExplicitProjectionsBlock | ImplicitProjectionsBlock
    kpoints: Kpoints | GammaOnlyKpoints


class Wannier90Outputs(CalculatorOutputs):
    """Outputs of calculations with wannier90.x."""

    nnkp_file: File | None = None
    job_done: bool
    converged: bool | None = None
    walltime: float
    centers: list[list[float]] = Field(default_factory=list)
    spreads: list[float] = Field(default_factory=list)
    band_structure: BandStructure | None = None
    imre_ratio: list[float] = Field(default_factory=list)


class Wannier90Calculator(CalculatorProcess[Wannier90Inputs, Wannier90Outputs]):
    """Calculator process for performing calculations with wannier90.x."""

    ext_in = '.win'
    ext_out = '.wout'
    code = "wannier90"
    id = "wannier90"
    input_model = Wannier90Inputs
    output_model = Wannier90Outputs
