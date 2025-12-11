"""pw calculator module for koopmans."""

import numpy as np
from ase_koopmans import Atoms
from ase_koopmans.calculators.espresso import Espresso
from ase_koopmans.dft.kpoints import BandPath
from ase_koopmans.spectrum.band_structure import BandStructure

from koopmans.cell import cell_follows_qe_conventions, cell_to_parameters
from koopmans.parameters.kcp import KCPInputParameters

from .calculator import CalculatorInputs, CalculatorOutputs, CalculatorProcess


class KCPInputs(CalculatorInputs):
    parameters: KCPInputParameters
    pseudopotentials: dict[str, str]
    pseudo_family: str


class KCPOutputs(CalculatorOutputs):
    total_energy: float | None = None


class KCPCalculator(CalculatorProcess[KCPInputs, KCPOutputs]):
    """Calculator process for performing calculations with kcp.x."""

    ext_in = '.cpi'
    ext_out = '.cpo'
    code = "kcp"
    id = "kcp"
    input_model = KCPInputs
    output_model = KCPOutputs
