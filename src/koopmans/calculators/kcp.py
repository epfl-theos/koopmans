"""kcp calculator module for koopmans."""

from koopmans.parameters.kcp import KCPInputParameters

from .calculator import CalculatorInputs, CalculatorOutputs, CalculatorProcess


class KCPInputs(CalculatorInputs):
    """Inputs for calculations with kcp.x."""

    parameters: KCPInputParameters
    pseudopotentials: dict[str, str]
    pseudo_family: str


class KCPOutputs(CalculatorOutputs):
    """Outputs of calculations with kcp.x."""

    total_energy: float | None = None


class KCPCalculator(CalculatorProcess[KCPInputs, KCPOutputs]):
    """Calculator process for performing calculations with kcp.x."""

    ext_in = '.cpi'
    ext_out = '.cpo'
    code = "kcp"
    id = "kcp"
    input_model = KCPInputs
    output_model = KCPOutputs
