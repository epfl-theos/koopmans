"""pw calculator module for koopmans."""

from koopmans.parameters.wannier90 import Wannier90InputParameters

from .calculator import CalculatorInputs, CalculatorOutputs, CalculatorProcess


class Wannier90Inputs(CalculatorInputs):
    parameters: Wannier90InputParameters


class Wannier90Outputs(CalculatorOutputs):
    pass


class Wannier90Calculator(CalculatorProcess[Wannier90Inputs, Wannier90Outputs]):
    """Calculator process for performing calculations with wannier90.x."""

    ext_in = '.win'
    ext_out = '.wout'
    code = "wannier90"
    id = "wannier90"
    input_model = Wannier90Inputs
    output_model = Wannier90Outputs
