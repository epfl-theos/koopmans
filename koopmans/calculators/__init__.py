'Calculators for use with koopmans'
from typing import Type, Union

from ._environ import EnvironCalculator
from ._koopmans_cp import KoopmansCPCalculator, convert_flat_alphas_for_kcp
from ._koopmans_ham import KoopmansHamCalculator
from ._koopmans_screen import KoopmansScreenCalculator
from ._ph import PhCalculator
from ._projwfc import ProjwfcCalculator
from ._pw import PWCalculator
from ._pw2wannier import PW2WannierCalculator
from ._ui import UnfoldAndInterpolateCalculator
from ._utils import (CalculatorCanEnforceSpinSym, CalculatorExt,
                     ReturnsBandStructure)
from ._wann2kc import Wann2KCCalculator
from ._wann2kcp import Wann2KCPCalculator
from ._wannier90 import Wannier90Calculator

Calc = Union[EnvironCalculator,
             KoopmansCPCalculator,
             KoopmansHamCalculator,
             KoopmansScreenCalculator,
             ProjwfcCalculator,
             PWCalculator,
             PW2WannierCalculator,
             Wann2KCPCalculator,
             Wannier90Calculator,
             Wann2KCCalculator,
             UnfoldAndInterpolateCalculator]

CalcType = Union[Type[EnvironCalculator],
                 Type[KoopmansCPCalculator],
                 Type[KoopmansHamCalculator],
                 Type[KoopmansScreenCalculator],
                 Type[ProjwfcCalculator],
                 Type[PWCalculator],
                 Type[PW2WannierCalculator],
                 Type[Wann2KCPCalculator],
                 Type[Wannier90Calculator],
                 Type[Wann2KCCalculator],
                 Type[UnfoldAndInterpolateCalculator]]
