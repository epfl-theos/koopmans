'''

Outward-facing module for dealing with settings

Written by Edward Linscott May 2020

'''

from ._koopmans_cp import KoopmansCPSettingsDict
from ._koopmans_ham import KoopmansHamSettingsDict
from ._koopmans_screen import KoopmansScreenSettingsDict
from ._ml import MLSettingsDict
from ._ph import PhSettingsDict
from ._plot import PlotSettingsDict
from ._projwfc import ProjwfcSettingsDict
from ._pw import PWSettingsDict
from ._pw2wannier import PW2WannierSettingsDict
from ._ui import UnfoldAndInterpolateSettingsDict
from ._utils import Setting, SettingsDict, SettingsDictWithChecks
from ._wann2kc import Wann2KCSettingsDict
from ._wann2kcp import Wann2KCPSettingsDict
from ._wannier90 import Wannier90SettingsDict
from ._workflow import WorkflowSettingsDict
