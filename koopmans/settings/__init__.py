'''

Outward-facing module for dealing with settings

Written by Edward Linscott May 2020

'''

from ._koopmans_cp import KoopmansCPSettingsDict
from ._koopmans_ham import KoopmansHamSettingsDict
from ._koopmans_screen import KoopmansScreenSettingsDict
from ._pw import PWSettingsDict
from ._pw2wannier import PW2WannierSettingsDict
from ._ui import UnfoldAndInterpolateSettingsDict
from ._utils import Setting, SettingsDict, SettingsDictWithChecks
from ._wann2kc import Wann2KCSettingsDict
from ._wannier90 import Wannier90SettingsDict
from ._workflow import WorkflowSettingsDict

# Dictionary to be used as the default value for 'master_calc_params' when initialising a workflow
default_master_calc_params = {'kcp': KoopmansCPSettingsDict(),
                              'kc_ham': KoopmansHamSettingsDict(),
                              'kc_screen': KoopmansScreenSettingsDict(),
                              'pw': PWSettingsDict(),
                              'pw2wannier': PW2WannierSettingsDict(),
                              'ui': UnfoldAndInterpolateSettingsDict(),
                              'ui_occ': UnfoldAndInterpolateSettingsDict(),
                              'ui_emp': UnfoldAndInterpolateSettingsDict(),
                              'wann2kc': Wann2KCSettingsDict(),
                              'w90_occ': Wannier90SettingsDict(),
                              'w90_emp': Wannier90SettingsDict()}
