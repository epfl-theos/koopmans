from ase.io.espresso import koopmans_ham as kch_io
from ._utils import SettingsDict, kc_wann_defaults


class KoopmansHamSettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:
        super().__init__(valid=[k for block in kch_io.KEYS.values() for k in block],
                         defaults={'do_bands': True, 'use_ws_distance': True, 'write_hr': True,
                                   'l_alpha_corr': False, 'lrpa': False, **kc_wann_defaults},
                         are_paths=['outdir', 'pseudo_dir'],
                         to_not_parse=['assume_isolated'],
                         **kwargs)
