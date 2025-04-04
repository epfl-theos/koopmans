from typing import Any

from ase_koopmans.io.espresso import kch_keys

from ._utils import SettingsDict, kcw_defaults


class KoopmansHamSettingsDict(SettingsDict):
    """Settings for a Quantum ESPRESSO kcw.x Hamiltonian calculation."""

    def __init__(self, **kwargs) -> None:
        super().__init__(valid=[k for block in kch_keys.values() for k in block],
                         defaults={'calculation': 'ham', 'do_bands': True, 'use_ws_distance': True, 'write_hr': True,
                                   'l_alpha_corr': False, **kcw_defaults},
                         are_paths=['outdir', 'pseudo_dir'],
                         **kwargs)

    @property
    def _other_valid_keywords(self):
        # kpts should correspond to a BandPath object
        # kgrid is used to set mp1-3
        return ['kpts', 'kgrid']

    def __setitem__(self, key: str, value: Any):
        if key == 'kgrid':
            assert isinstance(value, list)
            assert len(value) == 3
            for i, x in enumerate(value):
                self.__setitem__(f'mp{i + 1}', x)
        else:
            return super().__setitem__(key, value)
