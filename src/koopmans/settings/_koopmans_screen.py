from typing import Any

from ase.io.espresso import kcs_keys

from ._utils import SettingsDict, kc_wann_defaults


class KoopmansScreenSettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:
        super().__init__(valid=[k for block in kcs_keys.values() for k in block],
                         defaults={'calculation': 'screen', 'tr2': 1.0e-18, 'nmix': 4, 'niter': 33,
                                   'check_spread': True, **kc_wann_defaults},
                         are_paths=['outdir', 'pseudo_dir'],
                         to_not_parse=['assume_isolated'],
                         ** kwargs)

    @property
    def _other_valid_keywords(self):
        # Accept only kpts -- and use that to set mp1-3
        return ['kpts']

    def __setitem__(self, key: str, value: Any):
        if key == 'kpts':
            assert isinstance(value, list)
            assert len(value) == 3
            for i, x in enumerate(value):
                self.__setitem__(f'mp{i+1}', x)
        else:
            return super().__setitem__(key, value)
