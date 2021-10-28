from typing import Any, Dict
from ase.io.espresso import pw as pw_io
from ase.dft.kpoints import BandPath
from ._utils import SettingsDict


class PWSettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:
        # Get rid of any nested kwargs
        flattened_kwargs: Dict[str, Any] = {}
        for k, v in kwargs.items():
            if isinstance(v, dict) and k != 'pseudopotentials':
                flattened_kwargs.update(**v)
            else:
                flattened_kwargs[k] = v

        super().__init__(valid=[k for block in pw_io.KEYS.values() for k in block],
                         defaults={'calculation': 'scf', 'outdir': './TMP/', 'prefix': 'kc',
                                   'conv_thr': '2.0e-9*nelec'},
                         are_paths=['outdir', 'pseudo_dir'],
                         to_not_parse=['assume_isolated'],
                         **flattened_kwargs)

    def is_valid(self, name: str) -> bool:
        # Allow for keywords such as ion_radius(i) where i is some arbitrary index
        if 'celldm' in name:
            return super().is_valid(name)
        else:
            return super().is_valid(name.split('(')[0])

    def __setitem__(self, key: str, value: Any):
        if key == 'nspin' and value == 1:
            self.data.pop('tot_magnetization', None)
        if key == 'kpts' and isinstance(value, BandPath):
            # Make sure the kpoints are always printed out in full
            if len(value.kpts) >= 100:
                self.verbosity = 'high'
        return super().__setitem__(key, value)
