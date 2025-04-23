from typing import Any, Dict

from ase_koopmans.io.espresso import pw_keys

from ._utils import IbravDict, SettingsDict


class PWSettingsDict(IbravDict, SettingsDict):
    """Settings for a Quantum ESPRESSO pw.x calculator."""

    def __init__(self, **kwargs) -> None:
        # Get rid of any nested kwargs
        flattened_kwargs: Dict[str, Any] = {}
        for k, v in kwargs.items():
            if isinstance(v, dict) and k != 'pseudopotentials':
                flattened_kwargs.update(**v)
            else:
                flattened_kwargs[k] = v

        super().__init__(valid=[k for block in pw_keys.values() for k in block],
                         defaults={'calculation': 'scf', 'outdir': './TMP/', 'prefix': 'kc',
                                   'conv_thr': '2.0e-9*nelec', 'verbosity': 'high'},
                         are_paths=['outdir', 'pseudo_dir'],
                         **flattened_kwargs)

    def is_valid(self, name: str) -> bool:
        """Check if a keyword is valid for this calculator.

        Allow for keywords such as ion_radius(i) where i is some arbitrary index
        """
        if 'celldm' in name:
            return super().is_valid(name)
        else:
            return super().is_valid(name.split('(')[0])

    def __setitem__(self, key: str, value: Any):
        if key == 'nspin' and value == 1:
            self.data.pop('tot_magnetization', None)
        return super().__setitem__(key, value)
