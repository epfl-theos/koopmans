from pathlib import Path
from typing import Any, Dict

from ase.io.espresso import kcp_keys

from ._utils import IbravDict, SettingsDict


class KoopmansCPSettingsDict(IbravDict, SettingsDict):
    def __init__(self, **kwargs) -> None:
        defaults = {'calculation': 'cp',
                    'outdir': Path('./TMP-CP/'),
                    'iprint': 1,
                    'prefix': 'kc',
                    'verbosity': 'low',
                    'disk_io': 'high',
                    'write_hr': False,
                    'do_wf_cmplx': True,
                    'do_ee': True,
                    'electron_dynamics': 'cg',
                    'nspin': 2,
                    'ortho_para': 1,
                    'passop': 2.0,
                    'ion_dynamics': 'none',
                    'ion_nstepe': 5,
                    'ion_radius(1)': 1.0,
                    'ion_radius(2)': 1.0,
                    'ion_radius(3)': 1.0,
                    'ion_radius(4)': 1.0,
                    'do_innerloop_cg': True,
                    'innerloop_cg_nreset': 20,
                    'innerloop_cg_nsd': 2,
                    'innerloop_init_n': 3,
                    'innerloop_nmax': 100,
                    'hartree_only_sic': False,
                    'conv_thr': '1.0e-9*nelec',
                    'esic_conv_thr': '1.0e-9*nelec',
                    'print_real_space_density': False
                    }

        # Get rid of any nested kwargs
        flattened_kwargs: Dict[str, Any] = {}
        for k, v in kwargs.items():
            if isinstance(v, dict) and k != 'pseudopotentials':
                flattened_kwargs.update(**v)
            else:
                flattened_kwargs[k] = v

        super().__init__(valid=[k for sublist in kcp_keys.values() for k in sublist],
                         defaults=defaults,
                         are_paths=['outdir', 'pseudo_dir'],
                         to_not_parse=['assume_isolated'],
                         **flattened_kwargs)

    @property
    def _other_valid_keywords(self):
        # Make k-point-related keywords invalid
        return ['pseudopotentials']

    def is_valid(self, name):
        # Allow for keywords such as ion_radius(i) where i is some arbitrary index
        if 'celldm' in name:
            return super().is_valid(name)
        else:
            return super().is_valid(name.split('(')[0])
