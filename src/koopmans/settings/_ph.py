from ._utils import SettingsDict


class PhSettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:
        super().__init__(valid=['outdir', 'prefix', 'epsil', 'amass',
                                'fildyn', 'tr2_ph', 'trans'],
                         defaults={'outdir': 'TMP', 'prefix': 'kc', 'epsil': True, 'trans': False
                                   },
                         are_paths=['outdir'],
                         **kwargs)

    @property
    def _other_valid_keywords(self):
        return ['pseudopotentials']

    def is_valid(self, name: str) -> bool:
        # Allow for keywords such as amass(i) where i is some arbitrary index
        return super().is_valid(name.split('(')[0])
