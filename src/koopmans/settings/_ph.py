from ._utils import SettingsDict


class PhSettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:
        super().__init__(valid=['outdir', 'prefix', 'epsil', 'amass',
                                'fildyn', 'tr2_ph'],
                         defaults={'outdir': 'TMP', 'prefix': 'kc', 'epsil': True,
                                   },
                         are_paths=['outdir'],
                         **kwargs)

    @property
    def _other_valid_keywords(self):
        return ['pseudopotentials']

    def is_valid(self, name: str) -> bool:
        # Allow for keywords such as amass(i) where i is some arbitrary index
        return super().is_valid(name.split('(')[0])
