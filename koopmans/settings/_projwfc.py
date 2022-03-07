from ._utils import SettingsDict


class ProjwfcSettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:

        super().__init__(valid=['outdir', 'prefix', 'filpdos', 'deltae', 'degauss'],
                         defaults={'outdir': 'TMP', 'prefix': 'kc', 'deltae': 0.01, 'degauss': 0.01},
                         are_paths=['outdir'],
                         **kwargs)

    @property
    def _other_valid_keywords(self):
        # Don't accept either kpoints data or pseudopotentials
        return []
