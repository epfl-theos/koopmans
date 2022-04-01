## MARIJA

from ._utils import SettingsDict


class PhSettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:
        super().__init__(valid=['outdir', 'prefix', 'epsil', 'amass',
                                'fildyn', 'tr2_ph',
                         defaults={'outdir': 'TMP', 'prefix': 'kc', 'epsil': '.true.',
                                   'wan_mode': 'standalone'},
                         are_paths=['outdir'],
                         to_not_parse=['wan_mode'],
                         **kwargs)

    @property
    def _other_valid_keywords(self):
        # Don't accept either kpoints data or pseudpotentials
        return []
