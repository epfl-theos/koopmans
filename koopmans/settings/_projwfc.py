from ._utils import SettingsDict


class ProjwfcSettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:
        # Marija <- replace all of these lists with versions for projwfc
        super().__init__(valid=['outdir', 'prefix', 'filpdos', 'deltae'], 

                defaults={'outdir': 'TMP', 'prefix': 'kc', 'deltae':0.001},
                         are_paths=['outdir'],
                         **kwargs)

    @property
    def _other_valid_keywords(self):
        # Don't accept either kpoints data or pseudpotentials
        return []
