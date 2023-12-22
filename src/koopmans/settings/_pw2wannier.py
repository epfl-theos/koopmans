from ._utils import SettingsDict


class PW2WannierSettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:
        super().__init__(valid=['outdir', 'prefix', 'seedname', 'write_mmn',
                                'write_amn', 'write_uHu', 'uHu_formatted',
                                'write_unk', 'reduce_unk', 'wan_mode', 'spin_component',
                                'atom_proj', 'atom_proj_ext', 'atom_proj_dir'],
                         defaults={'outdir': 'TMP', 'prefix': 'kc', 'seedname': 'wann',
                                   'wan_mode': 'standalone'},
                         are_paths=['outdir', 'atom_proj_dir'],
                         **kwargs)

    @property
    def _other_valid_keywords(self):
        # Don't accept either kpoints data or pseudpotentials
        return []
