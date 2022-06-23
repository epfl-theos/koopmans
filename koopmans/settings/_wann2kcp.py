from ._utils import SettingsDict


class Wann2KCPSettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:
        super().__init__(valid=['outdir', 'prefix', 'seedname', 'wan_mode', 'spin_component',
                                'gamma_trick', 'print_rho', 'wannier_plot', 'wannier_plot_list'],
                         defaults={'outdir': 'TMP', 'prefix': 'kc', 'seedname': 'wann',
                                   'wan_mode': 'wannier2kcp'},
                         are_paths=['outdir'],
                         to_not_parse=['wan_mode', 'wannier_plot_list'],
                         **kwargs)

    @property
    def _other_valid_keywords(self):
        # Don't accept either kpoints data or pseudpotentials
        return []
