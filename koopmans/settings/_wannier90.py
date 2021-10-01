from ._utils import SettingsDict


class Wannier90SettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:

        super().__init__(valid=['num_bands', 'num_wann', 'exclude_bands',
                                'num_iter', 'conv_window', 'conv_tol', 'num_print_cycles',
                                'dis_froz_max', 'dis_num_iter', 'dis_win_max', 'guiding_centres',
                                'bands_plot', 'mp_grid', 'kpoint_path', 'projections', 'write_hr',
                                'write_u_matrices', 'write_xyz', 'wannier_plot', 'gamma_only'],
                         defaults={'num_iter': 10000, 'conv_tol': 1.e-10, 'conv_window': 5,
                                   'write_hr': True, 'guiding_centres': True, 'gamma_only': False},
                         to_not_parse=['exclude_bands'],
                         **kwargs)

    @property
    def _other_valid_keywords(self):
        return ['kpts', 'koffset', 'kpath']
