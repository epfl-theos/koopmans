from typing import Any

import numpy as np
from ase.dft.kpoints import BandPath
from ase.io.wannier90 import construct_kpoint_path, proj_string_to_dict

from ._utils import SettingsDict


class Wannier90SettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:

        super().__init__(valid=['num_bands', 'num_wann', 'exclude_bands', 'kpoints',
                                'num_iter', 'conv_window', 'conv_tol', 'num_print_cycles',
                                'dis_froz_max', 'dis_num_iter', 'dis_win_max', 'guiding_centres',
                                'bands_plot', 'mp_grid', 'kpoint_path', 'projections', 'write_hr',
                                'write_u_matrices', 'write_xyz', 'wannier_plot', 'wannier_plot_list',
                                'gamma_only', 'spin', 'use_ws_distance', 'translate_home_cell',
                                'translation_centre_frac'],
                         defaults={'num_iter': 10000, 'conv_tol': 1.e-10, 'conv_window': 5,
                                   'write_hr': True, 'guiding_centres': True, 'gamma_only': False},
                         to_not_parse=['exclude_bands', 'wannier_plot_list'],
                         **kwargs)

    def update(self, *args, **kwargs) -> None:
        kpath = kwargs.pop('kpath', None)
        super().update(*args, **kwargs)
        # Defer setting kpath until after bands_plot has had the chance to be updated
        if kpath is not None and self.get('bands_plot', False):
            self.kpath = kpath

    @property
    def _other_valid_keywords(self):
        return ['kgrid', 'kpath']

    def __setitem__(self, key: str, value: Any):
        if key == 'kgrid':
            if value is None:
                # When kgrid is None, i.e. we are performing a Γ-only calculation,
                # Wannier90 still wants the mp_grid to be defined
                value = [1, 1, 1]

            # Making sure the input is the correct format (a list of length 3)
            assert isinstance(value, list)
            assert len(value) == 3

            # Wannier90 has two keywords associated with the kpoints grid

            # The first, mp_grid, is exactly the same as kgrid itself
            self.mp_grid = value

            # The second, kpoints, is a grid of each individual k-point
            kpts: np.ndarray = np.indices(value, dtype=float).transpose(1, 2, 3, 0).reshape(-1, 3)
            kpts /= value
            kpts[kpts >= 0.5] -= 1
            self.kpoints = kpts
        elif key == 'kpath':
            # Wannier90 calls the kpath "kpoint_path'. Furthermore, in Wannier90 the length of this BandPath is
            # specified by bands_plot_num, so we must adjust the input accordingly
            assert isinstance(value, BandPath)
            self.kpoint_path = construct_kpoint_path(
                path=value.path, cell=value.cell, bands_point_num=self.get('bands_point_num', 100))
        elif key == 'bands_point_num':
            super().__setitem__(key, value)
            # Update the bandpath accordingly
            self.kpoint_path = construct_kpoint_path(
                path=self.kpoint_path.path, cell=self.kpoint_path.cell, bands_point_num=value)
        else:
            if key == 'wannier_plot_list':
                assert isinstance(value, str), 'wannier_plot_list must be a string, e.g. "1,3-5"'
            if key == 'kpoint_path':
                assert self.bands_plot, 'Do not try and set a kpoint_path for a Wannier90 calculation which does ' \
                    'not have bands_plot = True'
            if key == 'projections':
                for i, v in enumerate(value):
                    if isinstance(v, str):
                        v = proj_string_to_dict(v)
                    assert isinstance(v, dict)
                    for k in v.keys():
                        if k not in ['site', 'csite', 'fsite', 'ang_mtm', 'zaxis', 'xaxis', 'radial', 'zona']:
                            raise KeyError(f'Unrecognized key {k} in the w90 projections block')
                    value[i] = v
            return super().__setitem__(key, value)
