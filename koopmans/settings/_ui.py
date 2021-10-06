"""
Settings module for the UI calculator

Originally written by Riccardo De Gennaro as part of the standalone 'unfolding and interpolate' code
Integrated within koopmans by Edward Linscott Jan 2021

"""

import numpy as np
from typing import List, Any
from pathlib import Path
from ase.dft.kpoints import BandPath
from ._utils import Setting, SettingsDictWithChecks

valid_settings: List[Setting] = [
    Setting('kc_ham_file',
            'the name of the Hamiltonian file to read in',
            Path, None, None),
    Setting('w90_seedname',
            'w90_seedname must be equal to the seedname used in the previous Wannier90 calculation. The code '
            'will look for a file called w90_seedname.wout',
            Path, None, None),
    Setting('sc_dim',
            'units of repetition of the primitive cell withing the supercell along the three lattice directions. '
            'Equivalently this has to match the Monkhorst-Pack mesh of k-points.',
            list, None, None),
    Setting('w90_calc',
            'Specifies the type of PW/Wannier90 calculation preceding the koopmans calculation. If the latter '
            'is done in a supercell at Gamma then w90_calc must be equal to \'sc\', otherwise if it comes from '
            'a calculation with k-points it must be equal to \'pc\'.\n',
            str, 'pc', ('pc', 'sc')),
    Setting('do_map',
            'if True, it realizes the map |m> --> |Rn>, that connects the Wannier functions in the supercell to '
            'those in the primitive cell. This is basically the unfolding procedure. It can be activated only '
            'if w90_calc=\'sc\'',
            bool, False, (True, False)),
    Setting('use_ws_distance',
            'if True, the real Wigner-Seitz distance between the Wannier functions centers is considered as in '
            'the Wannier90 code. In particular, this accounts for the periodic boundary conditions and it is '
            'crucial for a good interpolation when using coarse MP meshes or, equivalently, small supercells',
            bool, True, (True, False)),
    Setting('kpath',
            'path in the Brillouin zone for generating the band structure, specified by a string e.g. "GXG"',
            (str, BandPath), None, None),
    Setting('smooth_int_factor',
            'if this is > 1 (or is a 3-element list with at least one entry > 1), the smooth interpolation '
            'method is used. This consists of removing the DFT part of the Hamiltonian from the full Koopmans '
            'Hamiltonian and adding the DFT Hamiltonian from a calculation with a denser k-points mesh, where '
            'this keyword defines how many times denser to make the mesh. (If this is set to a scalar a, the '
            'new k-grid will be [a*kx_old, a*ky_old, a*kz_old]. If it is a list [a, b, c], the dense k-grid '
            'will be [a*kx_old, b*ky_old, c*kz_old].) This works only for a non self-consistent Koopmans '
            'calculation using Wannier since, to be consistent, all the Hamiltonians must be in the same '
            'gauge, i.e. the Wannier gauge',
            (int, list), 1, None),
    Setting('dft_ham_file', '', Path, None, None),
    Setting('dft_smooth_ham_file', '', Path, None, None),
    Setting('do_dos',
            'if True, the density-of-states is interpolated along the input kpath. The DOS is written to a '
            'file called "dos_interpolated.dat"',
            bool, True, (True, False)),
    Setting('degauss',
            'gaussian broadening (in eV) for the DOS interpolation, as in QE',
            (str, float, int), 0.05, None),
    Setting('nstep',
            'number of steps for the plot of the interpolated DOS',
            int, 1000, None),
    Setting('Emin',
            'minimum energy for the plot of the interpolated DOS',
            (str, float, int), None, None),
    Setting('Emax',
            'maximum energy for the plot of the interpolated DOS',
            (str, float, int), None, None)]


class UnfoldAndInterpolateSettingsDict(SettingsDictWithChecks):

    def __init__(self, **kwargs):
        # First, rename a few settings internally
        # - sc_dim is nothing more than 'kpts', so let's internally rename this
        name_changes = {'sc_dim': 'kpts'}
        local_valid_settings = []
        for setting in valid_settings:
            if setting.name in name_changes.keys():
                setting_dict = setting._asdict()
                setting_dict['name'] = name_changes[setting.name]
                setting = Setting(**setting_dict)

            local_valid_settings.append(setting)

        super().__init__(settings=local_valid_settings,
                         to_not_parse=[],
                         physicals=['degauss', 'Emin', 'Emax'],
                         **kwargs)

    @property
    def _other_valid_keywords(self):
        return ['kpts', 'kpath']

    def _check_before_setitem(self, key, value):
        # Additional sanity checks
        if self.do_map and key == 'w90_calc' and value.lower() != 'sc':
            raise ValueError('do_map = True is incompatible with w90_calc = "pc"')

        return super()._check_before_setitem(key, value)

    def __setitem__(self, key: str, value: Any):
        if key == 'w90_calc':
            value = value.lower()
            if value == 'sc':
                self.w90_input_sc = True
            else:
                self.w90_input_sc = False

        if key == 'smooth_int_factor' and isinstance(value, int):
            value = [value for _ in range(3)]

        return super().__setitem__(key, value)

    @property
    def do_smooth_interpolation(self):
        return any([f > 1 for f in self.smooth_int_factor])
