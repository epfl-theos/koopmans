"""
Settings module for the UI calculator

Originally written by Riccardo De Gennaro as part of the standalone 'unfolding and interpolate' code
Integrated within python_KI by Edward Linscott Jan 2021

"""

import numpy as np
from koopmans import utils
from ._utils import generate_path, MP_mesh

valid_settings = [
    utils.Setting('kc_ham_file',
                  'the name of the Hamiltonian file to read in',
                  str, None, None),
    utils.Setting('w90_seedname',
                  'w90_seedname must be equal to the seedname used in the previous Wannier90 calculation. The code '
                  'will look for a file called w90_seedname.wout',
                  str, None, None),
    utils.Setting('alat_sc',
                  'the lattice parameter (in Bohr) of the supercell, as celldm(1) in QE. NB: it is important to put '
                  'the supercell (and not the primitive cell) lattice parameter, otherwise the result will be wrong',
                  (float, int), None, None),
    utils.Setting('sc_dim',
                  'units of repetition of the primitive cell withing the supercell along the three lattice directions. '
                  'Equivalently this has to match the Monkhorst-Pack mesh of k-points.',
                  list, [1, 1, 1], None),
    utils.Setting('w90_calc',
                  'Specifies the type of PW/Wannier90 calculation preceding the koopmans calculation. If the latter '
                  'is done in a supercell at Gamma then w90_calc must be equal to \'sc\', otherwise if it comes from '
                  'a calculation with k-points it must be equal to \'pc\'.\n',
                  str, 'pc', ('pc', 'sc')),
    utils.Setting('do_map',
                  'if True, it realizes the map |m> --> |Rn>, that connects the Wannier functions in the supercell to '
                  'those in the primitive cell. This is basically the unfolding procedure. It can be activated only '
                  'if w90_calc=\'sc\'',
                  bool, False, (True, False)),
    utils.Setting('use_ws_distance',
                  'if True, the real Wigner-Seitz distance between the Wannier functions centers is considered as in '
                  'the Wannier90 code. In particular, this accounts for the periodic boundary conditions and it is '
                  'crucial for a good interpolation when using coarse MP meshes or, equivalently, small supercells',
                  bool, True, (True, False)),
    utils.Setting('k_path',
                  'path in the Brillouin zone for the band structure. The logic is the crystal_b units in QE: the '
                  'path must be defined by providing the initial and final point crystal coordinates of each line, '
                  'followed by the number of points along the line.\n\nExample: path Gamma-X-Gamma with 10 points '
                  'along each line\n\n\tk_path :\t[ [ 0.000, 0.000, 0.000, 10 ],\n\t        \t  [ 0.500, 0.000, '
                  '0.000, 10 ],\n\t        \t  [ 0.000, 0.000, 0.000,  1 ] ]\n\nThe band structure is then written '
                  'into a file called \'bands_interpolated.dat\'.', list, None, None),
    utils.Setting('smooth_int_factor',
                  'if this is > 1 (or is a 3-element list with at least one entry > 1), the smooth interpolation '
                  'method is used. This consists of removing the DFT part of the Hamiltonian from the full Koopmans '
                  'Hamiltonian and adding the DFT Hamiltonian from a calculation with a denser k-points mesh, where '
                  'this keyword defines how many times denser to make the mesh. (If this is set to a scalar a, the '
                  'new k-grid will be [a*kx_old, a*ky_old, a*kz_old]. If it is a list [a, b, c], the dense k-grid '
                  'will be [a*kx_old, b*ky_old, c*kz_old].) This works only for a non self-consistent Koopmans '
                  'calculation using Wannier since, to be consistent, all the Hamiltonians must be in the same '
                  'gauge, i.e. the Wannier gauge',
                  (int, list), 1, None),
    utils.Setting('dft_ham_file', '', str, None, None),
    utils.Setting('dft_smooth_ham_file', '', str, None, None),
    utils.Setting('do_dos',
                  'if True, the density-of-states is interpolated along the input k_path. The DOS is written to a '
                  'file called "dos_interpolated.dat"',
                  bool, True, (True, False)),
    utils.Setting('degauss',
                  'gaussian broadening (in eV) for the DOS interpolation, as in QE',
                  (str, float, int), 0.05, None),
    utils.Setting('nstep',
                  'number of steps for the plot of the interpolated DOS',
                  int, 1000, None),
    utils.Setting('Emin',
                  'minimum energy for the plot of the interpolated DOS',
                  (str, float, int), None, None),
    utils.Setting('Emax',
                  'maximum energy for the plot of the interpolated DOS',
                  (str, float, int), None, None)]


def load_defaults(self):
    self.valid_settings = valid_settings

    # Load UI settings
    self.mandatory_settings = ['w90_seedname', 'kc_ham_file', 'alat_sc', 'sc_dim']
    physicals = ['alat_sc', 'degauss', 'Emin', 'Emax']
    checked_settings = utils.check_settings(self.calc.parameters, self.valid_settings,
                                            physicals=physicals, do_not_lower=self.settings_to_not_parse)
    for key, val in checked_settings.items():
        setattr(self, key, val)

    if self.alat_sc is not None:
        self.alat_sc *= utils.units.Bohr  # conversion to angstroms

    self.w90_calc = self.w90_calc.lower()

    if self.w90_calc == 'sc':
        self.w90_input_sc = True
    else:
        self.w90_input_sc = False
        if self.do_map:
            raise ValueError('do_map = True is incompatible with w90_calc = "pc"')

    if isinstance(self.smooth_int_factor, int):
        self.smooth_int_factor = [self.smooth_int_factor for _ in range(3)]

    if self.k_path is None:
        self.kvec = MP_mesh(self.sc_dim[0], self.sc_dim[1], self.sc_dim[2])
        utils.warn('"k_path" missing in input, the energies will be calculated on a commensurate Monkhorst-Pack mesh')
    else:
        self.kvec = generate_path(self.k_path)

    if self.do_smooth_interpolation:
        assert 'dft_ham_file' is not None, 'Missing file_hr_coarse for smooth interpolation'
        assert 'dft_smooth_ham_file' is not None, 'Missing dft_smooth_ham_file for smooth interpolation'
