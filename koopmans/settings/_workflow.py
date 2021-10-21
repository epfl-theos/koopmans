import numpy as np
from ._utils import Setting, SettingsDictWithChecks


class WorkflowSettingsDict(SettingsDictWithChecks):

    def __init__(self, **kwargs) -> None:
        settings = [
            Setting('task',
                    'Task to perform',
                    str, 'singlepoint', ('singlepoint', 'convergence', 'wannierise', 'environ_dscf', 'ui', 'dft_bands')),
            Setting('functional',
                    'orbital-density-dependent-functional/density-functional to use',
                    str, 'ki', ('ki', 'kipz', 'pkipz', 'dft', 'all')),
            Setting('calculate_alpha',
                    'whether or not to calculate the screening parameters ab-initio',
                    bool, True, (True, False)),
            Setting('method',
                    'the method to calculate the screening parameters: either with ΔSCF or DFPT',
                    str, 'dscf', ('dscf', 'dfpt')),
            Setting('init_orbitals',
                    'which orbitals to use as an initial guess for the variational orbitals',
                    str, 'pz', ('pz', 'kohn-sham', 'mlwfs', 'projwfs', 'from old ki')),
            Setting('init_empty_orbitals',
                    'which orbitals to use as an initial guess for the empty variational orbitals '
                    '(defaults to the same value as "init_orbitals")',
                    str, 'same', ('same', 'pz', 'kohn-sham', 'mlwfs', 'projwfs', 'from old ki')),
            Setting('frozen_orbitals',
                    "if True, freeze the variational orbitals for the duration of the calculation once they've been "
                    "initialised",
                    bool, False, (True, False)),
            Setting('periodic',
                    'whether or not the system is periodic',
                    bool, False, (True, False)),
            Setting('npool',
                    'Number of pools for parallelising over kpoints (should be commensurate with the k-point grid)',
                    int, None, None),
            Setting('gb_correction',
                    'if True, apply the Gygi-Baldereschi scheme to deal with the q->0 divergence of the Coulomb '
                    'interation for periodic systems',
                    bool, None, (True, False)),
            Setting('mp_correction',
                    'if True, apply the Makov-Payne correction for charged periodic systems',
                    bool, False, (True, False)),
            Setting('mt_correction',
                    'if True, apply the Martyna-Tuckerman correction for charged aperiodic systems',
                    bool, None, (True, False)),
            Setting('eps_inf',
                    'dielectric constant of the system used by the Gygi-Baldereschi and Makov-Payne corrections',
                    float, None, None),
            Setting('n_max_sc_steps',
                    'maximum number of self-consistency steps for calculating alpha',
                    int, 1, None),
            Setting('alpha_conv_thr',
                    'convergence threshold for |Delta E_i - epsilon_i|; if below this '
                    'threshold, the corresponding alpha value is not updated',
                    (float, str), 1e-3, None),
            Setting('alpha_guess',
                    'starting guess for alpha (overridden if alpha_from_file is true)',
                    (float, list), 0.6, None),
            Setting('alpha_from_file',
                    'if True, uses the file_alpharef.txt from the base directory as a '
                    'starting guess',
                    bool, False, (True, False)),
            Setting('print_qc',
                    'if True, prints out strings for the purposes of quality control',
                    bool, False, (True, False)),
            Setting('from_scratch',
                    'if True, will delete any preexisting workflow and start again; '
                    'if False, will resume a workflow from where it was last up to',
                    bool, False, (True, False)),
            Setting('orbital_groups',
                    'a list of integers the same length as the total number of bands, '
                    'denoting which bands to assign the same screening parameter to',
                    list, None, None),
            Setting('orbital_groups_self_hartree_tol',
                    'when calculating alpha parameters, the code will group orbitals '
                    'together only if their self-Hartree energy is within this '
                    'threshold',
                    float, None, None),
            Setting('enforce_spin_symmetry',
                    'if True, the spin-up and spin-down wavefunctions will be forced '
                    'to be the same',
                    bool, True, (True, False)),
            Setting('check_wannierisation',
                    'if True, checks the Im/Re ratio and generates a plot of the interpolated band structure',
                    bool, False, (True, False)),
            Setting('convergence_observable',
                    'System observable of interest which we converge',
                    str, 'total energy', ('homo energy', 'lumo energy', 'total energy')),
            Setting('convergence_threshold',
                    'Convergence threshold for the system observable of interest',
                    (str, float), None, None),
            Setting('convergence_parameters',
                    'The observable of interest will be converged with respect to this/these '
                    'simulation parameter(s)',
                    (list, str), ['ecutwfc'], None),
            Setting('eps_cavity',
                    'a list of epsilon_infinity values for the cavity in dscf calculations',
                    list, [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20], None)]
        super().__init__(settings=settings, physicals=['alpha_conv_thr', 'convergence_threshold'], **kwargs)

    @property
    def _other_valid_keywords(self):
        return []
