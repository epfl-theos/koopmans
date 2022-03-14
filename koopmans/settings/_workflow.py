import os
from pathlib import Path
from typing import Any
from koopmans import pseudopotentials
from ._utils import Setting, SettingsDictWithChecks


class WorkflowSettingsDict(SettingsDictWithChecks):

    def __init__(self, **kwargs) -> None:
        settings = [
            Setting('task',
                    'Task to perform',
                    str, 'singlepoint', ('singlepoint', 'convergence', 'wannierise', 'environ_dscf', 'ui',
                                         'dft_bands')),
            Setting('functional',
                    'orbital-density-dependent-functional/density-functional to use',
                    str, 'ki', ('ki', 'kipz', 'pkipz', 'dft', 'all')),
            Setting('base_functional',
                    'base functional to use',
                    str, 'pbe', ('lda', 'pbe', 'pbesol')),
            Setting('calculate_alpha',
                    'whether or not to calculate the screening parameters ab-initio',
                    bool, True, (True, False)),
            Setting('pseudo_library',
                    'the pseudopotential library to use (valid options depend on the value of base_functional)',
                    str, None, None),
            Setting('method',
                    'the method to calculate the screening parameters: either with ΔSCF or DFPT',
                    str, 'dscf', ('dscf', 'dfpt')),
            Setting('init_orbitals',
                    'which orbitals to use as an initial guess for the variational orbitals',
                    str, 'pz', ('pz', 'kohn-sham', 'mlwfs', 'projwfs')),
            Setting('init_empty_orbitals',
                    'which orbitals to use as an initial guess for the empty variational orbitals '
                    '(defaults to the same value as "init_orbitals")',
                    str, 'same', ('same', 'pz', 'kohn-sham', 'mlwfs', 'projwfs')),
            Setting('frozen_orbitals',
                    "if True, freeze the variational orbitals for the duration of the calculation once they've been "
                    "initialised",
                    bool, None, (True, False)),
            Setting('periodic',
                    'whether or not the system is periodic',
                    bool, None, (True, False)),
            Setting('calculate_bands',
                    'whether or not to calculate the band structure of the system (if relevant)',
                    bool, True, (True, False)),
            Setting('spin_polarised',
                    'if True, the system will be allowed to break spin symmetry i.e. n^{up}(r) != n^{down}(r)',
                    bool, False, (True, False)),
            Setting('fix_spin_contamination',
                    'if True, steps will be taken to try and avoid spin contamination. This is only sensible when '
                    'performing a non-spin-polarised calculation, and is turned on by default for such calculations',
                    bool, None, (True, False)),
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

    def __setitem__(self, key: str, value: Any):
        # Be forgiving to Americans
        if key == 'task' and value == 'wannierize':
            value = 'wannierise'

        # Make sure that orbital_groups is always stored as a list of lists
        if key == 'orbital_groups' and value is not None:
            if len(value) == 0 or not isinstance(value[0], list):
                value = [value]

        # Make sure that pseudo libraries shortcuts (e.g. "sg15") are converted to the explicit version
        # (e.g. "sg15_v1.2")
        if key == 'pseudo_library':
            full_path = pseudopotentials.pseudos_directory / value
            if full_path.is_symlink():
                value = Path(os.path.realpath(full_path)).name

        return super().__setitem__(key, value)
