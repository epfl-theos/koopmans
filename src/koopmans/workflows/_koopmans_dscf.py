"""

Workflow module for koopmans, containing the workflow for performing KI and KIPZ calculations

Written by Edward Linscott Jan 2020
Split off from workflow.py Oct 2020

"""

import shutil
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
from ase.dft import DOS

from koopmans import calculators, utils
from koopmans.bands import Band, Bands
from koopmans.settings import KoopmansCPSettingsDict

from ._workflow import Workflow


class KoopmansDSCFWorkflow(Workflow):

    def __init__(self, *args, redo_smooth_dft: Optional[bool] = None, restart_from_old_ki: bool = False, **kwargs) -> None:
        super().__init__(*args, **kwargs)

        # The following two additional keywords allow for some tweaking of the workflow when running a singlepoint
        # workflow with functional == 'all'

        # By default, we don't override self.parameters.from_scratch when we arrive at the higher-res DFT calculations.
        # We change this flag to False for workflows where this is unnecessary (such as pKIPZ) to skip this step.
        self._redo_smooth_dft = redo_smooth_dft

        # For the KIPZ calculation we restart from the old KI calculation
        self._restart_from_old_ki = restart_from_old_ki

        # If periodic, convert the kcp calculation into a Γ-only supercell calculation
        kcp_params = self.calculator_parameters['kcp']
        if all(self.atoms.pbc):
            spins: List[Optional[str]]
            if self.parameters.spin_polarized:
                spins = ['up', 'down']
                nelecs = [kcp_params.nelup, kcp_params.neldw]
            else:
                spins = [None]
                nelecs = [kcp_params.nelec // 2]

            for spin, nelec in zip(spins, nelecs):
                # Check that we have wannierized every filled orbital
                if self.projections:
                    nbands_occ = self.projections.num_wann(occ=True, spin=spin)

                    if nbands_occ != nelec:
                        raise ValueError('You have configured this calculation to only wannierize a subset of the '
                                         'occupied bands:\n'
                                         f' number of occupied bands = {nelec}\n'
                                         f' number of occupied Wannier functions = {nbands_occ}\n'
                                         'This is incompatible with the subsequent Koopmans '
                                         'calculation.\nPlease modify the wannier90 settings in order to wannierize '
                                         'all of the occupied bands. (You may want to consider taking advantage of the '
                                         '"projections_blocks" functionality if your system has a lot of electrons.)')

                    nbands_emp = self.projections.num_wann(occ=False, spin=spin)
                else:
                    nbands_occ = nelec
                    nbands_emp = self.calculator_parameters['pw'].nbnd - nbands_occ

                # Check the number of empty states has been correctly configured
                spin_info = f'spin {spin} ' if self.parameters.spin_polarized else ''
                if kcp_params.nbnd is None:
                    if nbands_emp != 0:
                        kcp_params.nbnd = nbands_occ + nbands_emp
                elif nbands_occ > kcp_params.nbnd:
                    raise ValueError(f'The value you have provided for nbnd is less than the number of {spin_info}'
                                     f'electrons. Please increase nbnd to at least {nbands_occ}')
                elif kcp_params.nbnd != nbands_occ + nbands_emp:
                    raise ValueError(f'The number of {spin_info}empty states are inconsistent:\n'
                                     f' number of empty bands = {kcp_params.nbnd - nbands_occ}\n'
                                     f' number of empty Wannier functions = {nbands_emp}\n'
                                     'If you have provided "nbnd" explicitly to the kcp calculator, check that it '
                                     'matches with the number of empty projections/bands in your system.')

            # Populating self.parameters.orbital_groups if needed
            # N.B. self.bands.groups is guaranteed to be 2 x num_wann, but self.parameters.orbital_groups
            # is either 1- or 2- long, depending on if we are spin-polarized or not
            if self.parameters.orbital_groups is None:
                orbital_groups: List[List[int]] = []
                i_start = 0
                for nelec in nelecs:
                    i_end = i_start + kcp_params.get('nbnd', nelec) - 1
                    orbital_groups.append(list(range(i_start, i_end + 1)))
                    i_start = i_end + 1
                self.parameters.orbital_groups = orbital_groups

            if not self.kpoints.gamma_only:
                # Update the KCP settings to correspond to a supercell (leaving self.atoms unchanged for the moment)
                self.convert_kcp_to_supercell()

                # Expanding self.parameters.orbital_groups to account for the supercell, grouping equivalent wannier
                # functions together
                for i_spin, nelec in enumerate(nelecs):
                    assert self.kpoints.grid is not None
                    self.parameters.orbital_groups[i_spin] = [i for _ in range(np.prod(self.kpoints.grid))
                                                              for i in self.parameters.orbital_groups[i_spin][:nelec]] \
                        + [i for _ in range(np.prod(self.kpoints.grid)) for i in
                           self.parameters.orbital_groups[i_spin][nelec:]]

        # Check the shape of self.parameters.orbital_groups is as expected
        if self.parameters.spin_polarized:
            target_length = 2
        else:
            target_length = 1
        if self.parameters.orbital_groups is not None:
            assert len(self.parameters.orbital_groups) == target_length

        # Constructing the arrays required to initialize a Bands object
        if self.parameters.spin_polarized:
            if 'nbnd' in kcp_params:
                n_emp_up = kcp_params.nbnd - kcp_params.nelup
                n_emp_dw = kcp_params.nbnd - kcp_params.neldw
            else:
                n_emp_up = 0
                n_emp_dw = 0
            filling = [[True for _ in range(kcp_params.nelup)] + [False for _ in range(n_emp_up)],
                       [True for _ in range(kcp_params.neldw)] + [False for _ in range(n_emp_dw)]]
            groups = self.parameters.orbital_groups
        else:
            if 'nbnd' in kcp_params:
                n_emp = kcp_params.nbnd - kcp_params.nelec // 2
            else:
                n_emp = 0
            filling = [[True for _ in range(kcp_params.nelec // 2)]
                       + [False for _ in range(n_emp)] for _ in range(2)]
            # self.parameters.orbital_groups does not have a spin index
            if self.parameters.orbital_groups is None:
                groups = None
            else:
                groups = [self.parameters.orbital_groups[0] for _ in range(2)]

        # Checking groups and filling are the same dimensions
        if groups is not None:
            for g, f in zip(groups, filling):
                assert len(g) == len(f), 'orbital_groups is the wrong dimension; its length should match the number ' \
                    'of bands'

        # Initialize the bands object
        self.bands = Bands(n_bands=[len(f) for f in filling], n_spin=2, spin_polarized=self.parameters.spin_polarized,
                           filling=filling, groups=groups,
                           self_hartree_tol=self.parameters.orbital_groups_self_hartree_tol)

        if self.parameters.alpha_from_file:
            # Reading alpha values from file
            self.bands.alphas = self.read_alphas_from_file()
        else:
            # Initializing alpha with a guess
            self.bands.alphas = self.parameters.alpha_guess

        # Raise errors if any UI keywords are provided but will be overwritten by the workflow
        for ui_keyword in ['kc_ham_file', 'w90_seedname', 'dft_ham_file', 'dft_smooth_ham_file']:
            for ui_kind in ['occ', 'emp']:
                value = getattr(self.calculator_parameters[f'ui_{ui_kind}'], ui_keyword)
                [default_value] = [s.default for s in self.calculator_parameters['ui'].settings if s.name == ui_keyword]
                if value != default_value:
                    raise ValueError(f'UI keyword {ui_keyword} has been set in the input file, but this will be '
                                     'automatically set by the Koopmans workflow. Remove this keyword from the input '
                                     'file')

        # Check self.init_empty_orbitals
        if self.parameters.init_empty_orbitals != self.parameters.init_orbitals:
            raise NotImplementedError(f'The combination init_orbitals = {self.parameters.init_orbitals} '
                                      f'and init_empty_orbitals = {self.parameters.init_empty_orbitals} '
                                      'has not yet been implemented')

    def convert_kcp_to_supercell(self):
        # Multiply all extensive KCP settings by the appropriate prefactor
        prefactor = np.prod(self.kpoints.grid)
        for attr in ['nelec', 'nelup', 'neldw', 'nbnd', 'conv_thr', 'esic_conv_thr', 'tot_charge', 'tot_magnetization']:
            value = getattr(self.calculator_parameters['kcp'], attr, None)
            if value is not None:
                setattr(self.calculator_parameters['kcp'], attr, prefactor * value)

    def read_alphas_from_file(self, directory: Path = Path()):
        '''
        This routine reads in the contents of file_alpharef.txt and file_alpharef_empty.txt

        Since utils.read_alpha_file provides a flattened list of alphas so we must convert this
        to a nested list using convert_flat_alphas_for_kcp()
        '''

        flat_alphas = utils.read_alpha_file(directory)
        params = self.calculator_parameters['kcp']
        assert isinstance(params, KoopmansCPSettingsDict)
        alphas = calculators.convert_flat_alphas_for_kcp(flat_alphas, params)

        if self.parameters.spin_polarized:
            raise NotImplementedError('Need to check implementation')

        return alphas

    def _run(self) -> None:
        '''
        This function runs a KI/pKIPZ/KIPZ workflow from start to finish

        Running this function will generate several directories:
            init/                 -- the density and manifold initialization calculations
            calc_alpha/orbital_#/ -- calculations where we have fixed a particular orbital
                                     in order to calculate alpha
            final/                -- the final KI/KIPZ calculation
            postproc/             -- the unfolding and interpolation of the final band structure
        '''

        # Removing old directories
        if self.parameters.from_scratch:
            if not self._restart_from_old_ki:
                # if self._restart_from_old_ki we don't want to delete the directory containing
                # the KI calculation we're reading the manifold from, or the TMP files
                utils.system_call('rm -r init 2>/dev/null', False)
                utils.system_call(f'rm -r {self.calculator_parameters["kcp"].outdir} 2>/dev/null', False)
            utils.system_call('rm -r calc_alpha 2>/dev/null', False)
            utils.system_call('rm -r final 2>/dev/null', False)
            if self._redo_smooth_dft in [None, True]:
                utils.system_call('rm -r postproc 2>/dev/null', False)

        self.print('Initialization of density and variational orbitals', style='heading')
        self.perform_initialization()

        if self.parameters.from_scratch and not self._restart_from_old_ki \
                and self.parameters.fix_spin_contamination \
                and self.parameters.init_orbitals not in ['mlwfs', 'projwfs'] \
                and not (all(self.atoms.pbc) and self.parameters.init_orbitals == 'kohn-sham'):
            self.print('Copying the spin-up variational orbitals over to the spin-down channel')
            calc = self.calculations[-1]
            savedir = f'{calc.parameters.outdir}/{calc.parameters.prefix}_{calc.parameters.ndw}.save/K00001'
            utils.system_call(f'cp {savedir}/evc01.dat {savedir}/evc02.dat')

            assert isinstance(calc, calculators.KoopmansCPCalculator)
            if calc.has_empty_states():
                utils.system_call(f'cp {savedir}/evc0_empty1.dat {savedir}/evc0_empty2.dat')

        self.print('Calculating screening parameters', style='heading')
        if self.parameters.calculate_alpha:
            self.perform_alpha_calculations()
        else:
            self.print('Skipping calculation of screening parameters', end='')
            if len(self.bands.alpha_history) == 0:
                self.print('; reading values from file')
                self.bands.alphas = self.read_alphas_from_file()
            else:
                self.print()
            self.bands.print_history(indent=self.print_indent + 1)

        # Final calculation
        self.print(f'Final {self.parameters.functional.upper().replace("PK","pK")} calculation', style='heading')
        self.perform_final_calculations()

        # Postprocessing
        if all(self.atoms.pbc):
            if self.parameters.calculate_bands in [None, True] and self.projections and self.kpoints.path is not None:
                # Calculate interpolated band structure and DOS with UI
                from koopmans import workflows
                self.print(f'\nPostprocessing', style='heading')
                ui_workflow = workflows.UnfoldAndInterpolateWorkflow.fromparent(
                    self, redo_smooth_dft=self._redo_smooth_dft)
                ui_workflow.run(subdirectory='postproc')
            else:
                # Generate the DOS only
                dos = DOS(self.calculations[-1], width=self.plotting.degauss, npts=self.plotting.nstep + 1)
                self.calculations[-1].results['dos'] = dos

    def perform_initialization(self) -> None:
        # Import these here so that if these have been monkey-patched, we get the monkey-patched version
        from koopmans import workflows

        # The final calculation during the initialization, regardless of the workflow settings, should write to ndw = 51
        ndw_final = 51

        if self._restart_from_old_ki:
            self.print('Copying the density and orbitals from a pre-existing KI calculation')

            # Read the .cpi file to work out the value for ndw
            calc = calculators.KoopmansCPCalculator.fromfile('init/ki_init')

            # Move the old save directory to correspond to ndw_final, using pz_innerloop_init to work out where the
            # code will expect the tmp files to be
            old_savedir = Path(calc.parameters.outdir) / f'{calc.parameters.prefix}_{calc.parameters.ndw}.save'
            savedir = Path(f'{self.new_kcp_calculator("pz_innerloop_init").parameters.outdir}') / \
                f'{calc.parameters.prefix}_{ndw_final}.save'
            if not old_savedir.is_dir():
                raise ValueError(f'{old_savedir} does not exist; a previous '
                                 'and complete KI calculation is required '
                                 'to restart from an old KI calculation"')
            if savedir.is_dir():
                shutil.rmtree(savedir.as_posix())

            # Use the chdir construct in order to create the directory savedir if it does not already exist and is
            # nested
            with utils.chdir(savedir):
                utils.system_call(f'rsync -a {old_savedir}/ .')

            # Check that the files defining the variational orbitals exist
            savedir /= 'K00001'
            files_to_check = [Path('init/ki_init.cpo'), savedir / 'evc01.dat', savedir / 'evc02.dat']

            for ispin in range(2):
                if calc.has_empty_states(ispin):
                    files_to_check.append(savedir / f'evc0_empty{ispin + 1}.dat')

            for fname in files_to_check:
                if not fname.is_file():
                    raise ValueError(f'Could not find {fname}')

            if not calc.is_complete():
                raise ValueError('init/ki_init.cpo is incomplete so cannot be used '
                                 'to initialize the density and orbitals')

            self.calculations.append(calc)

        elif self.parameters.init_orbitals in ['mlwfs', 'projwfs'] or \
                (all(self.atoms.pbc) and self.parameters.init_orbitals == 'kohn-sham'):
            # Wannier functions using pw.x, wannier90.x and pw2wannier90.x (pw.x only for Kohn-Sham states)
            wannier_workflow = workflows.WannierizeWorkflow.fromparent(self)
            if wannier_workflow.parameters.calculate_bands:
                wannier_workflow.parameters.calculate_bands = \
                    not self.calculator_parameters['ui'].do_smooth_interpolation

            # Perform the wannierization workflow within the init directory
            wannier_workflow.run(subdirectory='init')

            # Now, convert the files over from w90 format to (k)cp format
            fold_workflow = workflows.FoldToSupercellWorkflow.fromparent(self)

            # Do this in the same directory as the wannierization
            fold_workflow.run(subdirectory='init')

            # Convert self.atoms to the supercell
            self.primitive_to_supercell()

            # We need a dummy calc before the real dft_init in order
            # to copy the previously calculated Wannier functions
            calc = self.new_kcp_calculator('dft_dummy')
            calc.directory = Path('init')
            self.run_calculator(calc, enforce_ss=False)

            # DFT restarting from Wannier functions (after copying the Wannier functions)
            calc = self.new_kcp_calculator('dft_init', restart_mode='restart',
                                           restart_from_wannier_pwscf=True, do_outerloop=True, ndw=ndw_final)
            calc.directory = Path('init')
            restart_dir = Path(f'{calc.parameters.outdir}/{calc.parameters.prefix}_{calc.parameters.ndr}.save/K00001')

            for filling in ['occ', 'emp']:
                for i_spin, spin in enumerate(['up', 'down']):
                    # Skip if we don't have wannier functions to copy over
                    if self.parameters.init_orbitals != 'kohn-sham':
                        if self.parameters.spin_polarized:
                            if self.projections.num_wann(occ=(filling == 'occ'), spin=spin) == 0:
                                continue
                        else:
                            if self.projections.num_wann(occ=(filling == 'occ'), spin=None) == 0:
                                continue

                    if self.parameters.init_orbitals == 'kohn-sham':
                        if filling == 'occ':
                            evcw_file = Path(f'init/wannier/ks2kcp/evc_occupied{i_spin + 1}.dat')
                        else:
                            evcw_file = Path(f'init/wannier/ks2kcp/evc0_empty{i_spin + 1}.dat')
                    elif self.parameters.spin_polarized:
                        evcw_file = Path(f'init/wannier/{filling}_{spin}/evcw.dat')
                    else:
                        evcw_file = Path(f'init/wannier/{filling}/evcw{i_spin + 1}.dat')

                    if filling == 'occ':
                        dest_file = restart_dir / f'evc_occupied{i_spin + 1}.dat'
                    else:
                        dest_file = restart_dir / f'evc0_empty{i_spin + 1}.dat'
                    if evcw_file.is_file():
                        shutil.copy(evcw_file, dest_file)
                    else:
                        raise OSError(f'Could not find {evcw_file}')

            self.run_calculator(calc, enforce_ss=False)

            # Check the consistency between the PW and CP band gaps
            pw_calc = [c for c in self.calculations if isinstance(
                c, calculators.PWCalculator) and c.parameters.calculation == 'nscf'][-1]
            pw_gap = pw_calc.results['lumo_energy'] - pw_calc.results['homo_energy']
            cp_gap = calc.results['lumo_energy'] - calc.results['homo_energy']
            if abs(pw_gap - cp_gap) > 2e-2 * pw_gap:
                raise ValueError(f'PW and CP band gaps are not consistent: {pw_gap} {cp_gap}')

            # The CP restarting from Wannier functions must be already converged
            Eini = calc.results['convergence']['filled'][0]['Etot']
            Efin = calc.results['energy']
            if abs(Efin - Eini) > 1e-6 * abs(Efin):
                raise ValueError(f'Too much difference between the initial and final CP energies: {Eini} {Efin}')

            # Add to the outdir of dft_init a link to the files containing the Wannier functions
            dst = Path(f'{calc.parameters.outdir}/{calc.parameters.prefix}_{calc.parameters.ndw}.save/K00001/')
            for file in ['evc_occupied1.dat', 'evc_occupied2.dat', 'evc0_empty1.dat', 'evc0_empty2.dat']:
                utils.symlink(f'{restart_dir}/{file}', dst, force=True)

        elif self.parameters.functional in ['ki', 'pkipz']:
            calc = self.new_kcp_calculator('dft_init')
            calc.directory = Path('init')
            self.run_calculator(calc, enforce_ss=self.parameters.fix_spin_contamination)

            # Use the KS eigenfunctions as better guesses for the variational orbitals
            self._overwrite_canonical_with_variational_orbitals(calc)

            if self.parameters.init_orbitals == 'kohn-sham':
                self._copy_most_recent_calc_to_ndw(ndw_final)
            elif self.parameters.init_orbitals == 'pz':
                calc = self.new_kcp_calculator('pz_innerloop_init', alphas=self.bands.alphas, ndw=ndw_final)
                calc.directory = Path('init')
                if self.calculations[-1].parameters.nelec == 2:
                    # If we only have two electrons, then the filled manifold is trivially invariant under unitary
                    # transformations. Furthermore, the PZ functional is invariant w.r.t. unitary rotations of the
                    # empty states. Thus in this instance we can skip the initialization of the manifold entirely
                    self.print('Skipping the optimisation of the variational orbitals since they are invariant under '
                               'unitary transformations')
                    self._copy_most_recent_calc_to_ndw(ndw_final)
                else:
                    self.run_calculator(calc)
            else:
                raise ValueError('Should not arrive here')

        elif self.parameters.functional == 'kipz':
            # DFT from scratch
            calc = self.new_kcp_calculator('dft_init')
            calc.directory = Path('init')
            self.run_calculator(calc, enforce_ss=self.parameters.fix_spin_contamination)

            if self.parameters.init_orbitals == 'kohn-sham':
                # Initialize the density with DFT and use the KS eigenfunctions as guesses for the variational orbitals
                self._overwrite_canonical_with_variational_orbitals(calc)
                self._copy_most_recent_calc_to_ndw(ndw_final)
            elif self.parameters.init_orbitals == 'pz':
                # PZ from DFT (generating PZ density and PZ orbitals)
                calc = self.new_kcp_calculator('pz_init', ndw=ndw_final)
                calc.directory = Path('init')
                self.run_calculator(calc)
            else:
                raise ValueError('Should not arrive here')

        else:
            raise ValueError("Should not arrive here; there must be an inconsistency between the above code and \
                             workflow.valid_settings")

        return

    def _copy_most_recent_calc_to_ndw(self, ndw):
        calc = self.calculations[-1]
        if calc.parameters.ndw != ndw:
            assert calc.is_complete(), 'Cannot copy results of a previous calculation that is not itself complete'
            save_prefix = f'{calc.parameters.outdir}/{calc.parameters.prefix}'
            utils.system_call(f'cp -r {save_prefix}_{calc.parameters.ndw}.save {save_prefix}_{ndw}.save')

    def _overwrite_canonical_with_variational_orbitals(self, calc: calculators.KoopmansCPCalculator) -> None:
        self.print('Overwriting the variational orbitals with Kohn-Sham orbitals')
        savedir = calc.parameters.outdir / f'{calc.parameters.prefix}_{calc.parameters.ndw}.save/K00001'
        for ispin in range(2):
            shutil.copy(savedir / f'evc{ispin + 1}.dat', savedir / f'evc0{ispin + 1}.dat')
            if calc.has_empty_states(ispin):
                shutil.copy(savedir / f'evc_empty{ispin + 1}.dat', savedir / f'evc0_empty{ispin + 1}.dat')

    def perform_alpha_calculations(self) -> None:
        # Set up directories
        Path('calc_alpha').mkdir(exist_ok=True)

        converged = False
        i_sc = 0

        alpha_indep_calcs: List[calculators.KoopmansCPCalculator] = []

        while not converged and i_sc < self.parameters.n_max_sc_steps:
            i_sc += 1

            # Setting up directories
            iteration_directory = Path('calc_alpha')
            outdir = self.calculator_parameters['kcp'].outdir.name
            outdir = Path.cwd() / iteration_directory / outdir

            if not outdir.is_dir():
                outdir.mkdir()

            if self.parameters.n_max_sc_steps > 1:
                self.print('SC iteration {}'.format(i_sc), style='subheading')
                iteration_directory /= f'iteration_{i_sc}'
                if not iteration_directory.is_dir():
                    iteration_directory.mkdir()

            # Do a KI/KIPZ calculation with the updated alpha values
            restart_from_wannier_pwscf = True if self.parameters.init_orbitals in [
                'mlwfs', 'projwfs'] and not self._restart_from_old_ki and i_sc == 1 else None
            if self.parameters.task in ['trajectory', 'convergence_ml'] and self.ml.input_data_for_ml_model == 'orbital_density':
                print_real_space_density = True
            else:
                print_real_space_density = False
            trial_calc = self.new_kcp_calculator(calc_presets=self.parameters.functional.replace('pkipz', 'ki'),
                                                 print_real_space_density=print_real_space_density,
                                                 alphas=self.bands.alphas,
                                                 restart_from_wannier_pwscf=restart_from_wannier_pwscf)
            trial_calc.directory = iteration_directory

            if i_sc == 1:
                if self.parameters.functional == 'kipz' and not all(self.atoms.pbc):
                    # For the first KIPZ trial calculation, do the innerloop
                    trial_calc.parameters.do_innerloop = True
            else:
                # For later SC loops, read in the matching calculation from the
                # previous loop rather than the initialization calculations
                trial_calc.parameters.ndr = trial_calc.parameters.ndw

            # Run the calculation and store the result. Note that we only need to continue
            # enforcing the spin symmetry if the density will change
            self.run_calculator(trial_calc, enforce_ss=self.parameters.fix_spin_contamination and i_sc > 1)

            alpha_dep_calcs = [trial_calc]

            # Update the bands' self-Hartree and energies (assuming spin-symmetry)
            self.bands.self_hartrees = trial_calc.results['orbital_data']['self-Hartree']

            # Group the bands
            self.bands.assign_groups(allow_reassignment=True)

            skipped_orbitals = []
            first_band_of_each_channel = [self.bands.get(spin=spin)[0] for spin in range(2)]

            # Initialize the ML-model
            if self.ml.use_ml:
                from koopmans.workflows import MLFittingWorkflow
                mlfit = MLFittingWorkflow.fromparent(self, calc_that_produced_orbital_densities=trial_calc)
                mlfit.run()

            # Loop over removing/adding an electron from/to each orbital
            for band in self.bands:
                # For a KI calculation with only filled bands, we don't have any further calculations to
                # do, so in this case don't print any headings
                print_headings = self.parameters.functional != 'ki' \
                    or any([not b.filled for b in self.bands]) or i_sc == 1

                if self.parameters.spin_polarized and band in first_band_of_each_channel:
                    self.print(f'Spin {band.spin + 1}', style='subheading')

                # Working out what to print for the orbital heading (grouping skipped bands together)
                if band in self.bands.to_solve or band == self.bands.get(spin=band.spin)[-1]:
                    if band not in self.bands.to_solve and (self.parameters.spin_polarized or band.spin == 0):
                        skipped_orbitals.append(band.index)
                    if len(skipped_orbitals) > 0:
                        if len(skipped_orbitals) == 1:
                            if print_headings:
                                self.print(f'Orbital {skipped_orbitals[0]}', style='subheading')
                            else:
                                orb_range = f'{skipped_orbitals[0]}-{skipped_orbitals[-1]}'
                                if print_headings:
                                    self.print(f'Orbitals {orb_range}', style='subheading')
                            if print_headings:
                                self.print(f'Skipping; will use the screening parameter of an equivalent orbital')
                            skipped_orbitals = []
                    if band not in self.bands.to_solve:
                        continue
                elif not self.parameters.spin_polarized and band.spin == 1:
                    # In this case, skip over the bands entirely and don't include it in the printout about which
                    # bands we've skipped
                    continue
                else:
                    # Skip the bands which can copy the screening parameter from another
                    # calculation in the same orbital group
                    skipped_orbitals.append(band.index)
                    continue
                if print_headings:
                    self.print(f'Orbital {band.index}', style='subheading')

                # Set up directories
                if self.parameters.spin_polarized:
                    directory = Path(f'{iteration_directory}/spin_{band.spin + 1}/orbital_{band.index}')
                    outdir_band = outdir / f'spin_{band.spin + 1}/orbital_{band.index}'
                else:
                    directory = Path(f'{iteration_directory}/orbital_{band.index}')
                    outdir_band = outdir / f'orbital_{band.index}'
                if not directory.is_dir():
                    directory.mkdir(parents=True)

                # Link tmp files from band-independent calculations
                if not outdir_band.is_dir():
                    outdir_band.mkdir(parents=True)

                    utils.symlink(f'{trial_calc.parameters.outdir}/*.save', outdir_band)

                # Don't repeat if this particular alpha_i was converged
                if i_sc > 1 and abs(band.error) < self.parameters.alpha_conv_thr:
                    self.print(f'Skipping band {band.index} since this alpha is already converged')
                    # if self.parameters.from_scratch:
                    for b in self.bands:
                        if b == band or (band.group is not None and b.group == band.group):
                            b.alpha = band.alpha
                            b.error = band.error
                    continue

                # When we write/update the alpharef files in the work directory
                # make sure to include the fixed band alpha in file_alpharef.txt
                # rather than file_alpharef_empty.txt
                if band.filled:
                    index_empty_to_save = None
                else:
                    index_empty_to_save = band.index - self.bands.num(filled=True, spin=band.spin)
                    if self.parameters.spin_polarized and band.spin == 1:
                        index_empty_to_save += self.bands.num(filled=False, spin=0)

                # Make ML-prediction and decide whether we want to use this prediction
                if self.ml.use_ml:
                    alpha_predicted = mlfit.predict(band)
                    # Whether to use the ML-prediction
                    use_prediction = mlfit.use_prediction()
                if not self.ml.use_ml or not (use_prediction or self.ml.alphas_from_file):
                    self.perform_fixed_band_calculations(band, trial_calc, i_sc, alpha_dep_calcs, index_empty_to_save,
                                                         outdir_band, directory, alpha_indep_calcs)

                if self.ml.use_ml and use_prediction:
                    alpha = alpha_predicted
                    error = 0.0  # set the error for the predicted alphas to 0.0, because we don't want to make another
                    # scf-step because of predicted alphas
                else:
                    if self.ml.use_ml and self.ml.alphas_from_file:
                        # Dummy calculation to circumvent the fixed-band-calculation for debugging
                        alpha, error = mlfit.get_alpha_from_file_for_debugging(band)
                    else:
                        # Calculate an updated alpha and a measure of the error
                        # E(N) - E_i(N - 1) - lambda^alpha_ii(1)     (filled)
                        # E_i(N + 1) - E(N) - lambda^alpha_ii(0)     (empty)
                        #
                        # Note that we can do this even from calculations that have been skipped because
                        # we read in all the requisite information from the output files and .pkl files
                        # that do not get overwritten

                        calcs = [c for c in alpha_dep_calcs + alpha_indep_calcs if c.fixed_band == band]

                        alpha, error = self.calculate_alpha_from_list_of_calcs(
                            calcs, trial_calc, band, filled=band.filled)

                for b in self.bands:
                    if b == band or (b.group is not None and b.group == band.group):
                        b.alpha = alpha
                        b.error = error

                # add alpha to training data
                if self.ml.use_ml and not use_prediction:
                    mlfit.print_error_of_single_orbital(alpha_predicted, alpha, indent=self.print_indent+2)
                    mlfit.add_training_data(band)
                    # if the user wants to train on the fly, train the model after the calculation of each orbital
                    if self.ml.train_on_the_fly:
                        mlfit.train()

            self.bands.print_history(indent=self.print_indent + 1)

            converged = all([abs(b.error) < 1e-3 for b in self.bands])

            if self.ml.use_ml and not any(mlfit.use_predictions):
                # if the user don't wants to train on the fly, train the model at the end of each snapshot
                if not self.ml.train_on_the_fly:
                    mlfit.train()
                # Print summary of all predictions
                mlfit.print_error_of_all_orbitals(indent=self.print_indent + 1)

        if self.parameters.functional == 'ki' and self.bands.num(filled=False):
            # For this case the screening parameters are guaranteed to converge instantly
            if self.parameters.n_max_sc_steps == 1:
                # Print the "converged" message rather than the "determined but not necessarily converged" message
                converged = True
            else:
                # Do the subsequent loop
                utils.warn('The screening parameters for a KI calculation with no empty states will converge '
                           'instantly; to save computational time set n_max_sc_steps == 1')
        if converged:
            self.print('Screening parameters have been converged')
        else:
            self.print('Screening parameters have been determined but are not necessarily converged')

    def perform_fixed_band_calculations(self, band, trial_calc, i_sc, alpha_dep_calcs, index_empty_to_save, outdir_band, directory, alpha_indep_calcs) -> None:
        # Perform the fixed-band-dependent calculations
        if self.parameters.functional in ['ki', 'pkipz']:
            if band.filled:
                calc_types = ['dft_n-1']
            else:
                calc_types = ['pz_print', 'dft_n+1_dummy', 'dft_n+1']
        else:
            if band.filled:
                calc_types = ['kipz_n-1']
            else:
                calc_types = ['kipz_print', 'dft_n+1_dummy', 'kipz_n+1']

        for calc_type in calc_types:
            if self.parameters.functional in ['ki', 'pkipz']:
                # The calculations whose results change with alpha are...
                #  - the KI calculations
                #  - DFT calculations on empty variational orbitals
                # We don't need to redo any of the others
                if not trial_calc.has_empty_states() or band.filled:
                    if i_sc > 1 and 'ki' not in calc_type:
                        continue
            else:
                # No need to repeat the dummy calculation; all other
                # calculations are dependent on the screening parameters so
                # will need updating at each step
                if i_sc > 1 and calc_type == 'dft_n+1_dummy':
                    continue

            if 'print' in calc_type:
                # Note that the 'print' calculations for empty bands do not
                # in fact involve the fixing of that band (and thus for the
                # 'fixed' band the corresponding alpha should be in
                # file_alpharef_empty.txt)
                alphas = self.bands.alphas
                filling = self.bands.filling
            elif not band.filled:
                # In the case of empty orbitals, we gain an extra orbital in
                # the spin-up channel, so we explicitly construct both spin
                # channels for "alphas" and "filling"
                alphas = self.bands.alphas
                alphas[band.spin].append(alphas[band.spin][-1])
                filling = self.bands.filling
                filling[band.spin][band.index - 1] = True
                filling[band.spin].append(False)
            else:
                alphas = self.bands.alphas
                filling = self.bands.filling

            # Work out the index of the band that is fixed (noting that we will be throwing away all empty
            # bands)
            fixed_band = min(band.index, self.bands.num(filled=True, spin=band.spin) + 1)
            if self.parameters.spin_polarized and band.spin == 1:
                fixed_band += self.bands.num(filled=True, spin=0)

            # Set up calculator
            calc = self.new_kcp_calculator(calc_type, alphas=alphas, filling=filling, fixed_band=fixed_band,
                                           index_empty_to_save=index_empty_to_save, outdir=outdir_band,
                                           add_to_spin_up=(band.spin == 0))
            calc.directory = directory

            # Run kcp.x
            self.run_calculator(calc)

            # Store the band that we've perturbed as calc.fixed_band. Note that we can't use
            # calc.parameters.fixed_band to keep track of which band we held fixed, because for empty
            # orbitals, calc.parameters.fixed_band is always set to the LUMO but in reality we're fixing
            # the band corresponding # to index_empty_to_save from an earlier calculation
            calc.fixed_band = band

            # Store the result
            # We store the results in one of two lists: alpha_indep_calcs and
            # alpha_dep_calcs. The latter is overwritten at each new self-
            # consistency loop.
            if 'ki' in calc_type and 'print' not in calc_type:
                alpha_dep_calcs.append(calc)
            elif 'dft' in calc_type and 'dummy' not in calc_type:
                if self.parameters.functional in ['ki', 'pkipz']:
                    # For KI, the results of the DFT calculations are typically independent of alpha so we
                    # store these in a list that is never overwritten

                    # The exception to this are KI calculations on empty states. When we update alpha, the
                    # empty manifold changes, which in turn affects the lambda values
                    if trial_calc.has_empty_states() and not band.filled:
                        alpha_dep_calcs.append(calc)
                    else:
                        alpha_indep_calcs.append(calc)
                else:
                    # For KIPZ, the DFT calculations are dependent on alpha via
                    # the definition of the variational orbitals. We only want to
                    # store the calculations that used the most recent value of alpha

                    alpha_dep_calcs.append(calc)

            # Copying of evcfixed_empty.dat to evc_occupied.dat
            if calc_type in ['pz_print', 'kipz_print']:
                evcempty_dir = outdir_band / f'{calc.parameters.prefix}_{calc.parameters.ndw}.save/K00001/'
            elif calc_type == 'dft_n+1_dummy':
                evcocc_dir = outdir_band / f'{calc.parameters.prefix}_{calc.parameters.ndr}.save/K00001/'
                for i_spin in range(1, 3):
                    src = evcempty_dir / f'evcfixed_empty{i_spin}.dat'
                    dest = evcocc_dir / f'evc_occupied{i_spin}.dat'
                    if src.is_file():
                        shutil.copy(src, dest)
                    else:
                        raise OSError(f'Could not find {src}')

    def perform_final_calculations(self) -> None:

        directory = Path('final')
        if not directory.is_dir():
            directory.mkdir()

        if self.parameters.functional == 'pkipz':
            final_calc_types = ['ki', 'pkipz']
        else:
            final_calc_types = [self.parameters.functional]

        for final_calc_type in final_calc_types:

            final_calc_type += '_final'

            # For pKIPZ, the appropriate ndr can change but it is always ndw of the previous
            # KI calculation
            if final_calc_type == 'pkipz_final':
                ndr = [c.parameters.ndw for c in self.calculations if c.prefix in [
                    'ki', 'ki_final'] and hasattr(c.parameters, 'ndw')][-1]
                calc = self.new_kcp_calculator(final_calc_type, ndr=ndr, write_hr=True)
            else:
                calc = self.new_kcp_calculator(final_calc_type, write_hr=True)
                if self.parameters.functional == 'ki' and self.parameters.init_orbitals in ['mlwfs', 'projwfs'] \
                        and not self.parameters.calculate_alpha:
                    calc.parameters.restart_from_wannier_pwscf = True

            calc.directory = directory

            self.run_calculator(calc)

    def new_kcp_calculator(self, calc_presets: str = 'dft_init',
                           alphas: Optional[List[List[float]]] = None,
                           filling: Optional[List[List[bool]]] = None,
                           add_to_spin_up: bool = True,
                           **kwargs) -> calculators.KoopmansCPCalculator:
        """

        Generates a new KCP calculator based on the self.calculator_parameters["kcp"]
        parameters, modifying the appropriate settings to match the
        chosen calc_presets, and altering any Quantum Espresso keywords
        specified as kwargs

        Arguments:

            calc_presets
                The set of preset values to use; must be one of the following strings:

                Initialization
                'dft_init'            DFT calculation from scratch
                'pz_init'             PZ calculation starting from DFT restart
                'pz_innerloop_init'   PZ calculation starting from DFT restart (innerloop only)
                'dft_dummy'           DFT dummy calculation that generate store files
                                      for a periodic calculation restarting from Wannier functions

                Trial calculations
                'ki'     KI calculation with N electrons and empty bands if specified
                'kipz'   As above, but for KIPZ

                For calculating alpha_i for filled orbitals.
                'dft_n-1'       DFT calculation with N-1 electrons via fixed_state
                'kipz_n-1'      KIPZ calculation with N-1 electrons via fixed_state

                For calculating alpha_i for empty orbitals
                'pz_print'         PZ calculation that generates evcfixed_empty.dat file
                'kipz_print'       KIPZ calculation that generates evcfixed_empty.dat file
                'dft_n+1_dummy'    DFT dummy calculation that generates store files of
                                   the correct dimensions
                'dft_n+1'          DFT calculation with N+1 electrons
                'kipz_n+1'         KIPZ calculation with N+1 electrons

                Note that when calculating alpha_i, all empty states (bar orbital_i if it is empty) are removed
                and the convergence criteria are loosened

                Final calculation
                'ki_final'    Final KI calculation with N electrons and empty bands if specified
                'kipz_final'  As above, but for KIPZ
                'pkipz_final' A KIPZ calculation that leaves the manifold unchanged (for the
                              purposes of performing KIPZ on top of the KI manifold)

            alphas
                an array of screening parameters

           **kwargs accepts any Quantum Espresso keywords as an argument, and will
                    apply these options to the returned ASE calculator

        Returns: a new KCP calculator object

        """

        # By default, use the last row in the alpha table for the screening parameters
        if alphas is None:
            alphas = self.bands.alphas

        # Generate a new kcp calculator copied from the master calculator
        calc: calculators.KoopmansCPCalculator = self.new_calculator('kcp', alphas=alphas, filling=filling)

        # Set up read/write indexes
        if calc_presets in ['dft_init', 'dft_dummy']:
            ndr = 50
            ndw = 50
        elif calc_presets in ['pz_init', 'pz_innerloop_init']:
            ndr = 50
            ndw = 51
        elif calc_presets in ['ki', 'kipz']:
            ndr = 51
            ndw = 60
        elif calc_presets in ['dft_n-1', 'kipz_n-1']:
            ndr = 60
            ndw = 63
        elif calc_presets in ['pz_print', 'kipz_print']:
            ndr = 60
            ndw = 64
        elif calc_presets == 'dft_n+1_dummy':
            ndr = 65
            ndw = 65
        elif calc_presets in ['dft_n+1', 'kipz_n+1']:
            ndr = 65
            ndw = 68
        elif calc_presets in ['ki_final', 'kipz_final']:
            if self.parameters.calculate_alpha:
                ndr = 60
            else:
                ndr = 51
            ndw = 70
        elif calc_presets in ['pkipz_final']:
            ndr = 70
            ndw = 71
        else:
            raise ValueError('Invalid calc_presets "{}"'.format(calc_presets))

        # KCP options
        # control
        calc.prefix = calc_presets
        calc.parameters.ndw = ndw
        calc.parameters.ndr = ndr
        if calc.prefix in ['dft_init', 'dft_dummy', 'dft_n+1_dummy']:
            calc.parameters.restart_mode = 'from_scratch'
        else:
            calc.parameters.restart_mode = 'restart'

        # system
        if 'pz' in calc.prefix or 'ki' in calc.prefix:
            calc.parameters.do_orbdep = True
            calc.parameters.do_bare_eigs = True
        else:
            calc.parameters.do_orbdep = False
        if calc.prefix in ['dft_init', 'pz_init', 'pz_innerloop_init', 'dft_dummy',
                           'ki', 'kipz', 'pz_print', 'kipz_print', 'dft_n+1_dummy',
                           'ki_final', 'kipz_final', 'pkipz_final']:
            calc.parameters.fixed_state = False
        else:
            calc.parameters.fixed_state = True
            if '-1' in calc.prefix:
                calc.parameters.f_cutoff = 1e-5
            else:
                calc.parameters.f_cutoff = 1.0
        if 'n+1' in calc.prefix:
            calc.parameters.nelec += 1
            if add_to_spin_up:
                calc.parameters.nelup += 1
            else:
                calc.parameters.neldw += 1
            if 'dummy' not in calc.prefix:
                calc.parameters.restart_from_wannier_pwscf = True

        # electrons
        # For all calculations calculating alpha, remove the empty states and
        # increase the energy thresholds
        if not any([s in calc.prefix for s in ['init', 'print', 'final']]) and \
           calc.prefix not in ['ki', 'kipz', 'dft_dummy']:
            calc.parameters.nbnd = None
            calc.parameters.conv_thr *= 100
            calc.parameters.esic_conv_thr *= 100

        # For the dft_dummy calculation, we don't need empty states because these will be overwritten by the w90
        # wavefunctions
        if calc.prefix == 'dft_dummy':
            calc.parameters.nbnd = None

        if all(self.atoms.pbc) and not any([s == calc.prefix for s in ['dft_init', 'dft_n-1', 'dft_n+1',
                                                                       'kipz', 'kipz_n-1', 'kipz_n+1']]):
            calc.parameters.do_outerloop = False
            calc.parameters.do_innerloop = False
        elif any([s in calc.prefix for s in ['frozen', 'dummy', 'print', 'innerloop']]) or calc.prefix == 'pkipz_final':
            calc.parameters.do_outerloop = False
            if calc.has_empty_states():
                calc.parameters.do_outerloop_empty = False
        elif calc.prefix in ['ki', 'ki_final']:
            calc.parameters.do_outerloop = False
            if calc.has_empty_states():
                if self.parameters.init_empty_orbitals == 'pz':
                    calc.parameters.do_outerloop_empty = True
                else:
                    calc.parameters.do_outerloop_empty = False
        else:
            calc.parameters.do_outerloop = True
            if calc.has_empty_states():
                calc.parameters.do_outerloop_empty = True

        if calc.parameters.maxiter is None and calc.parameters.do_outerloop:
            calc.parameters.maxiter = 300
        if calc.parameters.empty_states_maxstep is None and calc.parameters.do_outerloop_empty:
            calc.parameters.empty_states_maxstep = 300

        # No empty states minimization in the solids workflow for the moment
        if all(self.atoms.pbc) and calc.has_empty_states():
            calc.parameters.do_outerloop_empty = False
            calc.parameters.do_innerloop_empty = False

        # nksic
        if calc.parameters.do_orbdep:
            calc.parameters.odd_nkscalfact = True
            calc.parameters.odd_nkscalfact_empty = True
        calc.parameters.do_innerloop_cg = True
        if calc.prefix[:2] == 'pz' and 'print' not in calc.prefix:
            calc.parameters.do_innerloop = True
        else:
            calc.parameters.do_innerloop = False
        if calc.has_empty_states():
            calc.parameters.do_innerloop_empty = False
        if 'kipz' in calc.prefix:
            calc.parameters.which_orbdep = 'nkipz'
        elif 'pz' in calc.prefix:
            calc.parameters.which_orbdep = 'pz'
        elif 'ki' in calc.prefix:
            calc.parameters.which_orbdep = 'nki'
        if 'print' in calc.prefix:
            calc.parameters.print_wfc_anion = True

        if self.parameters.mt_correction:
            calc.parameters.which_compensation = 'tcc'
        else:
            calc.parameters.which_compensation = 'none'

        # If we are using frozen orbitals, we override the above logic and freeze the variational orbitals
        # post-initialization
        if self.parameters.frozen_orbitals and 'init' not in calc.prefix and not any([s == calc.prefix for s in
                                                                                      ['dft_n-1', 'dft_n+1', 'kipz_n-1',
                                                                                       'kipz_n+1']]):
            calc.parameters.do_outerloop = False
            calc.parameters.do_innerloop = False
            if calc.has_empty_states():
                calc.parameters.do_outerloop_empty = False
                calc.parameters.do_innerloop_empty = False

        # Handle any keywords provided by kwargs
        # Note that since this is performed after the above logic this can (deliberately
        # or accidentally) overwrite the above settings
        calc.parameters.update(**kwargs)

        # Sanity checking
        if calc.parameters.print_wfc_anion and calc.parameters.index_empty_to_save is None:
            raise ValueError('Error: print_wfc_anion is set to true but you have not selected '
                             'an index_empty_to_save. Provide this as an argument to new_cp_calculator')

        # don't print QC in some cases
        if 'dummy' in calc.prefix or calc.prefix[-2:] == '+1':
            calc.skip_qc = True

        return calc

    def calculate_alpha_from_list_of_calcs(self,
                                           calcs: List[calculators.KoopmansCPCalculator],
                                           trial_calc: calculators.KoopmansCPCalculator,
                                           band: Band,
                                           filled: bool = True) -> Tuple[float, float]:
        '''

        Calculates alpha via equation 10 of Nguyen et. al (2018) 10.1103/PhysRevX.8.021051
        If the band is filled, use s = 1; if the band is empty, use s = 0

        Arguments:
            calcs          -- a list of selected calculations from which to calculate alpha
            trial_calc     -- the N-electron Koopmans calculation
            filled         -- True if the orbital for which we're calculating alpha is filled

        '''

        # Extract the energy difference Delta E
        if self.parameters.functional == 'kipz':
            if filled:
                # KIPZ N-1
                [kipz_m1_calc] = [c for c in calcs if c.parameters.which_orbdep == 'nkipz'
                                  and c.parameters.do_orbdep and c.parameters.f_cutoff < 0.0001]
                kipz_m1 = kipz_m1_calc.results
                charge = 1 - kipz_m1_calc.parameters.f_cutoff

                dE = trial_calc.results['energy'] - kipz_m1['energy']
                mp1 = kipz_m1['mp1_energy']
                mp2 = kipz_m1['mp2_energy']

            else:
                # KIPZ N+1
                [kipz_p1_calc] = [c for c in calcs if c.parameters.which_orbdep == 'nkipz'
                                  and c.parameters.do_orbdep and c.parameters.f_cutoff == 1.0]
                kipz_p1 = kipz_p1_calc.results
                charge = - kipz_p1_calc.parameters.f_cutoff

                dE = kipz_p1['energy'] - trial_calc.results['energy']
                mp1 = kipz_p1['mp1_energy']
                mp2 = kipz_p1['mp2_energy']

        else:
            # self.functional in ['ki', 'pkipz']
            if filled:
                # DFT N-1
                [dft_m1_calc] = [c for c in calcs if not c.parameters.do_orbdep
                                 and c.parameters.restart_mode == 'restart' and c.parameters.f_cutoff < 0.0001]
                dft_m1 = dft_m1_calc.results
                charge = 1 - dft_m1_calc.parameters.f_cutoff

                dE = trial_calc.results['energy'] - dft_m1['energy']
                mp1 = dft_m1['mp1_energy']
                mp2 = dft_m1['mp2_energy']

            else:
                # DFT N+1
                [dft_p1_calc] = [c for c in calcs if not c.parameters.do_orbdep
                                 and c.parameters.restart_mode == 'restart' and c.parameters.f_cutoff == 1.0]
                dft_p1 = dft_p1_calc.results
                charge = - dft_p1_calc.parameters.f_cutoff

                dE = dft_p1['energy'] - trial_calc.results['energy']
                mp1 = dft_p1['mp1_energy']
                mp2 = dft_p1['mp2_energy']

        # Extract lambda from the base calculator
        assert band.index is not None
        iband = band.index - 1  # converting from 1-indexing to 0-indexing
        lambda_a = trial_calc.results['lambda'][band.spin][iband, iband].real
        lambda_0 = trial_calc.results['bare lambda'][band.spin][iband, iband].real

        # Obtaining alpha
        if (trial_calc.parameters.odd_nkscalfact and filled) \
                or (trial_calc.parameters.odd_nkscalfact_empty and not filled):
            alpha_guess = trial_calc.alphas[band.spin][iband]
        else:
            alpha_guess = trial_calc.parameters.nkscalfact

        # Checking Makov-Payne correction energies and applying them (if needed)
        if self.parameters.mp_correction:
            if mp1 is None:
                raise ValueError('Could not find 1st order Makov-Payne energy')
            if mp2 is None:
                # utils.warn('Could not find 2nd order Makov-Payne energy; applying first order only')
                mp_energy = mp1
            else:
                mp_energy = mp1 + mp2

            dE -= np.sign(charge) * mp_energy / self.parameters.eps_inf

        alpha = alpha_guess * (dE - lambda_0) / (lambda_a - lambda_0)

        # The error is lambda^alpha(1) - lambda^alpha_i(1)
        error = dE - lambda_a

        return alpha, error
