"""

Workflow module for koopmans, containing the workflow for performing KI and KIPZ calculations

Written by Edward Linscott Jan 2020
Split off from workflow.py Oct 2020

"""

import os
import copy
import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union
import pandas as pd
from ase.dft.dos import DOS
from ase.spectrum.band_structure import BandStructure
from koopmans import utils
from koopmans.bands import Bands
from koopmans import calculators
from ._generic import Workflow


class KoopmansDSCFWorkflow(Workflow):

    def __init__(self, workflow_settings: Dict[str, Any], calcs_dct: Dict[str, calculators.ExtendedCalculator]) -> None:
        super().__init__(workflow_settings, calcs_dct)

        if 'kcp' not in self.master_calcs:
            raise ValueError(
                'Performing a KC calculation requires a "kcp" block in the .json input file')
        kcp_calc = self.master_calcs['kcp']

        # If periodic, convert the kcp calculation into a Γ-only supercell calculation
        if self.periodic:
            kpts = self.master_calcs['kcp'].parameters.kpts
            self.master_calcs['kcp'].transform_to_supercell(np.diag(kpts))
            nocc = self.master_calcs['w90_occ'].parameters.num_wann
            nemp = self.master_calcs['w90_emp'].parameters.num_wann
            if not self.orbital_groups:
                self.orbital_groups = range(0, nocc + nemp)
            self.orbital_groups = [i for _ in range(np.prod(kpts)) for i in self.orbital_groups[:nocc]] + \
                                  [i for _ in range(np.prod(kpts)) for i in self.orbital_groups[nocc:]]

            # Check the number of empty states has been correctly configured
            w90_emp_calc = self.master_calcs['w90_emp']
            expected_empty_states_nbnd = w90_emp_calc.parameters.num_wann * np.prod(kcp_calc.parameters.kpts)
            if kcp_calc.parameters.empty_states_nbnd == 0:
                # 0 is the default value
                kcp_calc.parameters.empty_states_nbnd = expected_empty_states_nbnd
            elif kcp_calc.parameters.empty_states_nbnd != expected_empty_states_nbnd:
                raise ValueError('kcp empty_states_nbnd and wannier90 num_wann (emp) are inconsistent')

        # Initialise the bands object
        filling = kcp_calc.filling[0]
        self.bands = Bands(n_bands=len(filling), filling=filling, groups=self.orbital_groups,
                           self_hartree_tol=self.orbital_groups_self_hartree_tol)
        if self.alpha_from_file:
            # Reading alpha values from file
            self.bands.alphas = self.read_alphas_from_file()
        else:
            # Initialising alpha with a guess
            self.bands.alphas = self.alpha_guess

        # Raise errors if any UI keywords are provided but will be overwritten by the workflow
        for ui_keyword in ['kc_ham_file', 'w90_seedname', 'alat_sc', 'sc_dim', 'dft_ham_file', 'dft_smooth_ham_file']:
            for ui_kind in ['occ', 'emp']:
                value = getattr(self.master_calcs[f'ui_{ui_kind}'].parameters, ui_keyword)
                [default_value] = [setting.default for setting in self.master_calcs['ui'].parameters.valid
                                   if setting.name == ui_keyword]
                if value != default_value:
                    raise ValueError(f'UI keyword {ui_keyword} has been set in the input file, but this will be '
                                     'automatically set by the Koopmans workflow. Remove this keyword from the input '
                                     'file')

        # Initialise self.init_empty_orbitals has not been set
        if self.init_empty_orbitals == 'same':
            self.init_empty_orbitals = self.init_orbitals
        if self.init_empty_orbitals != self.init_orbitals:
            raise NotImplementedError(f'The combination init_orbitals = {self.init_orbitals} and init_empty_orbitals '
                                      f'= {self.init_empty_orbitals} has not yet been implemented')

    def read_alphas_from_file(self, directory='.'):
        '''
        This routine reads in the contents of file_alpharef.txt and file_alpharef_empty.txt and
        stores the result in self.bands.alphas

        Note that utils.read_alpha_file provides a flattened list of alphas so we must exclude
        duplicates if calc.nspin == 2
        '''

        params = self.master_calcs['kcp'].parameters

        if params.nspin == 2:
            i_alphas = list(range(0, params.nelec // 2)) + \
                list((range(params.nelec, params.nelec + params.empty_states_nbnd)))
        else:
            i_alphas = list((range(0, params.nelec // 2 + params.empty_states_nbnd)))

        alphas = utils.read_alpha_file(directory)
        return [a for i, a in enumerate(alphas) if i in i_alphas]

    def run(self) -> None:
        '''
        This function runs the KI/KIPZ workflow from start to finish

        Running this function will generate a number of files:
            init/                 -- the density and manifold initialisation calculations
            calc_alpha/orbital_#/ -- calculations where we have fixed a particular orbital
                                     in order to calculate alpha
            final/                -- the final KI/KIPZ calculation
            alphas.pkl            -- a python pickle file containing the alpha values
            errors.pkl            -- a python pickle file containing the errors in the alpha
                                     values
            tab_alpha_values.tex  -- a latex table of the alpha values
        '''

        # Removing old directories
        if self.from_scratch:
            if self.init_orbitals != 'from old ki':
                # if self.init_orbitals == "from old ki" we don't want to delete the directory containing
                # the KI calculation we're reading the manifold from, or the TMP files
                utils.system_call('rm -r init 2>/dev/null', False)
                utils.system_call(f'rm -r {self.master_calcs["kcp"].parameters.outdir} 2>/dev/null', False)
            utils.system_call('rm -r calc_alpha 2>/dev/null', False)
            utils.system_call('rm -r final 2>/dev/null', False)
            if getattr(self, 'redo_preexisting_smooth_dft_calcs', True):
                utils.system_call('rm -r postproc 2>/dev/null', False)

        self.print('Initialisation of density and variational orbitals', style='heading')
        self.perform_initialisation()

        if self.from_scratch and self.init_orbitals != 'from old ki' and self.enforce_spin_symmetry \
                and self.init_orbitals not in ['mlwfs', 'projwfs']:
            self.print('Copying the spin-up variational orbitals over to the spin-down channel')
            calc = self.all_calcs[-1]
            savedir = f'{calc.parameters.outdir}/{calc.parameters.prefix}_{calc.parameters.ndw}.save/K00001'
            utils.system_call(f'cp {savedir}/evc01.dat {savedir}/evc02.dat')
            if calc.parameters.empty_states_nbnd is not None and calc.parameters.empty_states_nbnd > 0:
                utils.system_call(
                    f'cp {savedir}/evc0_empty1.dat {savedir}/evc0_empty2.dat')

        self.print('Calculating screening parameters', style='heading')
        if self.calculate_alpha:
            self.perform_alpha_calculations()
        else:
            self.print('Skipping calculation of screening parameters', end='')
            if len(self.bands.alpha_history) == 0:
                self.print('; reading values from file')
                self.bands.alphas = self.read_alphas_from_file()
            else:
                self.print()
            self.bands.print_history()

        # Final calculation
        self.print(f'Final {self.functional.upper().replace("PK","pK")} calculation', style='heading')
        self.perform_final_calculations()
        # Print out additional quality control data for final calculation
        if self.print_qc:
            calc = self.all_calcs[-1]
            for b in self.bands.to_solve:
                self.print_qc_keyval(f'alpha({b.index})', b.alpha, calc)
            for isp, orbs_self_hartree in enumerate(calc.results['orbital_data']['self-Hartree']):
                for i, orb_sh in enumerate(orbs_self_hartree):
                    self.print_qc_keyval(f'orbital_self_Hartree(orb={i+1},sigma={isp+1})', orb_sh, calc)

        # Postprocessing
        if self.periodic and self.master_calcs['ui'].parameters.kpath is not None:
            self.print(f'\nPostprocessing', style='heading')
            self.perform_postprocessing()

    def perform_initialisation(self) -> None:
        # Import these here so that if these have been monkey-patched, we get the monkey-patched version
        from koopmans.workflows import WannierizeWorkflow, FoldToSupercellWorkflow

        # The final calculation during the initialisation, regardless of the workflow settings, should write to ndw = 51
        ndw_final = 51

        if self.init_orbitals in ['mlwfs', 'projwfs']:
            # Wannier functions using pw.x, wannier90.x and pw2wannier90.x
            wannier_workflow = WannierizeWorkflow(self.settings, self.master_calcs)

            # Perform the wannierisation workflow within the init directory
            self.run_subworkflow(wannier_workflow, subdirectory='init')

            # Now, convert the files over from w90 format to (k)cp format
            fold_workflow = FoldToSupercellWorkflow(self.settings, self.master_calcs)

            # Do this in the same directory as the wannierisation
            self.run_subworkflow(fold_workflow, subdirectory='init/wannier')

            # We need a dummy calc before the real dft_init in order
            # to copy the previously calculated Wannier functions
            calc = self.new_calculator('kcp', 'dft_dummy')
            calc.directory = 'init'
            self.run_calculator(calc, enforce_ss=False)

            # DFT restarting from Wannier functions (after copying the Wannier functions)
            calc = self.new_calculator('kcp', 'dft_init', restart_mode='restart',
                                       restart_from_wannier_pwscf=True, do_outerloop=True, ndw=ndw_final)
            calc.directory = 'init'
            restart_dir = f'{calc.parameters.outdir}/{calc.parameters.prefix}_{calc.parameters.ndr}.save/K00001'
            for typ in ['occ', 'emp']:
                if typ == 'occ':
                    evcw_file = f'init/wannier/occ/evcw.dat'
                    dest_file = 'evc_occupied.dat'
                    if os.path.isfile(f'{evcw_file}'):
                        utils.system_call(f'cp {evcw_file} {restart_dir}/{dest_file}')
                    else:
                        raise OSError(f'Could not find {evcw_file}')
                if typ == 'emp':
                    for i_spin in ['1', '2']:
                        evcw_file = 'init/wannier/emp/evcw' + i_spin + '.dat'
                        dest_file = 'evc0_empty' + i_spin + '.dat'
                        if os.path.isfile(f'{evcw_file}'):
                            utils.system_call(f'cp {evcw_file} {restart_dir}/{dest_file}')
                        else:
                            raise OSError(f'Could not find {evcw_file}')

            self.run_calculator(calc, enforce_ss=False)

            # Check the consistency between the PW and CP band gaps
            pw_calc = [c for c in self.all_calcs if isinstance(
                c, calculators.PWCalculator) and c.parameters.calculation == 'nscf'][-1]
            pw_gap = pw_calc.results['lumo_ene'] - pw_calc.results['homo_ene']
            cp_gap = calc.results['lumo_energy'] - calc.results['homo_energy']
            if abs(pw_gap - cp_gap) > 2e-2 * pw_gap:
                raise ValueError(f'PW and CP band gaps are not consistent: {pw_gap} {cp_gap}')

            # The CP restarting from Wannier functions must be already converged
            Eini = calc.results['convergence']['filled'][0]['Etot']
            Efin = calc.results['energy']
            if abs(Efin - Eini) > 1e-6 * abs(Efin):
                raise ValueError(f'Too much difference between the initial and final CP energies: {Eini} {Efin}')

        elif self.init_orbitals == 'from old ki':
            self.print('Copying the density and orbitals from a pre-existing KI calculation')

            # Read the .cpi file to work out the value for ndw
            calc = calculators.KoopmansCPCalculator(qe_files='init/ki_init')

            # Move the old save directory to correspond to ndw_final, using pz_innerloop_init to work out where the
            # code will expect the tmp files to be
            old_savedir = f'{calc.parameters.outdir}/{calc.parameters.prefix}_{calc.parameters.ndw}.save'
            savedir = f'{self.new_calculator("kcp", "pz_innerloop_init").parameters.outdir}/{calc.parameters.prefix}_{ndw_final}.save'
            if not os.path.isdir(old_savedir):
                raise ValueError(f'{old_savedir} does not exist; a previous '
                                 'and complete KI calculation is required '
                                 'if init_orbitals="from old ki"')
            if os.path.isdir(savedir):
                utils.system_call(f'rm -r {savedir}')

            # Use the chdir construct in order to create the directory savedir if it does not already exist and is
            # nested
            with utils.chdir(savedir):
                utils.system_call(f'rsync -a {old_savedir}/ .')

            # Check that the files defining the variational orbitals exist
            savedir += '/K00001'
            files_to_check = ['init/ki_init.cpo', f'{savedir}/evc01.dat', f'{savedir}/evc02.dat']

            if calc.parameters.empty_states_nbnd is not None and calc.parameters.empty_states_nbnd > 0:
                files_to_check += [f'{savedir}/evc0_empty{i}.dat' for i in [1, 2]]

            for fname in files_to_check:
                if not os.path.isfile(fname):
                    raise ValueError(f'Could not find {fname}')

            if not calc.is_complete():
                raise ValueError('init/ki_init.cpo is incomplete so cannot be used '
                                 'to initialise the density and orbitals')

            self.all_calcs.append(calc)

        elif self.functional in ['ki', 'pkipz']:
            calc = self.new_calculator('kcp', 'dft_init')
            calc.directory = 'init'
            self.run_calculator(calc, enforce_ss=self.enforce_spin_symmetry)

            # Use the KS eigenfunctions as better guesses for the variational orbitals
            self._overwrite_canonical_with_variational_orbitals(calc)

            if self.init_orbitals == 'kohn-sham':
                self._copy_most_recent_calc_to_ndw(ndw_final)
            elif self.init_orbitals == 'pz':
                calc = self.new_calculator('kcp', 'pz_innerloop_init', alphas=self.bands.alphas, ndw=ndw_final)
                calc.directory = 'init'
                if self.all_calcs[-1].parameters.nelec == 2:
                    # If we only have two electrons, then the filled manifold is trivially invariant under unitary
                    # transformations. Furthermore, the PZ functional is invariant w.r.t. unitary rotations of the
                    # empty states. Thus in this instance we can skip the initialisation of the manifold entirely
                    self.print('Skipping the optimisation of the variational orbitals since they are invariant under '
                               'unitary transformations')
                    self._copy_most_recent_calc_to_ndw(ndw_final)
                else:
                    self.run_calculator(calc)
            else:
                raise ValueError('Should not arrive here')

        elif self.functional == 'kipz':
            # DFT from scratch
            calc = self.new_calculator('kcp', 'dft_init')
            calc.directory = 'init'
            self.run_calculator(calc, enforce_ss=self.enforce_spin_symmetry)

            if self.init_orbitals == 'kohn-sham':
                # Initialise the density with DFT and use the KS eigenfunctions as guesses for the variational orbitals
                self._overwrite_canonical_with_variational_orbitals(calc)
                self._copy_most_recent_calc_to_ndw(ndw_final)
            elif self.init_orbitals == 'pz':
                # PZ from DFT (generating PZ density and PZ orbitals)
                calc = self.new_calculator('kcp', 'pz_init', ndw=ndw_final)
                calc.directory = 'init'
                self.run_calculator(calc)
            else:
                raise ValueError('Should not arrive here')

        else:
            raise ValueError("Should not arrive here; there must be an inconsistency between the above code and \
                             workflow.valid_settings")

        return

    def _copy_most_recent_calc_to_ndw(self, ndw):
        calc = self.all_calcs[-1]
        if calc.parameters.ndw != ndw:
            assert calc.is_complete(), 'Cannot copy results of a previous calculation that is not itself complete'
            save_prefix = f'{calc.parameters.outdir}/{calc.parameters.prefix}'
            utils.system_call(f'cp -r {save_prefix}_{calc.parameters.ndw}.save {save_prefix}_{ndw}.save')

    def _overwrite_canonical_with_variational_orbitals(self, calc: calculators.KoopmansCPCalculator) -> None:
        self.print('Overwriting the variational orbitals with Kohn-Sham orbitals')
        savedir = f'{calc.parameters.outdir}/{calc.parameters.prefix}_{calc.parameters.ndw}.save/K00001'
        utils.system_call(f'cp {savedir}/evc1.dat {savedir}/evc01.dat')
        utils.system_call(f'cp {savedir}/evc2.dat {savedir}/evc02.dat')
        if calc.parameters.empty_states_nbnd is not None and calc.parameters.empty_states_nbnd > 0:
            utils.system_call(f'cp {savedir}/evc_empty1.dat {savedir}/evc0_empty1.dat')
            utils.system_call(f'cp {savedir}/evc_empty2.dat {savedir}/evc0_empty2.dat')

    def perform_alpha_calculations(self) -> None:
        # Set up directories
        if not os.path.isdir('calc_alpha'):
            utils.system_call('mkdir calc_alpha')

        converged = False
        i_sc = 0

        if not self.from_scratch and os.path.isfile('alphas.pkl'):
            # Reloading alphas and errors from file
            self.print('Reloading alpha values from file')
            self.bands.alphas = pd.read_pickle('alphas.pkl')
            self.bands.errors = pd.read_pickle('errors.pkl')

        alpha_indep_calcs = []

        while not converged and i_sc < self.n_max_sc_steps:
            i_sc += 1
            iteration_directory = 'calc_alpha'
            _, outdir = self.master_calcs['kcp'].outdir.rsplit('/', 1)
            outdir = os.getcwd() + f'/{iteration_directory}/{outdir}'

            if not os.path.isdir(outdir):
                utils.system_call(f'mkdir {outdir}')

            # Setting up directories
            if self.n_max_sc_steps > 1:
                self.print('SC iteration {}'.format(i_sc), style='subheading')
                iteration_directory += f'/iteration_{i_sc}'
                if not os.path.isdir(iteration_directory):
                    utils.system_call(f'mkdir {iteration_directory}')

            # Do a KI/KIPZ calculation with the updated alpha values
            calc = self.new_calculator('kcp', calc_presets=self.functional.replace('pkipz', 'ki'),
                                       alphas=self.bands.alphas)
            calc.directory = iteration_directory

            if i_sc == 1:
                if self.functional == 'kipz' and not self.periodic:
                    # For the first KIPZ trial calculation, do the innerloop
                    calc.parameters.do_innerloop = True
            else:
                # For later SC loops, read in the matching calculation from the
                # previous loop rather than the initialisation calculations
                calc.parameters.ndr = calc.parameters.ndw

            # Run the calculation and store the result. Note that we only need to continue
            # enforcing the spin symmetry if the density will change
            self.run_calculator(calc, enforce_ss=self.enforce_spin_symmetry and i_sc > 1)
            alpha_dep_calcs = [calc]

            # Update the bands' self-Hartree and energies (assuming spin-symmetry)
            self.bands.self_hartrees = calc.results['orbital_data']['self-Hartree'][0]

            # Group the bands
            self.bands.assign_groups(allow_reassignment=True)

            skipped_orbitals = []
            # Loop over removing/adding an electron from/to each orbital
            for band in self.bands:
                # Working out what to print for the orbital heading (grouping skipped bands together)
                if band in self.bands.to_solve or band == self.bands[-1]:
                    if band not in self.bands.to_solve:
                        skipped_orbitals.append(band.index)
                    if len(skipped_orbitals) > 0:
                        if len(skipped_orbitals) == 1:
                            self.print(f'Orbital {skipped_orbitals[0]}', style='subheading')
                        else:
                            orb_range = f'{skipped_orbitals[0]}-{skipped_orbitals[-1]}'
                            self.print(f'Orbitals {orb_range}', style='subheading')
                        self.print(f'Skipping; will use the screening parameter of an equivalent orbital')
                        skipped_orbitals = []
                    if band not in self.bands.to_solve:
                        continue
                else:
                    # Skip the bands which can copy the screening parameter from another
                    # calculation in the same orbital group
                    skipped_orbitals.append(band.index)
                    continue
                self.print(f'Orbital {band.index}', style='subheading')

                # Set up directories
                directory = f'{iteration_directory}/orbital_{band.index}'
                if not os.path.isdir(directory):
                    utils.system_call(f'mkdir {directory}')
                outdir_band = f'{outdir}/orbital_{band.index}'

                # Link tmp files from band-independent calculations
                if not os.path.isdir(outdir_band):
                    utils.system_call(f'mkdir {outdir_band}')
                    utils.system_call(
                        f'ln -sr {self.master_calcs["kcp"].parameters.outdir}/*.save {outdir_band}')

                # Don't repeat if this particular alpha_i was converged
                if i_sc > 1 and abs(band.error) < self.alpha_conv_thr:
                    self.print(f'Skipping band {band.index} since this alpha is already converged')
                    for b in self.bands:
                        if b == band or (band.group is not None and b.group == band.group):
                            b.alpha = band.alpha
                            b.error = band.error
                    continue

                # When we write/update the alpharef files in the work directory
                # make sure to include the fixed band alpha in file_alpharef.txt
                # rather than file_alpharef_empty.txt
                band_filled_or_fixed = [b is band or b.filled for b in self.bands]
                if band.filled:
                    index_empty_to_save = None
                else:
                    index_empty_to_save = band.index - self.bands.num(filled=True)

                # Perform the fixed-band-dependent calculations
                if self.functional in ['ki', 'pkipz']:
                    if band.filled:
                        calc_types = ['ki_frozen', 'dft_frozen', 'dft_n-1']
                    else:
                        calc_types = ['pz_print', 'dft_n+1_dummy', 'dft_n+1',
                                      'dft_n+1-1_frozen', 'ki_n+1-1_frozen']
                else:
                    if band.filled:
                        calc_types = ['kipz_frozen', 'dft_frozen', 'kipz_n-1']
                    else:
                        calc_types = ['kipz_print', 'dft_n+1_dummy', 'kipz_n+1',
                                      'dft_n+1-1_frozen', 'kipz_n+1-1_frozen']

                for calc_type in calc_types:
                    if self.functional in ['ki', 'pkipz']:
                        # Only the KI calculations change with alpha; we don't need to
                        # redo any of the others
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
                        alphas = [alphas + [alphas[-1]], alphas]
                        filling = [band_filled_or_fixed + [False], self.bands.filling]
                    else:
                        alphas = self.bands.alphas
                        filling = band_filled_or_fixed

                    # Set up calculator
                    calc = self.new_calculator('kcp', calc_type, alphas=alphas, filling=filling, fixed_band=min(
                        band.index, self.bands.num(filled=True) + 1),
                        index_empty_to_save=index_empty_to_save, outdir=outdir_band)
                    calc.directory = directory

                    # Run kcp.x
                    self.run_calculator(calc)

                    # Reset the value of 'fixed_band' so we can keep track of which calculation
                    # is which. This is important for empty orbital calculations, where fixed_band
                    # is always set to the LUMO but in reality we're fixing the band corresponding
                    # to index_empty_to_save from an earlier calculation
                    calc.parameters.fixed_band = band.index

                    # Store the result
                    # We store the results in one of two lists: alpha_indep_calcs and
                    # alpha_dep_calcs. The latter is overwritten at each new self-
                    # consistency loop.
                    if 'ki' in calc_type and 'print' not in calc_type:
                        alpha_dep_calcs.append(calc)
                    elif 'dft' in calc_type and 'dummy' not in calc_type:
                        if self.functional in ['ki', 'pkipz']:
                            # For KI, the DFT calculations are independent of alpha so
                            # we store these in a list that is never overwritten
                            alpha_indep_calcs.append(calc)
                        else:
                            # For KIPZ, the DFT calculations are dependent on alpha via
                            # the definition of the variational orbitals. We only want to
                            # store the calculations that used the most recent value of alpha
                            alpha_dep_calcs.append(calc)

                    # Copying of evcfixed_empty.dat to evc_occupied.dat
                    if calc_type in ['pz_print', 'kipz_print']:
                        evcempty_dir = f'{outdir_band}/{calc.parameters.prefix}_{calc.parameters.ndw}.save/K00001/'
                    elif calc_type == 'dft_n+1_dummy':
                        evcocc_dir = f'{outdir_band}/{calc.parameters.prefix}_{calc.parameters.ndr}.save/K00001/'
                        if os.path.isfile(f'{evcempty_dir}/evcfixed_empty.dat'):
                            utils.system_call(f'cp {evcempty_dir}/evcfixed_empty.dat {evcocc_dir}/evc_occupied.dat')
                        else:
                            raise OSError(f'Could not find {evcempty_dir}/evcfixed_empty.dat')

                if self.from_scratch:

                    # Calculate an updated alpha and a measure of the error
                    # E(N) - E_i(N - 1) - lambda^alpha_ii(1)     (filled)
                    # E_i(N + 1) - E(N) - lambda^alpha_ii(0)     (empty)
                    #
                    # Note that we only do this if the calculation has not been skipped
                    # because calculate_alpha reads in alpha values from files which get
                    # overwritten by subsequent calculations

                    calcs = [c for calc_set in [alpha_dep_calcs, alpha_indep_calcs]
                             for c in calc_set if c.parameters.fixed_band == band.index]

                    alpha, error = self.calculate_alpha_from_list_of_calcs(calcs, filled=band.filled)

                    for b in self.bands:
                        if b == band or (b.group is not None and b.group == band.group):
                            b.alpha = alpha
                            b.error = error

                    self.bands.alpha_history.to_pickle('alphas.pkl')
                    self.bands.error_history.to_pickle('errors.pkl')

            self.bands.print_history()

            converged = all([abs(b.error) < 1e-3 for b in self.bands])

        if converged:
            self.print('Screening parameters have been converged')
        else:
            self.print('Screening parameters have been determined but are not necessarily converged')

    def perform_final_calculations(self) -> None:

        directory = 'final'
        if not os.path.isdir(directory):
            utils.system_call(f'mkdir {directory}')

        if self.functional == 'pkipz':
            final_calc_types = ['ki', 'pkipz']
        else:
            final_calc_types = [self.functional]

        for final_calc_type in final_calc_types:

            final_calc_type += '_final'

            # For pKIPZ, the appropriate ndr can change but it is always ndw of the previous
            # KI calculation
            if final_calc_type == 'pkipz_final':
                ndr = [c.parameters.ndw for c in self.all_calcs if c.prefix in [
                    'ki', 'ki_final'] and hasattr(c.parameters, 'ndw')][-1]
                calc = self.new_calculator('kcp', final_calc_type, ndr=ndr, write_hr=True)
            else:
                calc = self.new_calculator('kcp', final_calc_type, write_hr=True)

            calc.directory = directory

            self.run_calculator(calc)

    def perform_postprocessing(self) -> None:
        # Import these here so that if these have been monkey-patched, we get the monkey-patched version
        from koopmans.workflows import WannierizeWorkflow

        if self.master_calcs['ui'].parameters.do_smooth_interpolation:
            factors = self.master_calcs['ui'].parameters.smooth_int_factor
            # Run the PW + W90 for a much larger grid
            local_master_calcs = copy.deepcopy(self.master_calcs)
            for name in ['pw', 'pw2wannier', 'w90_occ', 'w90_emp']:
                calc = local_master_calcs[name]
                original_k_grid = calc.parameters.get('kpts', None)
                if original_k_grid is not None:
                    # For calculations with a kpts attribute, create a denser k-grid
                    new_k_grid = [x * y for x, y in zip(original_k_grid, factors)]
                    calc.parameters['kpts'] = new_k_grid
            wannier_workflow = WannierizeWorkflow(self.settings, local_master_calcs)

            # Here, we allow for skipping of the smooth dft calcs (assuming they have been already run)
            # This is achieved via the optional argument of from_scratch in run_subworkflow(), which
            # overrides the value of wannier_workflow.from_scratch, as well as preventing the inheritance of
            # self.from_scratch to wannier_workflow.from_scratch and back again after the subworkflow finishes
            from_scratch = getattr(self, 'redo_preexisting_smooth_dft_calcs', None)
            self.run_subworkflow(wannier_workflow, subdirectory='postproc', from_scratch=from_scratch)

        for calc_presets in ['occ', 'emp']:
            calc = self.new_calculator('ui', calc_presets)
            self.run_calculator(calc, enforce_ss=False)

        # Merge the two calculations to print out the DOS and bands
        calc = self.new_calculator('ui', 'merge')
        calc.alat_sc = self.all_calcs[-1].alat_sc

        # Merge the bands
        energies = [c.results['band structure'].energies for c in self.all_calcs[-2:]]
        reference = np.max(energies[0])
        calc.results['band structure'] = BandStructure(
            self.all_calcs[-1].kpath, np.concatenate(energies, axis=2) - reference)

        if calc.do_dos:
            # Generate the DOS
            calc.calc_dos()

        # Store the calculator in the workflow's list of all the calculators
        self.all_calcs.append(calc)

        # Print out the merged bands and DOS
        if self.from_scratch:
            with utils.chdir('postproc'):
                calc.write_results()

    def new_calculator(self, calc_type: str, calc_presets: str = 'dft_init', alphas: Optional[Union[List[float], List[List[float]]]] = None, filling: Optional[Union[List[List[bool]], List[bool]]] = None, **kwargs) -> calculators.ExtendedCalculator:
        """

        Generates a new calculator based on the self.master_calc[calc_type]
        template calculator

        """

        if calc_type == 'kcp':
            return self.new_kcp_calculator(calc_presets, alphas, filling, **kwargs)
        elif calc_type == 'ui':
            return self.new_ui_calculator(calc_presets, **kwargs)
        else:
            raise ValueError(f'Invalid calc type {calc_type}; must be "kcp"/"ui"')

    def new_kcp_calculator(self, calc_presets: str = 'dft_init', alphas: Optional[Union[List[float], List[List[float]]]] = None, filling: Optional[Union[List[List[bool]], List[bool]]] = None, **kwargs) -> calculators.KoopmansCPCalculator:
        """

        Generates a new KCP calculator based on the self.master_calc["kcp"]
        template calculator, modifying the appropriate settings to match the
        chosen calc_presets, and altering any Quantum Espresso keywords
        specified as kwargs

        Arguments:

            calc_presets
                The set of preset values to use; must be one of the following strings:

                Initialisation
                'dft_init'            DFT calculation from scratch
                'pz_init'             PZ calculation starting from DFT restart
                'pz_innerloop_init'   PZ calculation starting from DFT restart (innerloop only)
                'dft_dummy'           DFT dummy calculation that generate store files
                                      for a periodic calculation restarting from Wannier functions

                Trial calculations
                'ki'     KI calculation with N electrons and empty bands if specified
                'kipz'   As above, but for KIPZ

                For calculating alpha_i for filled orbitals.
                'dft_frozen'    DFT calculation leaving rho unchanged, for reporting energy and lambda
                'ki_frozen'     KI calculation with N electrons for generating lambda only (will not alter density)
                'kipz_frozen'   KIPZ calculation with N electrons for generating lambda only (will not alter density)
                'dft_n-1'       DFT calculation with N-1 electrons via fixed_state; rho is optimised
                'kipz_n-1'      KIPZ calculation with N-1 electrons via fixed_state; rho is optimised

                For calculating alpha_i for empty orbitals
                'pz_print'         PZ calculation that generates evcfixed_empty.dat file
                'kipz_print'       KIPZ calculation that generates evcfixed_empty.dat file
                'dft_n+1_dummy'    DFT dummy calculation that generates store files of
                                   the correct dimensions
                'dft_n+1'          DFT calculation with N+1 electrons; rho is optimised
                'kipz_n+1'         KIPZ calculation with N+1 electrons; rho is optimised
                'dft_n+1-1_frozen' DFT calculation with N electrons, starting from N+1
                                   but with f_cutoff = 0.00001
                'ki_n+1-1_frozen'  KI calculation with N electrons, starting from N+1
                                   but with f_cutoff = 0.00001
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
        calc = calculators.KoopmansCPCalculator(self.master_calcs["kcp"], alphas=alphas, filling=filling)

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
        elif calc_presets in ['ki_frozen', 'kipz_frozen']:
            ndr = 60
            ndw = 61
        elif calc_presets == 'dft_frozen':
            ndr = 60
            ndw = 62
        elif calc_presets in ['dft_n-1', 'kipz_n-1']:
            ndr = 60
            ndw = 63
        elif calc_presets in ['pz_print', 'kipz_print']:
            ndr = 60
            ndw = 64
        elif calc_presets == 'dft_n+1_dummy':
            ndr = 65
            ndw = 65
        elif calc_presets in ['ki_n+1-1_frozen', 'kipz_n+1-1_frozen']:
            ndr = 65
            ndw = 66
        elif calc_presets == 'dft_n+1-1_frozen':
            ndr = 65
            ndw = 67
        elif calc_presets in ['dft_n+1', 'kipz_n+1']:
            ndr = 65
            ndw = 68
        elif calc_presets in ['ki_final', 'kipz_final']:
            if self.calculate_alpha:
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
            calc.parameters.nelup += 1
            if 'dummy' not in calc.prefix:
                calc.parameters.restart_from_wannier_pwscf = True

        # electrons
        # For all calculations calculating alpha, remove the empty states and
        # increase the energy thresholds
        if not any([s in calc.prefix for s in ['init', 'print', 'final']]) and \
           calc.prefix not in ['ki', 'kipz', 'dft_dummy']:
            calc.parameters.empty_states_nbnd = 0
            calc.parameters.conv_thr *= 100
            calc.parameters.esic_conv_thr *= 100

        if self.periodic and not any([s == calc.prefix for s in ['dft_init', 'dft_n-1', 'dft_n+1',
                                                                 'kipz', 'kipz_n-1', 'kipz_n+1']]):
            calc.parameters.do_outerloop = False
            calc.parameters.do_innerloop = False
        elif any([s in calc.prefix for s in ['frozen', 'dummy', 'print', 'innerloop']]) or calc.prefix == 'pkipz_final':
            calc.parameters.do_outerloop = False
            if calc.parameters.empty_states_nbnd > 0:
                calc.parameters.do_outerloop_empty = False
        elif calc.prefix in ['ki', 'ki_final']:
            calc.parameters.do_outerloop = False
            if calc.parameters.empty_states_nbnd > 0:
                if self.init_empty_orbitals == 'pz':
                    calc.parameters.do_outerloop_empty = True
                else:
                    calc.parameters.do_outerloop_empty = False
        else:
            calc.parameters.do_outerloop = True
            if calc.parameters.empty_states_nbnd > 0:
                calc.parameters.do_outerloop_empty = True

        if calc.parameters.maxiter is None and calc.parameters.do_outerloop:
            calc.parameters.maxiter = 300
        if calc.parameters.empty_states_maxstep is None and calc.parameters.do_outerloop_empty:
            calc.parameters.empty_states_maxstep = 300

        # No empty states minimization in the solids workflow for the moment
        if self.periodic and calc.parameters.empty_states_nbnd > 0:
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
        if calc.parameters.empty_states_nbnd > 0:
            calc.parameters.do_innerloop_empty = False
        if 'kipz' in calc.prefix:
            calc.parameters.which_orbdep = 'nkipz'
        elif 'pz' in calc.prefix:
            calc.parameters.which_orbdep = 'pz'
        elif 'ki' in calc.prefix:
            calc.parameters.which_orbdep = 'nki'
        if 'print' in calc.prefix:
            calc.parameters.print_wfc_anion = True

        if self.mt_correction:
            calc.parameters.which_compensation = 'tcc'
        else:
            calc.parameters.which_compensation = 'none'

        # If we are using frozen orbitals, we override the above logic and freeze the variational orbitals
        # post-initialisation
        if self.frozen_orbitals and 'init' not in calc.prefix and not any([s == calc.prefix for s in
                                                                           ['dft_n-1', 'dft_n+1', 'kipz_n-1',
                                                                            'kipz_n+1']]):
            calc.parameters.do_outerloop = False
            calc.parameters.do_innerloop = False
            if calc.parameters.empty_states_nbnd > 0:
                calc.parameters.do_outerloop_empty = False
                calc.parameters.do_innerloop_empty = False

        # Handle any keywords provided by kwargs
        # Note that since this is performed after the above logic this can (deliberately
        # or accidentally) overwrite the above settings

        for keyword, value in kwargs.items():
            setattr(calc.parameters, keyword, value)

        # Sanity checking
        if calc.parameters.print_wfc_anion and calc.parameters.index_empty_to_save is None:
            raise ValueError('Error: print_wfc_anion is set to true but you have not selected '
                             'an index_empty_to_save. Provide this as an argument to new_cp_calculator')

        if calc.parameters.fixed_band is not None and calc.parameters.fixed_band > calc.parameters.nelup + 1:
            utils.warn('calc.fixed_band is higher than the LUMO; this should not happen')

        # avoid innerloops for one-orbital-manifolds
        if calc.parameters.nelup in [0, 1] and calc.parameters.neldw in [0, 1]:
            calc.parameters.do_innerloop = False
        if calc.parameters.empty_states_nbnd == 1:
            calc.parameters.do_innerloop_empty = False

        # don't print QC in some cases
        if 'dummy' in calc.prefix:
            calc.skip_qc = True
        elif calc.prefix[-2:] == '+1':
            # Don't check N+1 energies because they're known to be unreliable
            if 'energy' in calc.results_for_qc:
                calc.results_for_qc.remove('energy')

        return calc

    def new_ui_calculator(self, calc_presets: str, **kwargs) -> calculators.UnfoldAndInterpolateCalculator:

        valid_calc_presets = ['occ', 'emp', 'merge']
        assert calc_presets in valid_calc_presets, 'In KoopmansDSCFWorkflow.new_ui_calculator() calc_presets must be ' \
            '/'.join([f'"{s}"' for s in valid_calc_presets]) + ', but you have tried to set it equal to {calc_presets}'

        if calc_presets == 'merge':
            # Dummy calculator for merging bands and dos
            calc = calculators.UnfoldAndInterpolateCalculator(calc=self.master_calcs['ui'])
            pass
        else:
            calc = calculators.UnfoldAndInterpolateCalculator(calc=self.master_calcs[f'ui_{calc_presets}'])
            # Automatically generating UI calculator settings
            calc.directory = f'postproc/{calc_presets}'
            calc.parameters.kc_ham_file = os.path.abspath(f'final/ham_{calc_presets}_1.dat')
            calc.parameters.sc_dim = self.master_calcs['pw'].parameters.kpts
            calc.parameters.w90_seedname = os.path.abspath(f'init/wannier/{calc_presets}/wann')
            # The supercell can be obtained from the most recent CP calculation
            cell = [c.atoms.get_cell() for c in self.all_calcs if isinstance(c, calculators.KoopmansCPCalculator)][-1]
            calc.parameters.alat_sc = np.linalg.norm(cell[0])
            if calc.parameters.do_smooth_interpolation:
                calc.parameters.dft_smooth_ham_file = os.path.abspath(f'postproc/wannier/{calc_presets}/wann_hr.dat')
                calc.parameters.dft_ham_file = os.path.abspath(f'init/wannier/{calc_presets}/wann_hr.dat')
        calc.prefix = self.functional

        return calc

    def calculate_alpha_from_list_of_calcs(self, calcs: List[calculators.KoopmansCPCalculator], filled: bool = True) -> Union[Tuple[float, float]]:
        '''

        Calculates alpha via equation 10 of Nguyen et. al (2018) 10.1103/PhysRevX.8.021051
        If the band is filled, use s = 1; if the band is empty, use s = 0

        Arguments:
            calcs          -- a list of selected calculations from which to calculate alpha
            filled         -- True if the orbital for which we're calculating alpha is filled

        '''

        if self.functional == 'kipz':
            if filled:
                # KIPZ N
                [alpha_calc] = [c for c in calcs if c.parameters.which_orbdep == 'nkipz'
                                and c.parameters.do_orbdep and c.parameters.f_cutoff == 1.0]
                kipz = alpha_calc.results

                # KIPZ N-1
                [kipz_m1_calc] = [c for c in calcs if c.parameters.which_orbdep == 'nkipz'
                                  and c.parameters.do_orbdep and c.parameters.f_cutoff < 0.0001]
                kipz_m1 = kipz_m1_calc.results
                charge = 1 - kipz_m1_calc.f_cutoff

                # DFT N
                [dft] = [c.results for c in calcs if not c.parameters.do_orbdep
                         and c.parameters.restart_mode == 'restart' and c.parameters.f_cutoff == 1.0]

                dE = kipz['energy'] - kipz_m1['energy']
                lambda_a = kipz['lambda_ii']
                lambda_0 = dft['lambda_ii']
                mp1 = kipz_m1['mp1_energy']
                mp2 = kipz_m1['mp2_energy']

            else:
                # KIPZ N+1-1
                [alpha_calc] = [c for c in calcs if c.parameters.which_orbdep == 'nkipz'
                                and c.parameters.do_orbdep and c.parameters.f_cutoff < 0.0001]
                kipz = alpha_calc.results

                # KIPZ N+1
                [kipz_p1_calc] = [c for c in calcs if c.parameters.which_orbdep == 'nkipz'
                                  and c.parameters.do_orbdep and c.parameters.f_cutoff == 1.0]
                kipz_p1 = kipz_p1_calc.results
                charge = - kipz_p1_calc.parameters.f_cutoff

                # DFT N+1-1
                [dft] = [c.results for c in calcs if not c.parameters.do_orbdep
                         and c.parameters.restart_mode == 'restart' and c.parameters.f_cutoff < 0.0001]

                dE = kipz_p1['energy'] - kipz['energy']
                lambda_a = kipz['lambda_ii']
                lambda_0 = dft['lambda_ii']
                mp1 = kipz_p1['mp1_energy']
                mp2 = kipz_p1['mp2_energy']

        else:
            # self.functional in ['ki', 'pkipz']
            if filled:
                # KI N
                [alpha_calc] = [c for c in calcs if c.parameters.which_orbdep == 'nki'
                                and c.parameters.do_orbdep and c.parameters.f_cutoff == 1.0]
                ki = alpha_calc.results

                # DFT N-1
                [dft_m1_calc] = [c for c in calcs if not c.parameters.do_orbdep
                                 and c.parameters.restart_mode == 'restart' and c.parameters.f_cutoff < 0.0001]
                dft_m1 = dft_m1_calc.results
                charge = 1 - dft_m1_calc.parameters.f_cutoff

                # DFT N
                [dft] = [c.results for c in calcs if not c.parameters.do_orbdep
                         and c.parameters.restart_mode == 'restart' and c.parameters.f_cutoff == 1.0]

                dE = dft['energy'] - dft_m1['energy']
                lambda_a = ki['lambda_ii']
                lambda_0 = dft['lambda_ii']
                mp1 = dft_m1['mp1_energy']
                mp2 = dft_m1['mp2_energy']

            else:
                # KI N+1-1
                [alpha_calc] = [c for c in calcs if c.parameters.which_orbdep == 'nki'
                                and c.parameters.do_orbdep and c.parameters.f_cutoff < 0.0001]
                ki = alpha_calc.results

                # DFT N+1
                [dft_p1_calc] = [c for c in calcs if not c.parameters.do_orbdep
                                 and c.parameters.restart_mode == 'restart' and c.parameters.f_cutoff == 1.0]
                dft_p1 = dft_p1_calc.results
                charge = - dft_p1_calc.parameters.f_cutoff

                # DFT N+1-1
                [dft] = [c.results for c in calcs if not c.parameters.do_orbdep
                         and c.parameters.restart_mode == 'restart' and c.parameters.f_cutoff < 0.0001]

                dE = dft_p1['energy'] - dft['energy']
                lambda_a = ki['lambda_ii']
                lambda_0 = dft['lambda_ii']
                mp1 = dft_p1['mp1_energy']
                mp2 = dft_p1['mp2_energy']

            # Check total energy is unchanged
            dE_check = abs(ki['energy'] - dft['energy'])
            if dE_check > 1e-5:
                utils.warn('KI and DFT energies differ by {:.5f} eV'.format(dE_check))

        # Obtaining alpha
        if (alpha_calc.parameters.odd_nkscalfact and filled) or (alpha_calc.parameters.odd_nkscalfact_empty and not filled):
            alpha_guesses = alpha_calc.alphas[0]
            alpha_guess = alpha_guesses[alpha_calc.parameters.fixed_band - 1]
        else:
            alpha_guess = alpha_calc.parameters.nkscalfact

        # Checking Makov-Payne correction energies and applying them (if needed)
        if self.mp_correction:
            if mp1 is None:
                raise ValueError('Could not find 1st order Makov-Payne energy')
            if mp2 is None:
                utils.warn('Could not find 2nd order Makov-Payne energy; applying first order only')
                mp_energy = mp1
            else:
                mp_energy = mp1 + mp2

            dE -= np.sign(charge) * mp_energy / self.eps_inf

        alpha = alpha_guess * (dE - lambda_0) / (lambda_a - lambda_0)

        # The error is lambda^alpha(1) - lambda^alpha_i(1)
        error = dE - lambda_a

        return alpha, error
