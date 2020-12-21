"""

Workflow module for python_KI, containing the workflow for performing KI and KIPZ calculations

Written by Edward Linscott Jan 2020
Split off from workflow.py Oct 2020

"""

import os
import subprocess
import copy
import numpy as np
import pandas as pd
from koopmans import io, utils
from koopmans.calculators.cp import CP_calc
from koopmans.workflows.generic import Workflow
from koopmans.workflows.wf_with_w90 import WannierizeWorkflow


class KoopmansWorkflow(Workflow):

    def __init__(self, workflow_settings, calcs_dct, alphas=None):
        super().__init__(workflow_settings, calcs_dct)

        if 'cp' not in self.master_calcs:
            raise ValueError(
                'Performing a KC calculation requires a "cp" block in the .json input file')
        cp_calc = self.master_calcs['cp']

        # Sanitise outdir
        if cp_calc.outdir[0] != '/':
            cp_calc.outdir = os.getcwd() + '/' + cp_calc.outdir.strip('./')

        # If periodic, convert the cp calculation into a Γ-only supercell calculation
        if self.periodic:
            kpts = self.master_calcs['cp'].calc.parameters.kpts
            self.master_calcs['cp'].transform_to_supercell(np.diag(kpts))
            self.orbital_groups = [i for i in self.orbital_groups for _ in range(np.prod(kpts))]

        # Check the number of empty states has been correctly configured
        if self.init_variational_orbitals in ['mlwfs', 'projw']:
            w90_emp_calc = self.master_calcs['w90_emp']

            expected_empty_states_nbnd = w90_emp_calc.num_wann * np.prod(cp_calc.calc.parameters.kpts)
            if cp_calc.empty_states_nbnd is None:
                cp_calc.empty_states_nbnd = expected_empty_states_nbnd
            elif cp_calc.empty_states_nbnd != expected_empty_states_nbnd:
                raise ValueError('cp empty_states_nbnd and wannier90 num_wann (emp) are inconsistent')

        # Preparing panda dataframes in which to store results
        i_bands = range(1, len(cp_calc.filling[0]) + 1)
        self.alpha_df = pd.DataFrame(columns=i_bands)
        self.error_df = pd.DataFrame(columns=i_bands)

        if alphas is not None:
            print('Using alpha values provided from a previous calculation')
            self.alpha_df.loc[1] = alphas
        elif self.alpha_from_file:
            print('Reading alpha values from file')
            self.read_alphadf_from_file()
        else:
            print('Initialising alpha with a guess')
            self.alpha_df.loc[1] = [self.alpha_guess for _ in i_bands]

    def read_alphadf_from_file(self, directory='.'):
        '''
        This routine reads in the contents of file_alpharef.txt and file_alpharef_empty.txt and
        stores the result in self.alpha_df

        Note that io.read_alpha_file provides a flattened list of alphas so we must exclude
        duplicates if calc.nspin == 2
        '''

        calc = self.master_calcs['cp']

        if calc.nspin == 2:
            i_alphas = list(range(0, calc.nelec // 2)) + list((range(calc.nelec, calc.nelec + calc.empty_states_nbnd)))
        else:
            i_alphas = list((range(0, calc.nelec // 2 + calc.empty_states_nbnd)))

        alphas = io.read_alpha_file(directory)
        alphas = [a for i, a in enumerate(alphas) if i in i_alphas]
        alphas = pd.DataFrame([alphas], columns=self.alpha_df.columns)
        self.alpha_df = self.alpha_df.append(alphas)

    def run(self):
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
            if self.init_density != 'ki':
                # if self.init_density == "ki" we don't want to delete the directory containing
                # the KI calculation we're reading the manifold from, or the TMP files
                utils.system_call('rm -r init 2>/dev/null', False)
                utils.system_call(f'rm -r {self.master_calcs["cp"].outdir} 2>/dev/null', False)
            utils.system_call('rm -r calc_alpha 2>/dev/null', False)
            utils.system_call('rm -r final 2>/dev/null', False)

        print('\nINITIALISATION OF DENSITY')
        self.perform_density_initialisation()

        print('\nINITIALISATION OF MANIFOLD')
        self.perform_manifold_initialisation()

        if self.from_scratch and self.init_variational_orbitals != 'ki' and self.enforce_spin_symmetry \
                and self.init_variational_orbitals not in ['mlwfs', 'projw']:
            print('Copying the spin-up variational orbitals over to the spin-down channel')
            calc = self.all_calcs[-1]
            savedir = f'{calc.outdir}/{calc.prefix}_{calc.ndw}.save/K00001'
            utils.system_call(f'cp {savedir}/evc01.dat {savedir}/evc02.dat')
            if calc.empty_states_nbnd is not None and calc.empty_states_nbnd > 0:
                utils.system_call(
                    f'cp {savedir}/evc0_empty1.dat {savedir}/evc0_empty2.dat')

        print('\nDETERMINING ALPHA VALUES')
        if self.calculate_alpha:
            self.perform_alpha_calculations()
        else:
            print('Skipping calculation of screening parameters', end='')
            if self.alpha_df.empty:
                print('; reading values from file')
                self.read_alphadf_from_file()
            else:
                print('')
            self.print_summary()

        # Final calculation
        print(f'\nFINAL {self.functional.upper().replace("PK","pK")} CALCULATION')
        self.perform_final_calculations()

        # Print out additional quality control data for final calculation
        if self.print_qc:
            calc = self.all_calcs[-1]
            for i, alpha in enumerate(self.alpha_df.iloc[-1]):
                if not self.error_df.empty and not pd.isnull(self.error_df.iloc[-1, i]):
                    self.print_qc_keyval(f'alpha({i+1})', alpha)
            for isp, orbs_self_hartree in enumerate(calc.results['orbital_data']['self-Hartree']):
                for i, orb_sh in enumerate(orbs_self_hartree):
                    self.print_qc_keyval(f'orbital_self_Hartree(orb={i+1},sigma={isp+1})', orb_sh)

        print('\nWORKFLOW COMPLETE')

    def perform_density_initialisation(self):

        if self.init_variational_orbitals in ['mlwfs', 'projw']:
            # Wannier functions using pw.x, wannier90.x and pw2wannier90.x
            if self.init_density != 'pbe':
                raise ValueError('wannierize requires init_density to be pbe')
            wannier_workflow = WannierizeWorkflow(self.settings, self.master_calcs)

            # Perform the wannierisation workflow within the init directory
            if not os.path.isdir('init'):
                utils.system_call('mkdir init')
            os.chdir('init')
            wannier_workflow.run()
            os.chdir('..')

        if self.init_density == 'pbe':
            if self.init_variational_orbitals in ['mlwfs', 'projw']:
                # We need a dummy calc before the real pbe_init in order
                # to copy the previously calculated Wannier functions
                calc = self.new_calculator('cp', 'pbe_init', do_outerloop=False,
                                           do_outerloop_empty=False, do_innerloop=False,
                                           do_innerloop_empty=False)
                calc.directory = 'init'
                calc.name = 'pbe_dummy'
                calc.skip_qc = True
                self.run_calculator(calc, enforce_ss=False)

                # PBE restarting from Wannier functions (after copying the Wannier functions)
                calc = self.new_calculator('cp', 'pbe_init', restart_mode='restart',
                                           restart_from_wannier_pwscf=True, do_outerloop_empty=False)
                calc.directory = 'init'
                restart_dir = f'{calc.outdir}/{calc.prefix}_{calc.ndr}.save/K00001'
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

            else:
                calc = self.new_calculator('cp', 'pbe_init')
                calc.directory = 'init'
                self.run_calculator(calc, enforce_ss=self.enforce_spin_symmetry)

        elif self.init_density == 'pbe-pw':
            # PBE using pw.x
            raise ValueError('init_density: pbe-pw is not yet implemented')

        elif self.init_density == 'pz':
            # PBE from scratch
            calc = self.new_calculator('cp', 'pbe_init')
            calc.directory = 'init'
            self.run_calculator(calc)
            # PZ from PBE
            calc = self.new_calculator('cp', 'pz_init')
            calc.directory = 'init'
            self.run_calculator(calc)

        elif self.init_density == 'ki':
            print('Copying the density from a pre-existing KI calculation')

            # Read the .cpi file to work out the value for ndw
            calc = CP_calc(qe_files='init/ki_init')

            # Move the old save directory to correspond to ndw = 50
            old_savedir = f'{calc.outdir}/{calc.prefix}_{calc.ndw}.save'
            savedir = f'{calc.outdir}/{calc.prefix}_50.save'
            if not os.path.isdir(old_savedir):
                raise ValueError(f'{old_savedir} does not exist; a previous '
                                 'and complete KI calculation is required '
                                 'if init_density=ki')
            if os.path.isdir(savedir):
                utils.system_call(f'rm -r {savedir}')
            utils.system_call(f'rsync -a {old_savedir}/ {savedir}/')

            # Check that the files defining the variational orbitals exist
            savedir += '/K00001'
            files_to_check = ['init/ki_init.cpo',
                              f'{savedir}/evc01.dat', f'{savedir}/evc02.dat']

            if calc.empty_states_nbnd is not None and calc.empty_states_nbnd > 0:
                files_to_check += [f'{savedir}/evc0_empty{i}.dat' for i in [1, 2]]

            for fname in files_to_check:
                if not os.path.isfile(fname):
                    raise ValueError(f'Could not find {fname}')

            if not calc.is_complete():
                raise ValueError('init/ki_init.cpo is incomplete so cannot be used '
                                 'to initialise the density')

            self.all_calcs.append(calc)
        else:
            raise ValueError(
                "Should not arrive here; compare the above code with workflow.valid_settings")

    def perform_manifold_initialisation(self):
        calc = self.all_calcs[-1]
        if self.init_density == 'pbe' and self.init_variational_orbitals not in ['mlwfs', 'projw']:
            # Using KS eigenfunctions as guesses for the variational orbitals
            print('Overwriting the CP variational orbitals with Kohn-Sham orbitals')
            savedir = f'{calc.outdir}/{calc.prefix}_{calc.ndw}.save/K00001'
            utils.system_call(f'cp {savedir}/evc1.dat {savedir}/evc01.dat')
            utils.system_call(f'cp {savedir}/evc2.dat {savedir}/evc02.dat')
            if calc.empty_states_nbnd is not None and calc.empty_states_nbnd > 0:
                utils.system_call(
                    f'cp {savedir}/evc_empty1.dat {savedir}/evc0_empty1.dat')
                utils.system_call(
                    f'cp {savedir}/evc_empty2.dat {savedir}/evc0_empty2.dat')

        if self.init_variational_orbitals == 'pz':
            calc = self.new_calculator('cp', 'pz_innerloop_init', alphas=self.alpha_df.loc[1])
        elif self.init_variational_orbitals == 'ki':
            if self.init_density != 'ki':
                raise ValueError('Initialising variational orbitals with KI makes no sense unless '
                                 'reading from a pre-existing KI calculation')
            print('Copying the density from a pre-existing KI calculation')

        calc.directory = 'init'

        if self.init_variational_orbitals == 'ki':
            pass
        elif self.all_calcs[-1].nelec == 2 and (self.init_variational_orbitals == 'pz' or calc.empty_states_nbnd == 0):
            # If we only have two electrons, then the filled manifold is trivially invariant under unitary
            # transformations. Likewise, if we have no empty states or if we're using a functional which is
            # invariant w.r.t. unitary rotations of the empty states, then the empty manifold need not be minimised
            # In these instances, we can skip the initialisation of the manifold entirely
            print('Skipping the optimisation of the variational orbitals since they are invariant under unitary '
                  'transformations')
            save_prefix = f'{calc.outdir}/{calc.prefix}'
            utils.system_call(
                f'cp -r {save_prefix}_{calc.ndr}.save {save_prefix}_{calc.ndw}.save')
        elif self.init_variational_orbitals in ['mlwfs', 'projw']:
            print('The variational orbitals have already been initialised to Wannier functions during the density '
                  'initialisation')
        elif self.init_variational_orbitals in [self.init_density, 'skip']:
            if self.init_variational_orbitals == 'skip':
                print('Skipping the optimisation of the variational orbitals')
            else:
                print('Skipping the optimisation of the variational orbitals since they were already optimised during '
                      'the density initialisation')
            save_prefix = f'{calc.outdir}/{calc.prefix}'
            utils.system_call(
                f'cp -r {save_prefix}_{calc.ndr}.save {save_prefix}_{calc.ndw}.save')
        else:
            self.run_calculator(calc)

    def perform_alpha_calculations(self):
        # Group the bands by group, and work out which bands to solve explicitly
        band_filling = self.all_calcs[-1].filling[0]
        n_filled_bands = band_filling.count(True)
        n_bands = len(band_filling)
        i_bands = range(1, n_bands + 1)

        if self.orbital_groups is None:
            self.orbital_groups = range(len(band_filling))

        bands_to_solve = {}
        groups_found = set([])

        for i_band in list(range(1, n_filled_bands + 1)[::-1]) \
                + list(range(n_filled_bands + 1, n_bands + 1)):
            # Looping through the filled bands from highest to lowest, then empty bands from
            # lowest to highest
            group = self.orbital_groups[i_band - 1]
            if group not in groups_found:
                groups_found.add(group)
                bands_to_solve[i_band] = [i + 1 for i,
                                          g in enumerate(self.orbital_groups) if g == group]
        if groups_found != set(self.orbital_groups):
            raise ValueError('Splitting of orbitals into groups failed')

        # Set up directories
        if not os.path.isdir('calc_alpha'):
            utils.system_call('mkdir calc_alpha')

        converged = False
        i_sc = 0

        if not self.from_scratch and os.path.isfile('alphas.pkl'):
            # Reloading alphas and errors from file
            print('Reloading alpha values from file')
            self.alpha_df = pd.read_pickle('alphas.pkl')
            self.error_df = pd.read_pickle('errors.pkl')

        alpha_indep_calcs = []

        while not converged and i_sc < self.n_max_sc_steps:
            i_sc += 1
            iteration_directory = 'calc_alpha'
            _, outdir = self.master_calcs['cp'].outdir.rsplit('/', 1)
            outdir = os.getcwd() + f'/{iteration_directory}/{outdir}'

            if not os.path.isdir(outdir):
                utils.system_call(f'mkdir {outdir}')

            # Setting up directories
            if self.n_max_sc_steps > 1:
                print('\n== SC iteration {} ==================='.format(i_sc))
                iteration_directory += f'/iteration_{i_sc}'
                if not os.path.isdir(iteration_directory):
                    utils.system_call(f'mkdir {iteration_directory}')

            # Do a KI/KIPZ calculation with the updated alpha values
            calc = self.new_calculator('cp', calc_presets=self.functional.replace('pkipz', 'ki'),
                                       alphas=self.alpha_df.loc[i_sc])
            calc.directory = iteration_directory

            if self.init_variational_orbitals in ['mlwfs', 'projw']:
                calc.ndr = 50
                calc.do_outerloop = False
                calc.do_outerloop_empty = False
                calc.restart_from_wannier_pwscf = True

            if i_sc == 1:
                if self.functional == 'kipz':
                    # For the first KIPZ trial calculation, do the innerloop
                    calc.do_innerloop = True
            else:
                # For later SC loops, read in the matching calculation from the
                # previous loop rather than the initialisation calculations
                calc.ndr = calc.ndw

            # Run the calculation and store the result. Note that we only need to continue
            # enforcing the spin symmetry if the density will change
            self.run_calculator(calc, enforce_ss=self.enforce_spin_symmetry and i_sc > 1)
            alpha_dep_calcs = [calc]

            # Loop over removing/adding an electron from/to each orbital
            for fixed_band, filled in zip(i_bands, band_filling):
                print('-- Orbital {} ------------------------'.format(fixed_band))
                # Skip the bands which can copy the screening parameter from another
                # calculation in the same orbital group
                if fixed_band not in bands_to_solve:
                    print(
                        f'Skipping; will use the screening parameter of an equivalent orbital')
                    continue
                all_bands_in_group = bands_to_solve[fixed_band]

                # Set up directories
                directory = f'{iteration_directory}/orbital_{fixed_band}'
                if not os.path.isdir(directory):
                    utils.system_call(f'mkdir {directory}')
                outdir_band = f'{outdir}/orbital_{fixed_band}'

                # Link tmp files from band-independent calculations
                if not os.path.isdir(outdir_band):
                    utils.system_call(f'mkdir {outdir_band}')
                    utils.system_call(
                        f'ln -sr {self.master_calcs["cp"].outdir}/*.save {outdir_band}')

                # Don't repeat if this particular alpha_i was converged
                if i_sc > 1 and any([abs(e) < self.alpha_conv_thr for e in
                                     self.error_df.loc[:i_sc - 1, fixed_band]]):
                    print(
                        f'Skipping band {fixed_band} since this alpha is already '
                        'converged')
                    self.alpha_df.loc[i_sc + 1,
                                      all_bands_in_group] = self.alpha_df.loc[i_sc, fixed_band]
                    self.error_df.lor[i_sc,
                                      all_bands_in_group] = self.error_df.loc[i_sc - 1, fixed_band]
                    continue

                # When we write/update the alpharef files in the work directory
                # make sure to include the fixed band alpha in file_alpharef.txt
                # rather than file_alpharef_empty.txt
                band_filled_or_fixed = [
                    b or i == fixed_band - 1 for i, b in enumerate(band_filling)]
                if filled:
                    index_empty_to_save = None
                else:
                    index_empty_to_save = fixed_band - n_filled_bands

                # Perform the fixed-band-dependent calculations
                if self.functional in ['ki', 'pkipz']:
                    if filled:
                        calc_types = ['ki_frozen', 'pbe_frozen', 'pbe_n-1']
                    else:
                        calc_types = ['pz_print', 'pbe_n+1_dummy', 'pbe_n+1',
                                      'pbe_n+1-1_frozen', 'ki_n+1-1_frozen']
                else:
                    if filled:
                        calc_types = ['kipz_frozen', 'pbe_frozen', 'kipz_n-1']
                    else:
                        calc_types = ['kipz_print', 'pbe_n+1_dummy', 'kipz_n+1',
                                      'pbe_n+1-1_frozen', 'kipz_n+1-1_frozen']

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
                        if i_sc > 1 and calc_type == 'pbe_n+1_dummy':
                            continue

                    if 'print' in calc_type:
                        # Note that the 'print' calculations for empty bands do not
                        # in fact involve the fixing of that band (and thus for the
                        # 'fixed' band the corresponding alpha should be in
                        # file_alpharef_empty.txt)
                        alphas = self.alpha_df.loc[i_sc]
                        filling = band_filling
                    elif not filled:
                        # In the case of empty orbitals, we gain an extra orbital in
                        # the spin-up channel, so we explicitly construct both spin
                        # channels for "alphas" and "filling"
                        alphas = self.alpha_df.loc[i_sc].values.tolist()
                        alphas = [alphas + [alphas[-1]], alphas]
                        filling = [band_filled_or_fixed + [False], band_filling]
                    else:
                        alphas = self.alpha_df.loc[i_sc]
                        filling = band_filled_or_fixed

                    # Set up calculator
                    calc = self.new_calculator('cp', calc_type, alphas=alphas,
                                               filling=filling, fixed_band=min(
                                                   fixed_band, n_filled_bands + 1),
                                               index_empty_to_save=index_empty_to_save,
                                               outdir=outdir_band)
                    calc.directory = directory

                    # Run cp.x
                    self.run_calculator(calc)

                    # Reset the value of 'fixed_band' so we can keep track of which calculation
                    # is which. This is important for empty orbital calculations, where fixed_band
                    # is always set to the LUMO but in reality we're fixing the band corresponding
                    # to index_empty_to_save from an earlier calculation
                    calc.fixed_band = fixed_band

                    # Store the result
                    # We store the results in one of two lists: alpha_indep_calcs and
                    # alpha_dep_calcs. The latter is overwritten at each new self-
                    # consistency loop.
                    if 'ki' in calc_type and 'print' not in calc_type:
                        alpha_dep_calcs.append(calc)
                    elif 'pbe' in calc_type and 'dummy' not in calc_type:
                        if self.functional in ['ki', 'pkipz']:
                            # For KI, the PBE calculations are independent of alpha so
                            # we store these in a list that is never overwritten
                            alpha_indep_calcs.append(calc)
                        else:
                            # For KIPZ, the PBE calculations are dependent on alpha via
                            # the definition of the variational orbitals. We only want to
                            # store the calculations that used the most recent value of alpha
                            alpha_dep_calcs.append(calc)

                    # Copying of evcfixed_empty.dat to evc_occupied.dat
                    if calc_type in ['pz_print', 'kipz_print']:
                        evcempty_dir = f'{outdir_band}/' \
                            f'{calc.prefix}_{calc.ndw}.save/K00001/'
                    elif calc_type == 'pbe_n+1_dummy':
                        evcocc_dir = f'{outdir_band}/' \
                            f'{calc.prefix}_{calc.ndr}.save/K00001/'
                        if os.path.isfile(f'{evcempty_dir}/evcfixed_empty.dat'):
                            utils.system_call(
                                f'cp {evcempty_dir}/evcfixed_empty.dat {evcocc_dir}/evc_occupied.dat')
                        else:
                            raise OSError(
                                f'Could not find {evcempty_dir}/evcfixed_empty.dat')

                if self.from_scratch:

                    # Calculate an updated alpha and a measure of the error
                    # E(N) - E_i(N - 1) - lambda^alpha_ii(1)     (filled)
                    # E_i(N + 1) - E(N) - lambda^alpha_ii(0)     (empty)
                    #
                    # Note that we only do this if the calculation has not been skipped
                    # because calculate_alpha reads in alpha values from files which get
                    # overwritten by subsequent calculations

                    calcs = [c for calc_set in [alpha_dep_calcs, alpha_indep_calcs]
                             for c in calc_set if c.fixed_band == fixed_band]

                    if self.functional in ['ki', 'pkipz']:
                        alpha, error = calculate_alpha(
                            calcs, filled=filled, kipz=False)
                    else:
                        alpha, error = calculate_alpha(
                            calcs, filled=filled, kipz=True)

                    for band_in_group in all_bands_in_group:
                        self.alpha_df.loc[i_sc + 1, band_in_group] = alpha
                    self.error_df.loc[i_sc, fixed_band] = error
                    self.alpha_df.to_pickle('alphas.pkl')
                    self.error_df.to_pickle('errors.pkl')

            self.print_summary()

            converged = all([abs(e) < 1e-3 or pd.isnull(e)
                             for e in self.error_df.loc[i_sc, :]])

        if converged:
            print('Alpha values have been converged')
        else:
            print('Alpha values have been determined but are not necessarily converged')

        # Writing alphas to a .tex table
        latex_table = self.alpha_df.to_latex(column_format='l' + 'd' * n_bands,
                                             float_format='{:.3f}'.format, escape=False)
        with open('tab_alpha_values.tex', 'w') as f:
            f.write(latex_table)

    def perform_final_calculations(self):

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
                ndr = [c.ndw for c in self.all_calcs if c.name in ['ki', 'ki_final']][-1]
                calc = self.new_calculator('cp', final_calc_type, ndr=ndr)
            else:
                calc = self.new_calculator('cp', final_calc_type)

            calc.directory = directory

            if self.init_variational_orbitals in ['mlwfs', 'projw']:
                calc.do_outerloop = False
                calc.do_outerloop_empty = False
                calc.write_hr = True

            self.run_calculator(calc)

    def new_calculator(self, calc_type, calc_presets='pbe_init', alphas=None, filling=None, **kwargs):
        """

        Generates a new calculator based on the self.master_calc[calc_type]
        template calculator

        """

        if calc_type == 'cp':
            return self.new_cp_calculator(calc_presets, alphas, filling, **kwargs)
        else:
            raise ValueError(f'You should not be requesting calculators with {calc_type} != "cp"')

    def new_cp_calculator(self, calc_presets='pbe_init', alphas=None, filling=None, **kwargs):
        """

        Generates a new CP calculator based on the self.master_calc["cp"]
        template calculator, modifying the appropriate settings to match the
        chosen calc_presets, and altering any Quantum Espresso keywords
        specified as kwargs

        Arguments:

            calc_presets
                The set of preset values to use; must be one of the following strings:

                Initialisation
                'pbe_init'            PBE calculation from scratch
                'pz_init'             PZ calculation starting from PBE restart
                'pz_innerloop_init'   PZ calculation starting from PBE restart (innerloop only)

                Trial calculations
                'ki'     KI calculation with N electrons and empty bands if specified
                'kipz'   As above, but for KIPZ

                For calculating alpha_i for filled orbitals.
                'pbe_frozen'    PBE calculation leaving rho unchanged, for reporting energy and lambda
                'ki_frozen'     KI calculation with N electrons for generating lambda only (will not alter density)
                'kipz_frozen'   KIPZ calculation with N electrons for generating lambda only (will not alter density)
                'pbe_n-1'       PBE calculation with N-1 electrons via fixed_state; rho is optimised
                'kipz_n-1'      KIPZ calculation with N-1 electrons via fixed_state; rho is optimised

                For calculating alpha_i for empty orbitals
                'pz_print'         PZ calculation that generates evcempty_fixed.dat file
                'kipz_print'       KIPZ calculation that generates evcempty_fixed.dat file
                'pbe_n+1_dummy'    PBE dummy calculation that generates store files of
                                   the correct dimensions
                'pbe_n+1'          PBE calculation with N+1 electrons; rho is optimised
                'kipz_n+1'         KIPZ calculation with N+1 electrons; rho is optimised
                'pbe_n+1-1_frozen' PBE calculation with N electrons, starting from N+1
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

        Returns: a new CP calculator object

        """

        # By default, use the last row in the alpha table for the screening parameters
        if alphas is None:
            alphas = self.alpha_df.iloc[-1]

        # Generate a new cp calculator copied from the master calculator
        calc = CP_calc(self.master_calcs["cp"], alphas=alphas, filling=filling)

        # Set up read/write indexes
        if calc_presets == 'pbe_init':
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
        elif calc_presets == 'pbe_frozen':
            ndr = 60
            ndw = 62
        elif calc_presets in ['pbe_n-1', 'kipz_n-1']:
            ndr = 60
            ndw = 63
        elif calc_presets in ['pz_print', 'kipz_print']:
            ndr = 60
            ndw = 64
        elif calc_presets == 'pbe_n+1_dummy':
            ndr = 65
            ndw = 65
        elif calc_presets in ['ki_n+1-1_frozen', 'kipz_n+1-1_frozen']:
            ndr = 65
            ndw = 66
        elif calc_presets == 'pbe_n+1-1_frozen':
            ndr = 65
            ndw = 67
        elif calc_presets in ['pbe_n+1', 'kipz_n+1']:
            ndr = 65
            ndw = 68
        elif calc_presets in ['ki_final', 'kipz_final']:
            ndr = 60
            ndw = 70
        elif calc_presets in ['pkipz_final']:
            ndr = 70
            ndw = 71
        else:
            raise ValueError('Invalid calc_presets "{}"'.format(calc_presets))

        # CP options
        # control
        calc.name = calc_presets
        calc.ndw = ndw
        calc.ndr = ndr
        if calc.name in ['pbe_init', 'pbe_n+1_dummy']:
            calc.restart_mode = 'from_scratch'
        else:
            calc.restart_mode = 'restart'

        # system
        if 'pz' in calc.name or 'ki' in calc.name:
            calc.do_orbdep = True
        else:
            calc.do_orbdep = False
        if calc.name in ['pbe_init', 'pz_init', 'pz_innerloop_init', 'ki',
                         'kipz', 'pz_print', 'kipz_print', 'pbe_n+1_dummy',
                         'ki_final', 'kipz_final', 'pkipz_final']:
            calc.fixed_state = False
        else:
            calc.fixed_state = True
            if '-1' in calc.name:
                calc.f_cutoff = 1e-5
            else:
                calc.f_cutoff = 1.0
        if 'n+1' in calc.name:
            calc.nelec += 1
            calc.nelup += 1
            if 'dummy' not in calc.name:
                calc.restart_from_wannier_pwscf = True
        # Peridoc calculations do not use complex wavefunctions (by default the wavefunctions
        # are real since they come from a Γ-only calculation and no counter-charge corrections
        # are applied), whereas aperiodic calculations benefit from using complex wavefunctions
        calc.do_wf_cmplx = not self.periodic

        # electrons
        # For all calculations calculating alpha, remove the empty states and
        # increase the energy thresholds
        if not any([s in calc.name for s in ['init', 'print', 'final']]) and calc.name not in ['ki', 'kipz']:
            calc.empty_states_nbnd = 0
            calc.conv_thr *= 100
            calc.esic_conv_thr *= 100

        if any([s in calc.name for s in ['frozen', 'dummy', 'print', 'innerloop']]) or calc.name == 'pkipz_final':
            calc.do_outerloop = False
            if calc.empty_states_nbnd > 0:
                calc.do_outerloop_empty = False
        elif calc.name in ['ki', 'ki_final']:
            calc.do_outerloop = False
            if calc.empty_states_nbnd > 0:
                calc.do_outerloop_empty = True
        else:
            calc.do_outerloop = True
            if calc.empty_states_nbnd > 0:
                calc.do_outerloop_empty = True
        if calc.maxiter is None and calc.do_outerloop:
            calc.maxiter = 300
        if calc.empty_states_maxstep is None and calc.do_outerloop_empty:
            calc.empty_states_maxstep = 300

        # nksic
        if calc.do_orbdep:
            calc.odd_nkscalfact = True
            calc.odd_nkscalfact_empty = True
        calc.do_innerloop_cg = True
        if calc.name[:2] == 'pz' and 'print' not in calc.name:
            calc.do_innerloop = True
        else:
            calc.do_innerloop = False
        if calc.empty_states_nbnd > 0:
            calc.do_innerloop_empty = False
        if 'kipz' in calc.name:
            calc.which_orbdep = 'nkipz'
        elif 'pz' in calc.name:
            calc.which_orbdep = 'pz'
        elif 'ki' in calc.name:
            calc.which_orbdep = 'nki'
        if 'print' in calc.name:
            calc.print_wfc_anion = True

        if self.periodic:
            calc.which_compensation = 'none'
        else:
            calc.which_compensation = 'tcc'

        # Handle any keywords provided by kwargs
        # Note that since this is performed after the above logic this can (deliberately
        # or accidentally) overwrite the above settings

        for keyword, value in kwargs.items():
            setattr(calc, keyword, value)

        # Sanity checking
        if calc.print_wfc_anion and calc.index_empty_to_save is None:
            raise ValueError('Error: print_wfc_anion is set to true but you have not selected '
                             ' an index_empty_to_save. Provide this as an argument to set_calculator')

        if calc.fixed_band is not None and calc.fixed_band > calc.nelup + 1:
            utils.warn(
                'calc.fixed_band is higher than the LUMO; this should not happen')

        # avoid innerloops for one-orbital-manifolds
        if calc.nelup in [0, 1] and calc.neldw in [0, 1]:
            calc.do_innerloop = False
        if calc.empty_states_nbnd == 1:
            calc.do_innerloop_empty = False

        # don't print QC in some cases
        if 'dummy' in calc.name:
            calc.skip_qc = True

        return calc

    def print_summary(self):
        # Printing out a progress summary
        print('\nalpha')
        print(self.alpha_df)
        if not self.error_df.empty:
            print('\ndelta E - lambda^alpha_ii (eV)')
            print(self.error_df)
        print('')


def calculate_alpha(calcs, filled=True, kipz=False):
    '''

    Calculates alpha via equation 10 of Nguyen et. al (2018) 10.1103/PhysRevX.8.021051
    If the band is filled, use s = 1; if the band is empty, use s = 0

    Arguments:
        calcs  -- a list of calculations
        filled -- True if the orbital for which we're calculating alpha is filled
        kipz   -- True if KIPZ; False if KI

    '''

    if kipz and filled:
        # KIPZ N
        [alpha_calc] = [c for c in calcs if c.which_orbdep == 'nkipz'
                        and c.do_orbdep and c.f_cutoff == 1.0]
        kipz = alpha_calc.results

        # KIPZ N-1
        [kipz_m1] = [c.results for c in calcs if c.which_orbdep == 'nkipz'
                     and c.do_orbdep and c.f_cutoff < 0.0001]

        # PBE N
        [pbe] = [c.results for c in calcs if not c.do_orbdep
                 and c.restart_mode == 'restart' and c.f_cutoff == 1.0]

        dE = kipz['energy'] - kipz_m1['energy']
        lambda_a = kipz['lambda_ii']
        lambda_0 = pbe['lambda_ii']

    elif kipz and not filled:
        # KIPZ N+1-1
        [alpha_calc] = [c for c in calcs if c.which_orbdep == 'nkipz'
                        and c.do_orbdep and c.f_cutoff < 0.0001]
        kipz = alpha_calc.results

        # KIPZ N+1
        [kipz_p1] = [c.results for c in calcs if c.which_orbdep == 'nkipz'
                     and c.do_orbdep and c.f_cutoff == 1.0]

        # PBE N+1-1
        [pbe] = [c.results for c in calcs if not c.do_orbdep
                 and c.restart_mode == 'restart' and c.f_cutoff < 0.0001]

        dE = kipz_p1['energy'] - kipz['energy']
        lambda_a = kipz['lambda_ii']
        lambda_0 = pbe['lambda_ii']

    elif not kipz and filled:
        # KI N
        [alpha_calc] = [c for c in calcs if c.which_orbdep == 'nki'
                        and c.do_orbdep and c.f_cutoff == 1.0]
        ki = alpha_calc.results

        # PBE N-1
        [pbe_m1] = [c.results for c in calcs if not c.do_orbdep
                    and c.restart_mode == 'restart' and c.f_cutoff < 0.0001]

        # PBE N
        [pbe] = [c.results for c in calcs if not c.do_orbdep
                 and c.restart_mode == 'restart' and c.f_cutoff == 1.0]

        dE = pbe['energy'] - pbe_m1['energy']
        lambda_a = ki['lambda_ii']
        lambda_0 = pbe['lambda_ii']

    else:
        # KI N+1-1
        [alpha_calc] = [c for c in calcs if c.which_orbdep == 'nki'
                        and c.do_orbdep and c.f_cutoff < 0.0001]
        ki = alpha_calc.results

        # PBE N+1
        [pbe_p1] = [c.results for c in calcs if not c.do_orbdep
                    and c.restart_mode == 'restart' and c.f_cutoff == 1.0]

        # PBE N+1-1
        [pbe] = [c.results for c in calcs if not c.do_orbdep
                 and c.restart_mode == 'restart' and c.f_cutoff < 0.0001]

        dE = pbe_p1['energy'] - pbe['energy']
        lambda_a = ki['lambda_ii']
        lambda_0 = pbe['lambda_ii']

    if not kipz:
        # Check total energy is unchanged
        dE_check = abs(ki['energy'] - pbe['energy'])
        if dE_check > 1e-5:
            warn('KI and PBE energies differ by {:.5f} eV'.format(dE_check))

    # Obtaining alpha
    if (alpha_calc.odd_nkscalfact and filled) or (alpha_calc.odd_nkscalfact_empty and not filled):
        alpha_guesses = alpha_calc.alphas[0]
        alpha_guess = alpha_guesses[alpha_calc.fixed_band - 1]
    else:
        alpha_guess = alpha_calc.nkscalfact

    alpha = alpha_guess * (dE - lambda_0) / (lambda_a - lambda_0)

    # The error is lambda^alpha(1) - lambda^alpha_i(1)
    error = dE - lambda_a

    return alpha, error
