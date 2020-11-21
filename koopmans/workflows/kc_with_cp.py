"""

Workflow module for python_KI, containing the workflow for performing KI and KIPZ calculations

Written by Edward Linscott Jan 2020
Split off from workflow.py Oct 2020

"""

import os
import subprocess
import copy
import pandas as pd
import koopmans.utils as utils
from koopmans.io import write_alpharef, read_alpharef, print_summary, print_qc
from koopmans.calculators.calculator import run_qe, calculate_alpha
from koopmans.calculators.cp import CP_calc
import ipdb


def run(workflow_settings, calcs_dct):
    '''
    This function runs the KI/KIPZ workflow from start to finish

    Arguments:
        workflow_settings -- a dictionary containing workflow settings
        calcs_dct         -- a dictionary of calculators (one per code e.g. cp.x, w90, ...)

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

    from koopmans.config import from_scratch

    functional = workflow_settings['functional']

    if 'cp' not in calcs_dct:
        raise ValueError(
            'Performing a KC calculation requires a "cp" block in the .json input file')
    master_calc = copy.deepcopy(calcs_dct['cp'])

    # Sanitise outdir
    master_calc.outdir = os.getcwd() + '/' + master_calc.outdir.strip('./')

    # Removing old directories
    if from_scratch:
        if workflow_settings['init_density'] != 'ki':
            # if init_density == "ki" we don't want to delete the directory containing
            # the KI calculation we're reading the manifold from, or the TMP files
            utils.system_call('rm -r init 2>/dev/null', False)
            utils.system_call(f'rm -r {master_calc.outdir} 2>/dev/null', False)
        utils.system_call('rm -r calc_alpha 2>/dev/null', False)
        utils.system_call('rm -r final 2>/dev/null', False)

    # Counting the number of bands
    n_filled_bands = master_calc.nelup
    n_empty_bands = master_calc.empty_states_nbnd
    if n_empty_bands is None:
        n_empty_bands = 0
    n_bands = n_filled_bands + n_empty_bands
    band_filling = [True for _ in range(
        n_filled_bands)] + [False for _ in range(n_empty_bands)]
    i_bands = range(1, n_bands + 1)

    # Preparing panda dataframes in which to store results
    alpha_df = pd.DataFrame(columns=i_bands)
    error_df = pd.DataFrame(columns=i_bands)
    if workflow_settings['alpha_from_file']:
        print('Reading alpha values from file')
        alpha_df.loc[1] = read_alpharef(directory='.')
    else:
        alpha_df.loc[1] = [workflow_settings['alpha_guess']
                           for _ in range(n_bands)]

    print('\nINITIALISATION OF DENSITY')

    init_density = workflow_settings['init_density']
    if init_density == 'pbe':
        calc = set_up_calculator(master_calc, 'pbe_init')
        calc.directory = 'init'
        run_qe(calc, silent=False,
               enforce_ss=workflow_settings['enforce_spin_symmetry'])

    elif init_density == 'pbe-pw':
        # PBE using pw.x
        raise ValueError('init_denisty: pbe-pw is not yet implemented')

    elif init_density == 'pz':
        # PBE from scratch
        calc = set_up_calculator(master_calc, 'pbe_init')
        calc.directory = 'init'
        run_qe(calc, silent=False)
        # PZ from PBE
        calc = set_up_calculator(master_calc, 'pz_init')
        calc.directory = 'init'
        write_alpharef(alpha_df.loc[1], calc)
        run_qe(calc, silent=False)

    elif init_density == 'ki':
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
        utils.system_call(f'mv {old_savedir} {savedir}')

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
    else:
        raise ValueError(
            "Should not arrive here; compare the above code with workflow.valid_settings")

    if init_density == 'pbe':
        # Using KS eigenfunctions as guess variational orbitals
        print('Overwriting the CP variational orbitals with Kohn-Sham orbitals')
        savedir = f'{calc.outdir}/{calc.prefix}_{calc.ndw}.save/K00001'
        utils.system_call(f'cp {savedir}/evc1.dat {savedir}/evc01.dat')
        utils.system_call(f'cp {savedir}/evc2.dat {savedir}/evc02.dat')
        if calc.empty_states_nbnd is not None and calc.empty_states_nbnd > 0:
            utils.system_call(
                f'cp {savedir}/evc_empty1.dat {savedir}/evc0_empty1.dat')
            utils.system_call(
                f'cp {savedir}/evc_empty2.dat {savedir}/evc0_empty2.dat')

    print('\nINITIALISATION OF MANIFOLD')

    init_manifold = workflow_settings['init_manifold']
    # Since nelec is defined automatically, record whether or not the previous calculation had 2 electrons
    # before calc gets overwritten
    calc_has_two_electrons = (calc.nelec == 2)
    if init_manifold == 'pz':
        calc = set_up_calculator(master_calc, 'pz_innerloop_init')
    elif init_manifold == 'ki':
        if init_density != 'ki':
            raise ValueError('Initialising manifold with KI makes no sense unless '
                             'reading from a pre-existing KI calculation')
        print('Copying the density from a pre-existing KI calculation')
    else:
        raise ValueError(
            f'Unrecognised option "{init_manifold}" for init_manifold. Should be one of "pz"/"ki"/"skip"')

    calc.directory = 'init'
    write_alpharef(alpha_df.loc[1], calc)

    if init_manifold == 'ki':
        pass
    elif calc_has_two_electrons and calc.one_innerloop_only and (init_manifold == 'pz' or calc.empty_states_nbnd == 0):
        # If we only have two electrons, then the filled manifold is trivially invariant under unitary
        # transformations. Likewise, if we have no empty states or if we're using a functional which is
        # invariant w.r.t. unitary rotations of the empty states, then the empty manifold need not be minimised
        # In these instances, we can skip the initialisation of the manifold entirely
        print('Skipping the optimisation of the manifold since it is invariant under unitary transformations')
        save_prefix = f'{calc.outdir}/{calc.prefix}'
        utils.system_call(
            f'cp -r {save_prefix}_{calc.ndr}.save {save_prefix}_{calc.ndw}.save')
    elif init_manifold in [init_density, 'skip']:
        if init_manifold == 'skip':
            print('Skipping the optimisation of the manifold')
        else:
            print('Skipping the optimisation of the manifold since it was already optimised during the density initialisation')
        save_prefix = f'{calc.outdir}/{calc.prefix}'
        utils.system_call(
            f'cp -r {save_prefix}_{calc.ndr}.save {save_prefix}_{calc.ndw}.save')
    else:
        run_qe(calc, silent=False)

    if from_scratch and init_manifold != 'ki' and workflow_settings['enforce_spin_symmetry']:
        print('Copying the spin-up variational orbitals over to the spin-down channel')
        savedir = f'{calc.outdir}/{calc.prefix}_{calc.ndw}.save/K00001'
        utils.system_call(f'cp {savedir}/evc01.dat {savedir}/evc02.dat')
        if calc.empty_states_nbnd is not None and calc.empty_states_nbnd > 0:
            utils.system_call(
                f'cp {savedir}/evc0_empty1.dat {savedir}/evc0_empty2.dat')

    print('\nDETERMINING ALPHA VALUES')
    if workflow_settings['calculate_alpha']:
        # Group the bands by group, and work out which bands to solve explicitly
        groups = workflow_settings['orbital_groups']
        if groups is None:
            groups = range(len(i_bands))
        bands_to_solve = {}
        groups_found = set([])
        for i_band in list(range(1, n_filled_bands + 1)[::-1]) \
                + list(range(n_filled_bands + 1, n_bands + 1)):
            # Looping through the filled bands from highest to lowest, then empty bands from
            # lowest to highest
            group = groups[i_band - 1]
            if group not in groups_found:
                groups_found.add(group)
                bands_to_solve[i_band] = [i + 1 for i,
                                          g in enumerate(groups) if g == group]
        if groups_found != set(groups):
            raise ValueError('Splitting of orbitals into groups failed')

        # Set up directories
        if not os.path.isdir('calc_alpha'):
            utils.system_call('mkdir calc_alpha')
        # for i_band in bands_to_solve:
        #     if not os.path.isdir(f'calc_alpha/orbital_{i_band}'):
        #         os.system('cp -r init calc_alpha/orbital_{}'.format(i_band))

        converged = False
        i_sc = 0

        if not from_scratch and os.path.isfile('alphas.pkl'):
            # Reloading alphas and errors from file
            print('Reloading alpha values from file')
            alpha_df = pd.read_pickle('alphas.pkl')
            error_df = pd.read_pickle('errors.pkl')

        alpha_indep_calcs = []

        while not converged and i_sc < workflow_settings['n_max_sc_steps']:
            i_sc += 1
            iteration_directory = 'calc_alpha'
            _, outdir = master_calc.outdir.rsplit('/', 1)
            outdir = os.getcwd() + f'/{iteration_directory}/{outdir}'

            if not os.path.isdir(outdir):
                utils.system_call(f'mkdir {outdir}')

            # Setting up directories
            if workflow_settings['n_max_sc_steps'] > 1:
                print('\n== SC iteration {} ==================='.format(i_sc))
                iteration_directory += f'/iteration_{i_sc}'
                if not os.path.isdir(iteration_directory):
                    utils.system_call(f'mkdir {iteration_directory}')

            # Do a KI/KIPZ calculation with the updated alpha values
            calc = set_up_calculator(
                master_calc, calc_type=functional.replace('pkipz', 'ki'))
            calc.directory = iteration_directory
            write_alpharef(alpha_df.loc[i_sc], calc)

            enforce_ss = workflow_settings['enforce_spin_symmetry']
            if i_sc == 1:
                enforce_ss = False
                if functional == 'kipz':
                    # For the first KIPZ trial calculation, do the innerloop
                    calc.do_innerloop = True
            else:
                # For later SC loops, read in the matching calculation from the
                # previous loop rather than the initialisation calculations
                calc.ndr = calc.ndw

            # Run the calculation and store the result
            run_qe(calc, silent=False)
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
                        f'ln -sr {master_calc.outdir}/*.save {outdir_band}')

                # Don't repeat if this particular alpha_i was converged
                if i_sc > 1 and any([abs(e) < workflow_settings['alpha_conv_thr'] for e in
                                     error_df.loc[:i_sc - 1, fixed_band]]):
                    print(
                        f'Skipping band {fixed_band} since this alpha is already '
                        'converged')
                    alpha_df.loc[i_sc + 1,
                                 all_bands_in_group] = alpha_df.loc[i_sc, fixed_band]
                    error_df.loc[i_sc,
                                 all_bands_in_group] = error_df.loc[i_sc - 1, fixed_band]
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
                if functional in ['ki', 'pkipz']:
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
                    if functional in ['ki', 'pkipz']:
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
                        write_alpharef(alpha_df.loc[i_sc], None,
                                       band_filling, directory)
                    elif not filled:
                        # In the case of empty orbitals, we gain an extra orbital in
                        # the spin-up channel
                        alpha_padded = list(alpha_df.loc[i_sc])
                        alpha_padded += [alpha_padded[-1]] + alpha_padded
                        write_alpharef(alpha_padded, None, band_filled_or_fixed + [False]
                                       + band_filling, directory, duplicate=False)
                    else:
                        write_alpharef(alpha_df.loc[i_sc], None,
                                       band_filled_or_fixed, directory)

                    # Set up calculator
                    calc = set_up_calculator(master_calc, calc_type,
                                             fixed_band=min(
                                                 fixed_band, n_filled_bands + 1),
                                             index_empty_to_save=index_empty_to_save,
                                             outdir=outdir_band)
                    calc.directory = directory

                    # Run cp.x
                    run_qe(calc, silent=False)

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
                        if functional in ['ki', 'pkipz']:
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

                if from_scratch:

                    # Calculate an updated alpha and a measure of the error
                    # E(N) - E_i(N - 1) - lambda^alpha_ii(1)     (filled)
                    # E_i(N + 1) - E(N) - lambda^alpha_ii(0)     (empty)
                    #
                    # Note that we only do this if the calculation has not been skipped
                    # because calculate_alpha reads in alpha values from files which get
                    # overwritten by subsequent calculations

                    calcs = [c for calc_set in [alpha_dep_calcs, alpha_indep_calcs]
                             for c in calc_set if c.fixed_band == fixed_band]

                    if functional in ['ki', 'pkipz']:
                        alpha, error = calculate_alpha(
                            calcs, filled=filled, kipz=False)
                    else:
                        alpha, error = calculate_alpha(
                            calcs, filled=filled, kipz=True)

                    for band_in_group in all_bands_in_group:
                        alpha_df.loc[i_sc + 1, band_in_group] = alpha
                    error_df.loc[i_sc, fixed_band] = error
                    alpha_df.to_pickle('alphas.pkl')
                    error_df.to_pickle('errors.pkl')

            print_summary(alpha_df, error_df)

            converged = all([abs(e) < 1e-3 or pd.isnull(e)
                             for e in error_df.loc[i_sc, :]])

        if converged:
            print('Alpha values have been converged')
        else:
            print('Alpha values have been determined but are not necessarily converged')

        # Writing alphas to a .tex table
        latex_table = alpha_df.to_latex(column_format='l' + 'd' * n_bands,
                                        float_format='{:.3f}'.format, escape=False)
        with open('tab_alpha_values.tex', 'w') as f:
            f.write(latex_table)
    else:
        print('Skipping calculation of screening parameters; reading values from file')
        alpha_df.loc[2] = read_alpharef(directory='.')
        print_summary(alpha_df)

    # Final calculation
    print(f'\nFINAL {functional.upper().replace("PK","pK")} CALCULATION')

    directory = 'final'
    if not os.path.isdir(directory):
        utils.system_call(f'mkdir {directory}')

    if functional == 'pkipz':
        final_calc_types = ['ki', 'pkipz']
    else:
        final_calc_types = [functional]

    for final_calc_type in final_calc_types:

        # If we performed the alpha calculations, direct the calculator
        # to restart from them
        if workflow_settings['calculate_alpha'] or final_calc_type == 'pkipz':
            final_calc_type += '_final'

        # For pKIPZ, the appropriate ndr can change but it is always ndw of the previous
        # KI calculation
        if final_calc_type == 'pkipz_final':
            ndr = calc.ndw
            calc = set_up_calculator(master_calc, final_calc_type, ndr=ndr,
                                     empty_states_nbnd=n_empty_bands)
        else:
            calc = set_up_calculator(master_calc, final_calc_type,
                                     empty_states_nbnd=n_empty_bands)

        calc.directory = directory
        write_alpharef(alpha_df.iloc[-1], calc)

        run_qe(calc, silent=False)

    # Print out data for quality control
    if workflow_settings['print_qc']:
        print('\nQUALITY CONTROL')
        print_qc('energy', calc.results['energy'])
        for i, alpha in enumerate(alpha_df.iloc[-1]):
            print_qc(f'alpha({i})', alpha)
        for isp, orbs_self_hartree in enumerate(calc.results['orbital_data']['self-Hartree']):
            for i, orb_sh in enumerate(orbs_self_hartree):
                print_qc(f'orbital_self_Hartree(orb={i},sigma={isp})', orb_sh)
        print_qc('HOMO', calc.results['homo_energy'])
        if calc.results['lumo_energy'] is not None:
            print_qc('LUMO', calc.results['lumo_energy'])

    print('\nWORKFLOW COMPLETE')

    return calc


def set_up_calculator(calc, calc_type='pbe_init', **kwargs):
    """
    Generates a new calculator based on a template calculator, modifying
    the appropriate settings to match the chosen calc_type, and altering any
    Quantum Espresso keywords specified as kwargs

    Arguments:

        calc: an template calculator. N.B. this object is not modified by
              this routine

        calc_type: the calculation type; must be one of

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

       **kwargs: accepts any Quantum Espresso keywords as an argument, and will
                 apply these options to the returned ASE calculator

    Returns: a new ASE calculator object

    """

    # Avoid modifying the master calculator
    calc = copy.deepcopy(calc)

    # Set up read/write indexes
    if calc_type == 'pbe_init':
        ndr = 50
        ndw = 50
    elif calc_type in ['pz_init', 'pz_innerloop_init']:
        ndr = 50
        ndw = 51
    elif calc_type in ['ki', 'kipz']:
        ndr = 51
        ndw = 60
    elif calc_type in ['ki_frozen', 'kipz_frozen']:
        ndr = 60
        ndw = 61
    elif calc_type == 'pbe_frozen':
        ndr = 60
        ndw = 62
    elif calc_type in ['pbe_n-1', 'kipz_n-1']:
        ndr = 60
        ndw = 63
    elif calc_type in ['pz_print', 'kipz_print']:
        ndr = 60
        ndw = 64
    elif calc_type == 'pbe_n+1_dummy':
        ndr = 65
        ndw = 65
    elif calc_type in ['ki_n+1-1_frozen', 'kipz_n+1-1_frozen']:
        ndr = 65
        ndw = 66
    elif calc_type == 'pbe_n+1-1_frozen':
        ndr = 65
        ndw = 67
    elif calc_type in ['pbe_n+1', 'kipz_n+1']:
        ndr = 65
        ndw = 68
    elif calc_type in ['ki_final', 'kipz_final']:
        ndr = 60
        ndw = 70
    elif calc_type in ['pkipz_final']:
        ndr = 70
        ndw = 71
    else:
        raise ValueError('Invalid calc_type "{}"'.format(calc_type))

    # CP options
    # control
    calc.name = calc_type
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

    # if periodic, set calc.which_compensation = 'none'
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

    return calc


keywords_altered_during_workflow = ['ndw', 'ndr', 'restart_mode', 'nelec', 'nelup',
                                    'neldw', 'do_orbdep', 'fixed_state', 'f_cutoff',
                                    'fixed_band', 'conv_thr', 'restart_from_wannier_pwscf',
                                    'maxiter', 'empty_states_maxstep', 'esic_conv_thr',
                                    'do_innerloop', 'freeze_density', 'which_orbdep',
                                    'print_wfc_anion', 'odd_nkscalfact',
                                    'odd_nkscalfact_empty', 'nkscalfact', 'do_innerloop_empty',
                                    'innerloop_nmax']
