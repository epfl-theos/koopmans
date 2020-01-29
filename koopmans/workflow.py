"""

Workflow module for python_KI, containing the workflow for performing KI and KIPZ calculations

Written by Edward Linscott Jan 2020

"""

import os
import copy
import pickle
import pandas as pd
import warnings
from ase.io import espresso_cp as cp_io
from koopmans.io import write_alpharef, read_alpharef, print_summary
from koopmans.calc import run_cp, calculate_alpha, Extended_Espresso_cp, calc_from_json


def set_up_calculator(calc, calc_type='pbe_init', **kwargs):
    """
    Generates a new ASE calculator based on a template ASE calculator, modifying
    the appropriate settings to match the chosen calc_type, and altering any 
    Quantum Espresso keywords specified as kwargs

    Arguments:

        calc: an template ASE calculator. N.B. this object is not modified by
              this routine

        calc_type: the calculation type; must be one of

            Initialisation
            'pbe_init'  PBE calculation from scratch
            'pz'        PZ calculation starting from PBE restart
            'kipz_init' KIPZ starting from PBE restart

            For calculating alpha_i for filled orbitals
            'pbe'      PBE calculation starting from restart
            'pbe_n-1'  PBE calculation with N-1 electrons via fixed_state
            'kipz_n-1' KIPZ calculation with N-1 electrons via fixed_state
            'ki'       KI calculation with N electrons
            'kipz'     KIPZ calculation with N electrons

            For calculating alpha_i for empty orbitals
            'pz_print' PZ calculation that generates evcempty_fixed.dat file
            'kipz_print'    KIPZ calculation that generates evcempty_fixed.dat file
            'pbe_n+1_dummy' PBE dummy calculation that generates store files of
                            the correct dimensions
            'pbe_n+1'       PBE calculation with N+1 electrons
            'kipz_n+1'      KIPZ calculation with N+1 electrons
            'pbe_n+1-1'     PBE calculation with N electrons, starting from N+1 
                            but with f_cutoff = 0.00001
            'ki_n+1-1'      KI calculation with N electrons, starting from N+1 
                            but with f_cutoff = 0.00001
            'kipz_n+1-1'    KIPZ calculation with N electrons, starting from N+1 
                            but with f_cutoff = 0.00001

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
    elif calc_type in ['pz', 'kipz_init']:
        ndr = 50
        ndw = 51
    elif calc_type == 'pbe':
        ndr = 51
        ndw = 52
    elif calc_type == 'pbe_n-1':
        ndr = 51
        ndw = 53
    elif calc_type in ['pz_print', 'kipz_print']:
        ndr = 51
        ndw = 54
    elif calc_type == 'pbe_n+1_dummy':
        ndr = 55
        ndw = 55
    elif calc_type == 'pbe_n+1-1':
        ndr = 55
        ndw = 56
    elif calc_type == 'pbe_n+1':
        ndr = 55
        ndw = 57
    elif calc_type in ['ki', 'kipz']:
        ndr = 51
        ndw = 60
    elif calc_type == 'kipz_n-1':
        ndr = 51
        ndw = 70
    elif calc_type == 'kipz_n+1':
        ndr = 55
        ndw = 80
    elif calc_type in ['ki_n+1-1', 'kipz_n+1-1']:
        ndr = 55
        ndw = 90
    else:
        raise ValueError('Invalid calc_type "{}"'.format(calc_type))

    # Pseudopotentials

    # This script will...
    #  1. try to locating the directory as currently specified by the calculator
    #  2. if that fails, it will check if $ESPRESSO_PSEUDO is set
    #  3. if that fails, it will raise an OS error

    # Therefore, if you want to use a pseudo directory that is different to your
    # $PSEUDO_DIR, reset pseudo_dir for the template calc object before entering
    # this routine.

    if calc.pseudo_dir is None or not os.path.isdir(calc.pseudo_dir):
        try:
            calc.pseudo_dir = os.environ.get('ESPRESSO_PSEUDO')
        except:
            raise OSError('Directory for pseudopotentials not found. Please define '
                          'the environment variable ESPRESSO_PSEUDO or provide the function run_cp '
                          'with a pseudo_dir argument')

    # CP options
    # control
    calc.setattr_only('prefix', calc_type)
    calc.ndw = ndw
    calc.ndr = ndr
    if calc.prefix in ['pbe_init', 'pbe_n+1_dummy']:
        calc.restart_mode = 'from_scratch'
    else:
        calc.restart_mode = 'restart'

    # system
    calc.nspin = 2
    if 'pz' in calc.prefix or 'ki' in calc.prefix:
        calc.do_orbdep = True
    else:
        calc.do_orbdep = False
    if calc.prefix in ['pbe_init', 'pz_print', 'pbe_n+1_dummy', 'pz', 'kipz_init', 'kipz_print']:
        calc.fixed_state = False
    else:
        calc.fixed_state = True
        if '-1' in calc.prefix:
            calc.f_cutoff = 1e-5
        else:
            calc.f_cutoff = 1.0
    if 'n+1' in calc.prefix:
        calc.nelec += 1
        calc.nelup += 1
    if calc.prefix in ['pbe_n+1', 'pbe_n+1-1', 'ki_n+1-1', 'kipz_n+1-1', 'kipz_n+1']:
        calc.restart_from_wannier_pwscf = True

    # electrons
    calc.conv_thr = calc.nelec*1e-9
    if 'print' in calc.prefix or (('pz' in calc.prefix or 'ki' in calc.prefix) and 'kipz' not in calc.prefix):
        calc.maxiter = 2
        calc.empty_states_maxstep = 1
    else:
        calc.maxiter = 200
        calc.empty_states_maxstep = 300

    # nksic
    calc.esic_conv_thr = calc.nelec*1e-9
    calc.do_innerloop_cg = True
    if calc.prefix[:2] == 'pz':
        calc.do_innerloop = True
        calc.one_innerloop_only = True
    else:
        calc.do_innerloop = False
        calc.one_innerloop_only = False
    if 'kipz' in calc.prefix:
        calc.which_orbdep = 'nkipz'
    elif calc.prefix == 'pz':
        calc.which_orbdep = 'pz'
    elif 'ki' in calc.prefix:
        calc.which_orbdep = 'nki'
    if 'print' in calc.prefix:
        calc.print_wfc_anion = True

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
        warnings.warn(
            'calc.fixed_band is higher than the LUMO; this should not happen')

    return calc


def run(workflow_type, json, alpha_guess=0.6, alpha_from_file=False, n_max_sc_steps=1, from_scratch=False):
    '''
    This function runs the KI/KIPZ workflow from start to finish

    Arguments:
        workflow_type       -- the workflow type: must be one of 'ki', 'kipz'
        json                -- the path to the json file containing the system data
                               (atomic positions, number of atoms etc.)
        alpha_guess         -- the initial guess for alpha (a single float for all orbitals)
        alpha_from_file     -- read alpha from pre-existing file_alpharef.txt
        n_max_sc_steps      -- the maximum number of self-consistent steps for 
                               determining {alpha_i}

    Running this function will generate a number of files:
       init/ -- the PBE and PZ calculations for initialisation
       calc_alpha/orbital_#/ -- PBE/PZ/KI calculations where we have fixed a particular orbital
       final/ -- the directory containing the final KI calculation
       alphas.pkl -- a python pickle file containing the alpha values
       errors.pkl -- a python pickle file containing the errors in the alpha values
       tab_alpha_values.tex -- a latex table of the alpha values
    '''

    workflow_type = workflow_type.lower()
    if workflow_type not in ['ki', 'kipz']:
        raise ValueError(
            f'Unrecognised calculation type {workflow_type} provided to koopmans.workflow.run()')

    # Removing old directories
    if from_scratch:
        os.system('rm -r init 2>/dev/null')
        os.system('rm -r calc_alpha 2>/dev/null')
        os.system('rm -r final 2>/dev/null')

    # Reading in JSON file
    master_calc = calc_from_json(json)

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

    print('\nINITIALISATION OF DENSITY')
    # PBE from scratch
    calc = set_up_calculator(master_calc, 'pbe_init')
    calc.directory = 'init'
    run_cp(calc, silent=False, from_scratch=from_scratch)

    # Moving orbitals
    print('Overwriting the CP variational orbitals with Kohn-Sham orbitals')
    savedir = f'{calc.directory}/TMP-CP/pbe_50.save/K00001'
    os.system(f'cp {savedir}/evc1.dat {savedir}/evc01.dat')
    os.system(f'cp {savedir}/evc2.dat {savedir}/evc02.dat')

    print('\nINITIALISATION OF MANIFOLD')
    # PZ/KIPZ reading in PBE to define manifold
    if alpha_from_file:
        print(r'Reading alpha values from file')
        alpha_df.loc[1] = read_alpharef(directory='.')
    else:
        alpha_df.loc[1] = [alpha_guess for _ in range(n_bands)]

    if workflow_type == 'ki':
        calc = set_up_calculator(
            master_calc, 'pz', empty_states_nbnd=n_empty_bands, from_scratch=from_scratch)
    else:
        calc = set_up_calculator(master_calc, 'kipz_init',
                                 odd_nkscalfact=True, odd_nkscalfact_empty=True)

    calc.directory = 'init'
    write_alpharef(alpha_df.loc[1], band_filling, calc.directory)
    prev_calc_not_skipped = run_cp(calc, silent=False)

    if prev_calc_not_skipped:
        prefix = calc.parameters['input_data']['control']['prefix']
        print('Copying the spin-up variational orbitals over to the spin-down channel')
        os.system(f'cp {calc.directory}/{calc.outdir}/{prefix}_{calc.ndw}.save/K00001/evc01.dat '
                  f'{calc.directory}/{calc.outdir}/{prefix}_{calc.ndw}.save/K00001/evc02.dat')

    print('\nDETERMINING ALPHA VALUES')

    # Set up directories
    if not os.path.isdir(f'calc_alpha'):
        os.system('mkdir calc_alpha')
    for i_band in i_bands:
        if not os.path.isdir(f'calc_alpha/orbital_{i_band}'):
            os.system('cp -r init calc_alpha/orbital_{}'.format(i_band))

    converged = False
    i_sc = 0

    if not prev_calc_not_skipped:
        # Reloading from file
        alpha_df = pd.read_pickle('alphas.pkl')
        error_df = pd.read_pickle('errors.pkl')
        i_sc = len(error_df)
        print_summary(alpha_df, error_df)
        converged = all([abs(e) < 1e-3 for e in error_df.loc[i_sc, :]])

    alpha_indep_calcs = []

    while not converged and i_sc < n_max_sc_steps:
        i_sc += 1
        alpha_dep_calcs = []

        if n_max_sc_steps > 1:
            print('\n== SC iteration {} ==================='.format(i_sc))

        # Loop over removing an electron from each band
        for fixed_band, filled in zip(i_bands, band_filling):
            print('-- Orbital {} ------------------------'.format(fixed_band))
            directory = 'calc_alpha/orbital_{}'.format(fixed_band)

            # Don't repeat if this particular alpha_i was converged
            if any([abs(e) < 1e-3 for e in error_df.loc[:, fixed_band]]):
                print(
                    'Skipping band {} since this alpha is already converged'.format(fixed_band))
                alpha_df.loc[i_sc + 1,
                             fixed_band] = alpha_df.loc[i_sc, fixed_band]
                error_df.loc[i_sc,
                             fixed_band] = error_df.loc[i_sc - 1, fixed_band]
                continue

            # Write/update the alpharef files in the work directory
            # Make sure to include the fixed band alpha in file_alpharef.txt
            # rather than file_alpharef_empty.txt
            band_filled_or_fixed = [
                b or i == fixed_band - 1 for i, b in enumerate(band_filling)]

            # Perform the fixed-band-dependent calculations
            if workflow_type == 'ki':
                if filled:
                    calc_types = ['pbe', 'pbe_n-1', 'ki']
                else:
                    calc_types = ['pz_print', 'pbe_n+1_dummy', 'pbe_n+1',
                                  'pbe_n+1-1', 'ki_n+1-1']
            else:
                if filled:
                    calc_types = ['kipz', 'pbe', 'kipz_n-1']
                else:
                    calc_types = ['kipz_print', 'pbe_n+1_dummy', 'kipz_n+1',
                                  'pbe_n+1-1', 'kipz_n+1-1']

            for calc_type in calc_types:
                if workflow_type == 'ki':
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

                if prev_calc_not_skipped:
                    if 'print' in calc_type:
                        # Note that the 'print' calculations for empty bands do not
                        # in fact involve the fixing of that band (and thus for the
                        # 'fixed' band the corresponding alpha should be in
                        # file_alpharef_empty.txt)
                        write_alpharef(alpha_df.loc[i_sc],
                                       band_filling, directory)
                    else:
                        write_alpharef(alpha_df.loc[i_sc],
                                       band_filled_or_fixed, directory)

                if filled:
                    index_empty_to_save = None
                else:
                    index_empty_to_save = fixed_band - n_filled_bands

                if 'ki' in calc_type:
                    odd_nkscalfact = True
                    odd_nkscalfact_empty = True
                else:
                    odd_nkscalfact = False
                    odd_nkscalfact_empty = False

                # Set up calculator
                calc = set_up_calculator(master_calc, calc_type,
                                         odd_nkscalfact=odd_nkscalfact,
                                         odd_nkscalfact_empty=odd_nkscalfact_empty,
                                         empty_states_nbnd=n_empty_bands,
                                         fixed_band=min(
                                             fixed_band, n_filled_bands + 1),
                                         index_empty_to_save=index_empty_to_save)
                calc.directory = directory

                # Ensure we don't overwrite KI results
                if 'ki' in calc_type:
                    calc.ndw += i_sc - 1
                    calc.setattr_only('prefix', calc.prefix + f'_it{i_sc}')

                # Run cp.x
                prev_calc_not_skipped = run_cp(
                    calc, silent=False, from_scratch=from_scratch)

                if not prev_calc_not_skipped:
                    continue

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
                    if calc.fixed_band == n_filled_bands and calc.f_cutoff == 1.00:
                        ki_calc_for_final_restart = calc
                elif 'pbe' in calc_type and 'dummy' not in calc_type:
                    if workflow_type == 'ki':
                        # For KI, the PBE calculations are independent of alpha so
                        # we store these in a list that is never overwritten
                        alpha_indep_calcs.append(calc)
                    else:
                        # For KIPZ, the PBE calculations are dependent on alpha via
                        # the definition of the variational orbitals. We only want to
                        # store the calculations that used the most recent value of alpha
                        alpha_dep_calcs.append(calc)

                # Copying of evcfixed_empty.dat to evc_occupied.dat
                prefix = calc.parameters['input_data']['control']['prefix']
                if calc_type in ['pz_print', 'kipz_print']:
                    evcempty_dir = f'calc_alpha/orbital_{fixed_band}/{calc.outdir}/' \
                                   f'{prefix}_{calc.ndw}.save/K00001/'
                elif calc_type == 'pbe_n+1_dummy':
                    evcocc_dir = f'calc_alpha/orbital_{fixed_band}/{calc.outdir}/' \
                                 f'{prefix}_{calc.ndr}.save/K00001/'
                    if os.path.isfile(f'{evcempty_dir}/evcfixed_empty.dat'):
                        os.system(
                            f'cp {evcempty_dir}/evcfixed_empty.dat {evcocc_dir}/evc_occupied.dat')
                    else:
                        raise OSError(
                            f'Could not find {evcempty_dir}/evcfixed_empty.dat')

            if not prev_calc_not_skipped:
                continue

            # Calculate an updated alpha and a measure of the error
            # E(N) - E_i(N - 1) - lambda^alpha_ii(1)     (filled)
            # E_i(N + 1) - E(N) - lambda^alpha_ii(0)     (empty)
            calcs = [c for calc_set in [alpha_dep_calcs, alpha_indep_calcs]
                     for c in calc_set if c.fixed_band == fixed_band]

            if workflow_type == 'ki':
                alpha, error = calculate_alpha(
                    calcs, filled=filled, kipz=False)
            else:
                alpha, error = calculate_alpha(
                    calcs, filled=filled, kipz=True)

            alpha_df.loc[i_sc + 1, fixed_band] = alpha
            error_df.loc[i_sc, fixed_band] = error

        if prev_calc_not_skipped:
            print_summary(alpha_df, error_df)

            # Writing alphas and errors to python-readable files and generating a .tex table
            alpha_df.to_pickle('alphas.pkl')
            error_df.to_pickle('errors.pkl')
            latex_table = alpha_df.to_latex(column_format='l' + 'd'*n_bands,
                                            float_format='{:.3f}'.format, escape=False)
            with open('tab_alpha_values.tex', 'w') as f:
                f.write(latex_table)

        converged = all([abs(e) < 1e-3 for e in error_df.loc[i_sc, :]])

    if converged:
        print('Alpha values have been converged')
    else:
        print('Alpha values have been determined but are not necessarily converged')

    # Final calculation
    print('\nFINAL KI CALCULATION')

    directory = 'final'
    calc.directory = directory
    if not os.path.isdir(directory):
        os.system(f'mkdir {directory}')

    write_alpharef(alpha_df.loc[i_sc + 1], band_filling, directory)

    if prev_calc_not_skipped:
        ndr = ki_calc_for_final_restart.ndw
        ndw = ki_calc_for_final_restart.ndw + 1
    else:
        ndr = None
        ndw = None

    calc = set_up_calculator(master_calc, workflow_type, odd_nkscalfact=True,
                             odd_nkscalfact_empty=True, empty_states_nbnd=n_empty_bands, ndw=ndw, ndr=ndr)
    calc.directory = directory

    outdir = f'{directory}/{calc.outdir}'
    if prev_calc_not_skipped:
        if not os.path.isdir(outdir):
            os.system(f'mkdir {outdir}')
        prefix = calc.parameters['input_data']['control']['prefix']
        savedir = f'{directory}/{calc.outdir}/{prefix}_{ndr}.save'
        if not os.path.isdir(savedir):
            os.system(f'cp -r calc_alpha/orbital_{n_filled_bands}/{ki_calc_for_final_restart.outdir}'
                      f'/*{ki_calc_for_final_restart.ndw}.save {savedir}')

    run_cp(calc, silent=False)

    print('\nWORKFLOW COMPLETE')
