import argparse
import os
import subprocess
import pickle
import pandas as pd
from ase.io import espresso_cp as cp_io
from koopmans_utils import run_cp, calculate_alpha, set_up_calculator, \
    Extended_Espresso_cp, write_alpharef, read_alpharef

'''
Perform KI calculations
'''


def print_summary(alpha_df, error_df):
    # Printing out a progress summary
    print('\nalpha')
    print(alpha_df)
    print('\ndelta E - lambda^alpha_ii (eV)')
    print(error_df)
    print('')


def run_ki(master_cpi, alpha_guess=0.6, alpha_from_file=False, n_max_sc_steps=1, from_scratch=False):
    '''
    This function runs the KI workflow from start to finish

    Arguments:
         master_cpi     -- the path to the master cp.x input file from which to take
                         all settings from
        alpha_guess     -- the initial guess for alpha (a single float for all orbitals)
        alpha_from_file -- read alpha from pre-existing file_alpharef.txt
       n_max_sc_steps -- the maximum number of self-consistent steps for 
                         determining {alpha_i}

    Running this function will generate a number of files:
       init/ -- the PBE and PZ calculations for initialisation
       fix_orbital_#/ -- PBE/PZ/KI calculations where we have fixed a particular orbital
       final/ -- the directory containing the final KI calculation
       alphas.pkl -- a python pickle file containing the alpha values
       errors.pkl -- a python pickle file containing the errors in the alpha values
       tab_alpha_values.tex -- a latex table of the alpha values
    '''

    # Removing old directories
    if from_scratch:
        os.system('rm -r init 2>/dev/null')
        os.system('rm -r fix_orbital* 2>/dev/null')
        os.system('rm -r final 2>/dev/null')

    # Reading in template cp input file
    master_atoms = cp_io.read_espresso_cp_in(open(master_cpi, 'r'))
    master_calc = Extended_Espresso_cp(master_atoms.calc)

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

    print('\nINITIALISATION OF DENSITY AND VARIATIONAL ORBITALS')
    # PBE from scratch
    calc = set_up_calculator(master_calc, 'pbe_init')
    calc.directory = 'init'
    run_cp(calc, silent=False, from_scratch=from_scratch)

    # # Moving orbitals
    # os.system('cp TMP-CP/pbe_50.save/K00001/evc1.dat TMP-CP/pbe_50.save/K00001/evc01.dat')
    # os.system('cp TMP-CP/pbe_50.save/K00001/evc2.dat TMP-CP/pbe_50.save/K00001/evc02.dat')

    # PZ reading in PBE to define manifold
    if alpha_from_file:
        print(r'Reading alpha values from file')
        alpha_df.loc[1] = read_alpharef(directory='.')
    else:
        alpha_df.loc[1] = [alpha_guess for _ in range(n_bands)]

    calc = set_up_calculator(
        master_calc, 'pz', empty_states_nbnd=n_empty_bands, from_scratch=from_scratch)
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
    for i_band in i_bands:
        if not os.path.isdir(f'fix_orbital_{i_band}'):
            os.system('cp -r init fix_orbital_{}'.format(i_band))

    converged = False
    i_sc = 0
    pbe_calcs = []

    if not prev_calc_not_skipped:
        # Reloading from file
        alpha_df = pd.read_pickle('alphas.pkl')
        error_df = pd.read_pickle('errors.pkl')
        i_sc = len(error_df)
        print_summary(alpha_df, error_df)
        converged = all([abs(e) < 1e-3 for e in error_df.loc[i_sc, :]])

    while not converged and i_sc < n_max_sc_steps:
        i_sc += 1

        if n_max_sc_steps > 1:
            print('\n== SC iteration {} ==================='.format(i_sc))

        # Loop over removing an electron from each band
        for fixed_band, filled in zip(i_bands, band_filling):
            print('-- Orbital {} ------------------------'.format(fixed_band))
            directory = 'fix_orbital_{}'.format(fixed_band)

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

            if prev_calc_not_skipped:
                write_alpharef(alpha_df.loc[i_sc],
                               band_filled_or_fixed, directory)

            # Perform the fixed-band-dependent calculations: PBE, PBE_N-1/PBE_N+1, and KI
            if filled:
                calc_types = ['pbe', 'pbe_n-1', 'ki']
            else:
                calc_types = ['pz_print', 'pbe_n+1_dummy', 'pbe_n+1',
                              'pbe_n+1-1', 'ki_n+1-1']

            for calc_type in calc_types:
                # Only the KI calculations change with alpha; we don't need to
                # redo any of the others
                if i_sc > 1 and 'ki' not in calc_type:
                    continue

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
                if 'ki' in calc_type:
                    ki_calc = calc
                    if calc.fixed_band == n_filled_bands and calc.f_cutoff == 1.00:
                        ki_calc_for_final_restart = calc
                elif 'pbe' in calc_type and 'dummy' not in calc_type:
                    pbe_calcs.append(calc)

                # Copying of evcfixed_empty.dat to evc_occupied.dat
                prefix = calc.parameters['input_data']['control']['prefix']
                if calc_type == 'pz_print':
                    evcempty_dir = f'fix_orbital_{fixed_band}/{calc.outdir}/' \
                                   f'{prefix}_{calc.ndw}.save/K00001/'
                elif calc_type == 'pbe_n+1_dummy':
                    evcocc_dir = f'fix_orbital_{fixed_band}/{calc.outdir}/' \
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
            calcs = [c for c in pbe_calcs if c.fixed_band ==
                     fixed_band] + [ki_calc]
            alpha, error = calculate_alpha(calcs, filled=filled)
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

    calc = set_up_calculator(master_calc, 'ki', odd_nkscalfact=True,
                             odd_nkscalfact_empty=True, empty_states_nbnd=n_empty_bands, ndw=ndw, ndr=ndr)
    calc.directory = directory

    outdir = f'{directory}/{calc.outdir}'
    if prev_calc_not_skipped:
        if not os.path.isdir(outdir):
            os.system(f'mkdir {outdir}')
        prefix = calc.parameters['input_data']['control']['prefix']
        savedir = f'{directory}/{calc.outdir}/{prefix}_{ndr}.save'
        if not os.path.isdir(savedir):
            os.system(f'cp -r fix_orbital_{n_filled_bands}/{ki_calc_for_final_restart.outdir}'
                      f'/*{ki_calc_for_final_restart.ndw}.save {savedir}')

    run_cp(calc, silent=False)

    print('\nWORKFLOW COMPLETE')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Perform a KI calculation using cp.x')
    parser.add_argument('template', metavar='template.cpi', type=str,
                        help='a template cp input file')
    parser.add_argument('-a', '--alpha', default=0.6, type=float,
                        help='starting guess for alpha as a single float')
    parser.add_argument('-f', '--alpha_from_file', action='store_true',
                        help='read in starting guess for alpha from file')
    parser.add_argument('-i', '--maxit', default=1, type=int,
                        help='maximum number of self-consistent iterations')
    parser.add_argument('-c', '--cont', action='store_true',
                        help='continue from a previous calculation')

    args = parser.parse_args()

    run_ki(args.template, args.alpha,
           args.alpha_from_file, args.maxit, not args.cont)
