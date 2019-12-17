#!/usr/bin/python3

import argparse
import os
import subprocess
import pickle
import pandas as pd
from ase.io import espresso_cp as cp_io
from koopmans_utils import run_cp, calculate_alpha, set_up_calculator, \
    Extended_Espresso_cp, write_alpharef


'''
Perform KI calculations
'''


def run_ki(master_cpi, alpha_guess=0.6, n_max_sc_steps=1):
    '''
    This function runs the KI workflow from start to finish

    Arguments:
       master_cpi     -- the path to the master cp.x input file from which to take
                         all settings from
       alpha_guess    -- the initial guess for alpha
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

    print('\nINITIALISATION OF DENSITY AND VARIATIONAL ORBITALS')
    # PBE from scratch
    calc = set_up_calculator(master_calc, 'pbe_init',
                             empty_states_nbnd=n_empty_bands)
    calc.directory = 'init'
    run_cp(calc, silent=False)

    # PZ reading in PBE to define manifold
    calc = set_up_calculator(
        master_calc, 'pz', empty_states_nbnd=n_empty_bands)
    calc.directory = 'init'
    run_cp(calc, silent=False)

    print('\nDETERMINING ALPHA VALUES')
    # Preparing panda dataframes in which to store results
    alpha_df = pd.DataFrame(columns=i_bands)
    error_df = pd.DataFrame(columns=i_bands)

    # Set up directories
    for i_band in i_bands:
        os.system('cp -r init fix_orbital_{}'.format(i_band))

    converged = False
    i_sc = 0
    alpha_df.loc[1] = [alpha_guess for _ in range(n_bands)]
    pbe_calcs = []

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
            write_alpharef(alpha_df.loc[i_sc], band_filled_or_fixed, directory)

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
                run_cp(calc, silent=False)

                # Reset the value of 'fixed_band' so we can keep track of which calculation
                # is which. This is important for empty orbital calculations, where fixed_band
                # is always set to the LUMO but in reality we're fixing the band corresponding
                # to index_empty_to_save from an earlier calculation
                calc.fixed_band = fixed_band

                # Store the result
                if 'ki' in calc_type:
                    ki_calc = calc
                    if calc.fixed_band == n_filled_bands:
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

            # Calculate an updated alpha and a measure of the error
            # E(N) - E_i(N - 1) - lambda^alpha_ii(1)     (filled)
            # E_i(N + 1) - E(N) - lambda^alpha_ii(0)     (empty)
            calcs = [c for c in pbe_calcs if c.fixed_band ==
                     fixed_band] + [ki_calc]
            alpha, error = calculate_alpha(calcs, filled=filled)
            alpha_df.loc[i_sc + 1, fixed_band] = alpha
            error_df.loc[i_sc, fixed_band] = error

        # Printing out a progress summary
        print('\nSummary\nalpha')
        print(alpha_df)
        print('\ndelta E - lambda^alpha_ii (eV)')
        print(error_df)
        print('')
        converged = all([abs(e) < 1e-3 for e in error_df.loc[i_sc, :]])

    if converged:
        print('Alpha values have been converged')
    else:
        print('Alpha values have been determined but are not necessarily converged')

    # Writing alphas and errors to python-readable files and generating a .tex table
    alpha_df.to_pickle('alphas.pkl')
    error_df.to_pickle('errors.pkl')
    latex_table = alpha_df.to_latex(column_format='l' + 'd'*n_bands,
                                    float_format='{:.3f}'.format, escape=False)
    with open('tab_alpha_values.tex', 'w') as f:
        f.write(latex_table)

    # Final calculation
    print('\nFINAL KI CALCULATION')

    directory = 'final'
    calc.directory = directory
    os.system(f'mkdir {directory}')

    write_alpharef(alpha_df.loc[i_sc + 1], band_filling, directory)
    calc = set_up_calculator(master_calc, 'ki', odd_nkscalfact=True,
                             odd_nkscalfact_empty=True, empty_states_nbnd=n_empty_bands,
                             ndr=ki_calc_for_final_restart.ndw,
                             ndw=ki_calc_for_final_restart.ndw + 1)
    calc.directory = directory

    os.system(f'mkdir {directory}/{calc.outdir}')
    os.system(f'cp -r fix_orbital_{n_filled_bands}/{ki_calc_for_final_restart.outdir}'
              f'/*{ki_calc_for_final_restart.ndw}.save {directory}/{calc.outdir}/')
    run_cp(calc, silent=False)

    print('\nWORKFLOW COMPLETE')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Perform a KI calculation using cp.x')
    parser.add_argument('template', metavar='template.cpi', type=str,
                        help='a template cp input file')
    parser.add_argument('-a', '--alpha', default=0.6, type=float,
                        help='Starting guess for alpha')
    parser.add_argument('-i', '--maxit', default=1, type=int,
                        help='Maximum number of self-consistent iterations')

    args = parser.parse_args()

    run_ki(args.template, args.alpha, args.maxit)
