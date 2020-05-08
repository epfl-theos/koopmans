"""

Workflow module for python_KI, containing the workflow for performing KI and KIPZ calculations

Written by Edward Linscott Jan 2020

"""

import os
import copy
import pandas as pd
from collections import namedtuple
from ase.io import espresso_cp as cp_io
from koopmans_cp.io import write_alpharef, read_alpharef, print_summary, read_json, \
    parse_algebraic_expressions, warn, print_qc
from koopmans_cp.defaults import load_defaults
from koopmans_cp.calc import run_cp, calculate_alpha, Extended_Espresso_cp

Setting = namedtuple(
    'Setting', ['name', 'description', 'type', 'default', 'options'])

valid_settings = [
    Setting('calc_type',
            'calculation type',
            str, 'ki', ('ki', 'kipz', 'pbe', 'both')),
    Setting('init_density',
            'how to initialise the density',
            str, 'pbe', ('pbe', 'pbe-pw', 'pz', 'ki')),
    Setting('init_manifold',
            'how to initialise the variational orbitals',
            str, 'pz', ('pz', 'ki', 'kipz', 'mwlf')),
    Setting('n_max_sc_steps',
            'maximum number of self-consistency steps for calculating alpha',
            int, 1, None),
    Setting('alpha_guess',
            'starting guess for alpha (overridden if alpha_from_file is true)',
            float, 0.6, None),
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
            bool, False, (True, False))]

valid_settings_dict = {s.name: s for s in valid_settings}


def check_settings(settings):
    '''
    Checks workflow settings against the list of valid settings, populates 
    missing keywords with their default values, and lowers any uppercases

    Arguments:
        settings -- a dictionary of workflow settings

    '''

    for key, value in settings.items():
        # Check key is a valid keyword
        if key in valid_settings_dict:
            valid_setting = valid_settings_dict[key]

            # Lowers any uppercase strings
            if isinstance(value, str):
                settings[key] = value.lower()
                value = value.lower()

            # Check value is the correct type
            if not isinstance(value, valid_setting.type):
                raise ValueError(
                    f'{type(value).__name__} is an invalid type for "{key}" (must be {valid_setting.type.__name__})')

            # Check value is among the valid options
            if valid_setting.options is not None and value not in valid_setting.options:
                raise ValueError(
                    f'"{value}" is an invalid value for "{key}" (options are {"/".join(valid_setting.options)})')
        else:
            raise ValueError(f'"{key}" is not a recognised workflow setting')

    # Populate missing settings with the default options
    for setting in valid_settings:
        if setting.name not in settings:
            settings[setting.name] = setting.default


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
            'pbe_init'            PBE calculation from scratch
            'pz_init'             PZ calculation starting from PBE restart
            'pz_innerloop_init'   PZ calculation starting from PBE restart (innerloop only)
            'kipz_init'           KIPZ starting from PBE restart

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

            Final calculation
            'ki_final'   Final KI calculation with N electrons and empty bands if specified
            'kipz_final' As above, but for KIPZ

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
    elif calc_type == 'pz_init':
        ndr = 50
        ndw = 51
    elif calc_type in ['pz_innerloop_init', 'kipz_init']:
        ndr = 50
        ndw = 52
    elif calc_type == 'pbe':
        ndr = 52
        ndw = 53
    elif calc_type == 'pbe_n-1':
        ndr = 52
        ndw = 54
    elif calc_type in ['pz_print', 'kipz_print']:
        ndr = 52
        ndw = 55
    elif calc_type == 'pbe_n+1_dummy':
        ndr = 56
        ndw = 56
    elif calc_type == 'pbe_n+1-1':
        ndr = 56
        ndw = 57
    elif calc_type == 'pbe_n+1':
        ndr = 56
        ndw = 58
    elif calc_type in ['ki', 'kipz', 'ki_final', 'kipz_final']:
        ndr = 52
        ndw = 60
    elif calc_type == 'kipz_n-1':
        ndr = 52
        ndw = 70
    elif calc_type == 'kipz_n+1':
        ndr = 56
        ndw = 80
    elif calc_type in ['ki_n+1-1', 'kipz_n+1-1']:
        ndr = 56
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
    if calc.prefix in ['pbe_init', 'pz_print', 'pbe_n+1_dummy', 'pz_init',
                       'kipz_init', 'kipz_print', 'ki_final', 'kipz_final',
                       'pz_innerloop_init']:
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
    if 'print' in calc.prefix or ('pz' in calc.prefix and 'kipz' not in calc.prefix and calc.prefix != 'pz_init'):
        calc.maxiter = 2
        calc.empty_states_maxstep = 1
    else:
        calc.maxiter = 300
        calc.empty_states_maxstep = 300
    # For all calculations calculating alpha, remove the empty states and
    # increase the energy thresholds
    if not any([s in calc.prefix for s in ['init', 'print', 'final']]):
        calc.empty_states_nbnd = 0
        calc.conv_thr *= 100
        calc.esic_conv_thr *= 100

    # nksic
    calc.odd_nkscalfact = calc.do_orbdep
    calc.odd_nkscalfact_empty = calc.do_orbdep
    calc.do_innerloop_cg = True
    if calc.prefix[:2] == 'pz' and calc.prefix != 'pz_init':
        calc.do_innerloop = True
        calc.one_innerloop_only = True
    elif 'kipz' in calc.prefix or calc.prefix == 'pz_init':
        calc.do_innerloop = True
        calc.one_innerloop_only = False
    else:
        calc.do_innerloop = False
        calc.one_innerloop_only = False
    if 'kipz' in calc.prefix:
        calc.which_orbdep = 'nkipz'
    elif 'pz' in calc.prefix:
        calc.which_orbdep = 'pz'
    elif 'ki' in calc.prefix:
        calc.which_orbdep = 'nki'
    if 'print' in calc.prefix:
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
        warn('calc.fixed_band is higher than the LUMO; this should not happen')

    return calc


keywords_altered_during_workflow = ['ndw', 'ndr', 'restart_mode', 'nspin', 'nelec', 'nelup',
                                    'neldw', 'do_orbdep', 'fixed_state', 'f_cutoff',
                                    'fixed_band', 'conv_thr', 'restart_from_wannier_pwscf',
                                    'maxiter', 'empty_states_maxstep', 'esic_conv_thr',
                                    'do_innerloop', 'one_innerloop_only', 'which_orbdep',
                                    'print_wfc_anion', 'directory', 'prefix', 'odd_nkscalfact',
                                    'odd_nkscalfact_empty', 'nkscalfact']


def run_from_json(json):
    '''
    This function wraps 'run' to allow reading of a json file rather than
    requiring an ASE calculator object and a workflow settings dictionary

    This is kept separate from 'run' to allow for calc_type == "both".

    Arguments:
        json -- the json file containing all the workflow and cp.x settings

    '''

    # Reading in JSON file
    master_calc, workflow_settings = read_json(json)

    # Check that the workflow settings are valid and populate missing
    # settings with default values
    check_settings(workflow_settings)

    # Convert the master_calc object to a the Extended_Espresso_cp class
    master_calc = Extended_Espresso_cp(master_calc)

    # Load cp.x default values from koopmans.defaults
    master_calc = load_defaults(master_calc)

    # Parse any algebraic expressions used for cp.x keywords
    master_calc = parse_algebraic_expressions(master_calc)

    if workflow_settings['calc_type'] == 'both':
        # if 'both', create subdirectories and run
        calc_types = ['ki', 'kipz']

        # Make separate directories for KI and KIPZ
        for calc_type in calc_types:
            if workflow_settings['from_scratch'] and os.path.isdir(calc_type):
                os.system(f'rm -r {calc_type}')
            if not os.path.isdir(calc_type):
                os.system(f'mkdir {calc_type}')

        for calc_type in calc_types:
            print(f'\n{calc_type.upper()} CALCULATION')

            workflow_settings['calc_type'] = calc_type

            # For KIPZ, use KI as a starting point
            if calc_type == 'kipz':
                workflow_settings['init_density'] = 'ki'
                workflow_settings['init_manifold'] = 'kipz'
                workflow_settings['alpha_from_file'] = True

            # Change to relevant subdirectory
            os.chdir(calc_type)

            # Run workflow for this particular functional
            run(master_calc, workflow_settings)

            # Return to the base directory
            os.chdir('..')

            # Provide the KIPZ calculation with a KI starting point
            if calc_type == 'ki':
                os.system('cp -r ki/final kipz/init')
                os.system('mv kipz/init/ki_final.cpi kipz/init/ki_init.cpi')
                os.system('mv kipz/init/ki_final.cpo kipz/init/ki_init.cpo')
                os.system('cp -r ki/final/file_alpharef* kipz/')
    else:
        # If not 'both', simply run the workflow once
        run(master_calc, workflow_settings)


def run(master_calc, workflow_settings):
    '''
    This function runs the KI/KIPZ workflow from start to finish

    Arguments:
        master_calc       -- the master ASE calculator object containing cp.x settings
        workflow_settings -- a dictionary containing workflow settings

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

    # Check that the workflow settings are valid and populate missing
    # settings with default values
    check_settings(workflow_settings)

    # Extracting workflow_type (not to be confused with the internal 'calc_type' for
    # individual calculations)
    workflow_type = workflow_settings['calc_type']

    # Removing old directories
    if workflow_settings['from_scratch']:
        if workflow_settings['init_density'] != 'ki':
            # if init_density == "ki" we don't want to delete the directory containing
            # the KI calculation we're reading the manifold from
            os.system('rm -r init 2>/dev/null')
        os.system('rm -r calc_alpha 2>/dev/null')
        os.system('rm -r final 2>/dev/null')

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
        print(r'Reading alpha values from file')
        alpha_df.loc[1] = read_alpharef(directory='.')
    else:
        alpha_df.loc[1] = [workflow_settings['alpha_guess']
                           for _ in range(n_bands)]

    prev_calc_not_skipped = workflow_settings['from_scratch']
    # Note since we always provide 'from_scratch=prev_calc_not_skipped' to run_cp, if
    # workflow_settings['from_scratch'] is False, the workflow will skip all calculations
    # until it reaches an incomplete calculation, at which stage prev_calc_not_skipped
    # will go from False to True and all subsequent calculations will be run with
    # 'from_scratch = True'.

    # On the other hand, if workflow_settings['from_scratch'] is True,
    # prev_calc_not_skipped will remain equal to True throughout the workflow and no
    # calculations will be skipped

    # If workflow_type is PBE, perform PBE and exit immediately
    if workflow_type == 'pbe':
        # PBE
        calc = set_up_calculator(master_calc, 'pbe_init')
        calc.directory = '.'
        calc.prefix = 'pbe'
        prev_calc_not_skipped = run_cp(
            calc, silent=False, from_scratch=prev_calc_not_skipped)
        return

    print('\nINITIALISATION OF DENSITY')

    init_density = workflow_settings.get('init_density', 'pbe').lower()
    if init_density == 'pbe':
        # PBE from scratch
        calc = set_up_calculator(master_calc, 'pbe_init')
        calc.directory = 'init'
        prev_calc_not_skipped = run_cp(
            calc, silent=False, from_scratch=prev_calc_not_skipped)

    elif init_density == 'pbe-pw':
        # PBE using pw.x
        raise ValueError('init_denisty: pbe-pw is not yet implemented')

    elif init_density == 'pz':
        # PBE from scratch
        calc = set_up_calculator(master_calc, 'pbe_init')
        calc.directory = 'init'
        prev_calc_not_skipped = run_cp(
            calc, silent=False, from_scratch=prev_calc_not_skipped)
        # PZ from PBE
        calc = set_up_calculator(master_calc, 'pz_init')
        calc.directory = 'init'
        write_alpharef(alpha_df.loc[1], band_filling, calc.directory)
        prev_calc_not_skipped = run_cp(
            calc, silent=False, from_scratch=prev_calc_not_skipped)

    elif init_density == 'ki':
        print('Initialising the density with a pre-existing KI calculation')
        # a pre-existing, complete KI calculation
        prefix = master_calc.parameters['input_data']['control']['prefix']

        # Read the .cpi file to work out the value for ndw
        atoms = cp_io.read_espresso_cp_in('init/ki_init.cpi')
        calc = Extended_Espresso_cp(atoms.calc)

        # Move the old save directory to correspond to ndw = 50
        old_savedir = f'init/{master_calc.outdir}/{prefix}_{calc.ndw}.save'
        savedir = f'init/{master_calc.outdir}/{prefix}_50.save'
        if not os.path.isdir(old_savedir):
            raise ValueError(f'{old_savedir} does not exist; a previous '
                             'and complete KI calculation is required '
                             'if init_density=ki')
        if os.path.isdir(savedir):
            raise ValueError(f'{savedir} should not already exist')
        os.system(f'mv {old_savedir} {savedir}')

        # Check that the files defining the variational orbitals exist
        savedir += '/K00001'
        files_to_check = ['init/ki_init.cpo',
                          f'{savedir}/evc01.dat', f'{savedir}/evc02.dat']

        if calc.empty_states_nbnd is not None and calc.empty_states_nbnd > 0:
            files_to_check += [f'{savedir}/evc0_empty{i}.dat' for i in [1, 2]]

        for fname in files_to_check:
            if not os.path.isfile(fname):
                raise ValueError(f'Could not find {fname}')

        calc = next(cp_io.read_espresso_cp_out('init/ki_init.cpo')).calc
        if not calc.results['job_done']:
            raise ValueError('init/ki_init.cpo is incomplete so cannot be used '
                             'to initialise the density')
    else:
        raise ValueError(
            "Should not arrive here; compare the above code with workflow.valid_settings")

    if init_density == 'pbe':
        # Using KS eigenfunctions as guess variational orbitals
        print('Overwriting the CP variational orbitals with Kohn-Sham orbitals')
        prefix = calc.parameters['input_data']['control']['prefix']
        savedir = f'{calc.directory}/{calc.outdir}/{prefix}_{calc.ndw}.save/K00001'
        os.system(f'cp {savedir}/evc1.dat {savedir}/evc01.dat')
        os.system(f'cp {savedir}/evc2.dat {savedir}/evc02.dat')
        if calc.empty_states_nbnd is not None and calc.empty_states_nbnd > 0:
            os.system(f'cp {savedir}/evc_empty1.dat {savedir}/evc0_empty1.dat')
            os.system(f'cp {savedir}/evc_empty2.dat {savedir}/evc0_empty2.dat')

    print('\nINITIALISATION OF MANIFOLD')

    init_manifold = workflow_settings['init_manifold']
    if init_manifold == 'pz':
        write_alpharef(alpha_df.loc[1], band_filling, calc.directory)
        calc = set_up_calculator(master_calc, 'pz_innerloop_init')
    elif init_manifold == 'kipz':
        write_alpharef(alpha_df.loc[1], band_filling, calc.directory)
        calc = set_up_calculator(master_calc, 'kipz_init')
    elif init_manifold == 'mlwf':
        raise ValueError('mlwf initialisation not yet implemented')
    else:
        raise ValueError(
            f'Unrecognised option "{init_manifold}" for init_manifold. Should be one of "pz"/"kipz"/"mlwf"')

    calc.directory = 'init'
    write_alpharef(alpha_df.loc[1], band_filling, calc.directory)
    prev_calc_not_skipped = run_cp(
        calc, silent=False, from_scratch=prev_calc_not_skipped)

    if prev_calc_not_skipped:
        print('Copying the spin-up variational orbitals over to the spin-down channel')
        prefix = calc.parameters['input_data']['control']['prefix']
        savedir = f'{calc.directory}/{calc.outdir}/{prefix}_{calc.ndw}.save/K00001'
        os.system(f'cp {savedir}/evc01.dat {savedir}/evc02.dat')
        if calc.empty_states_nbnd is not None and calc.empty_states_nbnd > 0:
            os.system(
                f'cp {savedir}/evc0_empty1.dat {savedir}/evc0_empty2.dat')

    print('\nDETERMINING ALPHA VALUES')

    # Set up directories
    if not os.path.isdir('calc_alpha'):
        os.system('mkdir calc_alpha')
    for i_band in i_bands:
        if not os.path.isdir(f'calc_alpha/orbital_{i_band}'):
            os.system('cp -r init calc_alpha/orbital_{}'.format(i_band))

    converged = False
    i_sc = 0

    if not prev_calc_not_skipped:
        # Reloading alphas and errors from file
        print('Reloading alpha values from file')
        alpha_df = pd.read_pickle('alphas.pkl')
        error_df = pd.read_pickle('errors.pkl')

    alpha_indep_calcs = []

    while not converged and i_sc < workflow_settings['n_max_sc_steps']:
        i_sc += 1
        alpha_dep_calcs = []

        if workflow_settings['n_max_sc_steps'] > 1:
            print('\n== SC iteration {} ==================='.format(i_sc))

        # Loop over removing an electron from each band
        for fixed_band, filled in zip(i_bands, band_filling):
            print('-- Orbital {} ------------------------'.format(fixed_band))
            directory = 'calc_alpha/orbital_{}'.format(fixed_band)

            # Don't repeat if this particular alpha_i was converged
            if i_sc > 1 and any([abs(e) < 1e-3 for e in 
                                 error_df.loc[:i_sc - 1,fixed_band]]):
                print(
                    f'Skipping band {fixed_band} since this alpha is already '
                    'converged')
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
                    elif not filled:
                        # In the case of empty orbitals, we gain an extra orbital in
                        # the spin-up channel
                        alpha_padded = list(alpha_df.loc[i_sc])
                        alpha_padded += [alpha_padded[-1]] + alpha_padded
                        write_alpharef(alpha_padded, band_filled_or_fixed + [False]
                                       + band_filling, directory, duplicate=False)
                    else:
                        write_alpharef(alpha_df.loc[i_sc],
                                       band_filled_or_fixed, directory)

                if filled:
                    index_empty_to_save = None
                else:
                    index_empty_to_save = fixed_band - n_filled_bands

                # Set up calculator
                calc = set_up_calculator(master_calc, calc_type,
                                         fixed_band=min(
                                             fixed_band, n_filled_bands + 1),
                                         index_empty_to_save=index_empty_to_save)
                calc.directory = directory

                # Ensure we don't overwrite varational-orbital-dependent results
                if workflow_type == 'kipz' or 'ki' in calc_type:
                    calc.setattr_only('prefix', calc.prefix + f'_it{i_sc}')

                # Run cp.x
                prev_calc_not_skipped = run_cp(
                    calc, silent=False, from_scratch=prev_calc_not_skipped)

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

            if prev_calc_not_skipped:

                # Calculate an updated alpha and a measure of the error
                # E(N) - E_i(N - 1) - lambda^alpha_ii(1)     (filled)
                # E_i(N + 1) - E(N) - lambda^alpha_ii(0)     (empty)
                #
                # Note that we only do this if the calculation has not been skipped
                # because calculate_alpha reads in alpha values from files which get
                # overwritten by subsequent calculations

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
                alpha_df.to_pickle('alphas.pkl')
                error_df.to_pickle('errors.pkl')

        print_summary(alpha_df, error_df)

        converged = all([abs(e) < 1e-3 for e in error_df.loc[i_sc, :]])

    if converged:
        print('Alpha values have been converged')
    else:
        print('Alpha values have been determined but are not necessarily converged')

    # Writing alphas to a .tex table
    latex_table = alpha_df.to_latex(column_format='l' + 'd'*n_bands,
                                    float_format='{:.3f}'.format, escape=False)
    with open('tab_alpha_values.tex', 'w') as f:
        f.write(latex_table)

    # Final calculation
    print('\nFINAL KI CALCULATION')

    directory = 'final'
    calc.directory = directory
    if not os.path.isdir(directory):
        os.system(f'mkdir {directory}')

    write_alpharef(alpha_df.loc[i_sc + 1], band_filling, directory)

    calc = set_up_calculator(master_calc, workflow_type + '_final',
                             empty_states_nbnd=n_empty_bands)
    calc.directory = directory

    # Read in from (note that if we have empty states only init/ and final/ use the full complement of orbitals)
    outdir = f'{directory}/{calc.outdir}'
    if prev_calc_not_skipped:
        if not os.path.isdir(outdir):
            os.system(f'mkdir {outdir}')
        prefix = calc.parameters['input_data']['control']['prefix']
        savedir = f'{calc.outdir}/{prefix}_{calc.ndr}.save'
        if not os.path.isdir(savedir):
            os.system(f'cp -r init/{savedir} final/{savedir}')

    run_cp(calc, silent=False, from_scratch=prev_calc_not_skipped)

    # Print out data for quality control
    if workflow_settings.get('print_qc', False):
        print('\nQUALITY CONTROL')
        print_qc('energy', calc.results['energy'])
        for i, alpha in enumerate(alpha_df.loc[i_sc + 1]):
            print_qc(f'alpha({i})', alpha)
        for isp, orbs_self_hartree in enumerate(calc.results['orbital_data']['self-Hartree']):
            for i, orb_sh in enumerate(orbs_self_hartree):
                print_qc(f'orbital_self_Hartree(orb={i},sigma={isp})', orb_sh)
        print_qc('HOMO', calc.results['homo_energy'])
        if calc.results['lumo_energy'] is not None:
            print_qc('LUMO', calc.results['lumo_energy'])

    print('\nWORKFLOW COMPLETE')
