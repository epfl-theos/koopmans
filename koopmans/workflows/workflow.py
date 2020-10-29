"""

Generic workflow functions for python_KI

Written by Edward Linscott Jan 2020
KI/pKIPZ/KIPZ moved to kc_with_cp.py Oct 2020

"""

import os
import copy
import pandas as pd
import numpy as np
import itertools
from collections import namedtuple
from koopmans import utils, io
from koopmans.ase import read_json, write_json
from koopmans.defaults import load_defaults
from koopmans.workflows import kc_with_cp, pbe_with_cp

Setting = namedtuple(
    'Setting', ['name', 'description', 'type', 'default', 'options'])

valid_settings = [
    Setting('task',
            'Task to perform',
            str, 'singlepoint', ('singlepoint', 'convergence')),
    Setting('functional',
            'Orbital-density-dependent-functional/density-functional to use',
            str, 'ki', ('ki', 'kipz', 'pkipz', 'pbe', 'all')),
    Setting('init_density',
            'how to initialise the density',
            str, 'pbe', ('pbe', 'pbe-pw', 'pz', 'ki')),
    Setting('init_manifold',
            'how to initialise the variational orbitals',
            str, 'pz', ('pz', 'ki', 'mwlf')),
    Setting('n_max_sc_steps',
            'maximum number of self-consistency steps for calculating alpha',
            int, 1, None),
    Setting('alpha_conv_thr',
            'convergence threshold for |delta E - lambda|; if below this '
            'threshold, the corresponding alpha value is not updated',
            (float, str), 1e-3, None),
    Setting('calculate_alpha',
            'if True, the screening parameters will be calculated; if False, '
            'they will be read directly from file',
            bool, True, (True, False)),
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
            bool, False, (True, False)),
    Setting('orbital_groups',
            'a list of integers the same length as the total number of bands, '
            'denoting which bands to assign the same screening parameter to',
            list, None, None),
    Setting('enforce_spin_symmetry',
            'if True, the spin-up and spin-down wavefunctions will be forced '
            'to be the same',
            bool, True, (True, False)),
    Setting('convergence_observable',
            'System observable of interest which we converge',
            str, 'total energy', ('homo energy', 'lumo energy', 'total energy')),
    Setting('convergence_threshold',
            'Convergence threshold for the system observable of interest',
            str, None, None),
    Setting('convergence_parameters',
            'The observable of interest will be converged with respect to this/these'
            'simulation parameter(s)',
            (list, str), ['ecutwfc'], None)]

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
            if not isinstance(value, valid_setting.type) and not value is None:
                if isinstance(valid_setting.type, tuple):
                    raise ValueError(
                        f'{type(value).__name__} is an invalid type for "{key}" (must be '
                        'one of ' + '/'.join([t.__name__ for t in valid_setting.type]) + ')')
                else:
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

    # Parse physicals
    for physical in ['alpha_conv_thr', 'convergence_threshold']:
        settings[physical] = io.parse_physical(settings[physical])

def run_from_json(json):
    '''
    This function reads the input json file, checks its contents, and
    then runs a workflow

    Arguments:
        json -- the json file containing all the workflow and cp.x settings

    '''

    # Reading in JSON file
    workflow_settings, calcs_dct = read_json(json)
    cp_master_calc = calcs_dct['cp']

    # Check that the workflow settings are valid and populate missing
    # settings with default values
    check_settings(workflow_settings)

    # Load cp.x default values from koopmans.defaults
    cp_master_calc = load_defaults(cp_master_calc)

    # Parse any algebraic expressions used for cp.x keywords
    cp_master_calc.parse_algebraic_settings()

    # Set up pseudopotentials, by...
    #  1. trying to locating the directory as currently specified by the calculator
    #  2. if that fails, checking if $ESPRESSO_PSEUDO is set
    #  3. if that fails, raising an OS error
    if cp_master_calc.pseudo_dir is None or not os.path.isdir(cp_master_calc.pseudo_dir):
        try:
            cp_master_calc.pseudo_dir = os.environ.get('ESPRESSO_PSEUDO')
        except:
            raise NotADirectoryError('Directory for pseudopotentials not found. Please define '
                          'the environment variable ESPRESSO_PSEUDO or provide a pseudo_dir in '
                          'the cp block of your json input file.')

    if workflow_settings['task'] == 'convergence':
        run_convergence(workflow_settings, cp_master_calc)
    elif workflow_settings['task'] == 'singlepoint':
        run_singlepoint(workflow_settings, cp_master_calc)

    return

def run_convergence(workflow_settings, cp_master_calc, initial_depth=3):

    increments = {'cell_size': 0.1, 'ecutwfc': 10, 'empty_states_nbnd': 1}

    if workflow_settings['from_scratch']:
        for key in increments.keys():
            utils.system_call(f'rm -r {key}* 2> /dev/null', False)

    # Populate the param_dict
    param_dict = {}
        
    # Convert to list
    convergence_parameters = workflow_settings['convergence_parameters']
    if isinstance(convergence_parameters, str):
        convergence_parameters = [convergence_parameters]

    for parameter in convergence_parameters:
        if parameter not in increments:
            raise NotImplementedError(f'Convergence wrt {parameter} has not yet been implemented')

        if parameter == 'cell_size':
            param_dict[parameter] = [1 + i*increments['cell_size'] for i in range(initial_depth)]
        else:
            attr = getattr(cp_master_calc, parameter, None)
            if attr is None:
                raise AttributeError(f'In order to converge wrt {parameter} specify a baseline value for it')
            param_dict[parameter] = [attr + i*increments[parameter] for i in range(initial_depth)]

    # Create array for storing calculation results
    results = np.empty([initial_depth for _ in param_dict])
    results[:] = np.nan

    # Continue to increment the convergence parameters until convergence in the convergence observable 
    # is achieved. A lot of the following code is quite obtuse, in order to make it work for an arbitrary
    # number of convergence parameters

    while True:

        # Loop over all possible permutations of convergence parameter settings
        for indices in itertools.product(*[range(len(x)) for x in param_dict.values()]):
            # If we already have results for this permutation, don't recalculate it
            if not np.isnan(results[tuple(indices)]):
                continue

            # Create duplicate calculator
            cp_calc = copy.deepcopy(cp_master_calc)

            # For each parameter we're converging wrt...
            header = ''
            for index, param, values in zip(indices, param_dict.keys(), param_dict.values()):
                value = values[index]
                if isinstance(value, int):
                    value_str = str(value)
                elif isinstance(value, float):
                    value_str = '{:.1f}'.format(value)
                header += f'{param} = {value_str}, '
                    
                # Create new working directory
                subdir = f'{param}_{value_str}'.replace(' ', '_').replace('.', 'd')
                cp_calc.directory += '/' + subdir

                if param == 'cell_size':
                    cp_calc._ase_calc.atoms.cell *= value
                else:
                    setattr(cp_calc, param, value)
                if param == 'ecutwfc':
                    setattr(cp_calc, 'ecutrho', 4*value)

            print(header.rstrip(', '))

            # Perform calculation
            solved_calc = run_singlepoint(workflow_settings, cp_calc)

            # Store the result
            obs = workflow_settings['convergence_observable']
            obs = obs.replace('total ', '').replace(' ', '_')
            if obs not in solved_calc.results:
                raise ValueError(f'{solved_calc.name} has not returned a value for {obs}')
            result = solved_calc.results[obs]
            results[indices] = result

            print()

        # Check convergence
        converged = np.empty([len(x) for x in param_dict.values()], dtype=bool)
        converged[:] = False

        converged_result = results[tuple([-1 for _ in param_dict])]
        threshold = workflow_settings['convergence_threshold']
        for indices in itertools.product(*[range(len(x)) for x in param_dict.values()]):
            subarray_slice = [slice(None) for _ in range(len(indices))]
            for i_index, index in enumerate(indices):
                subarray_slice[i_index] = slice(index, None)
            if np.all(np.abs(results[tuple(subarray_slice)] - converged_result) < threshold):
                converged[tuple(indices)] = True

        # Check convergence, ignoring the calculations with the most fine parameters because
        # we have nothing to compare them against
        subarray_slice = tuple([slice(0, -1) for _ in param_dict])
        if np.any(converged[subarray_slice]):
            # Work out which was the fastest calculation, and propose those parameters

            # First, find the indices of the converged array
            slice_where = [slice(None) for _ in param_dict]
            slice_where[-1] = 0
            indices = np.array(np.where(converged[subarray_slice]))[tuple(slice_where)]

            # Extract the corresponding parameters
            converged_parameters = {}
            for index, param in zip(indices, param_dict.keys()):
                converged_parameters[param] = param_dict[param][index]

            print('Converged parameters are ' + ', '.join([f'{k} = {v}' for k, v in converged_parameters.items()]))

            # Construct a calculator with the converged settings
            # Abuse the fact we don't routinely update the ase settings to ignore the various defaults that
            # have been set
            cp_master_calc._update_settings_dict()
            for param, value in converged_parameters.items():
                if param == 'cell_size':
                    cp_master_calc._ase_calc.atoms.cell *= value
                else:
                    setattr(cp_master_calc, param, value)
                if param == 'ecutwfc':
                    setattr(cp_master_calc, 'ecutrho', 4*value)

                if workflow_settings['print_qc']:
                    io.print_qc(param, value)
            cp_master_calc.ibrav = 0

            # Save converged settings to a .json file
            with open('converged.json', 'w') as fd:
                write_json(fd, cp_master_calc, {})

            return converged_parameters
        else:
            # Work out which parameters are yet to converge, and line up more calculations
            # for increased values of those parameters

            new_array_shape = list(np.shape(results))
            new_array_slice = [slice(None) for _ in indices]
            for index, param in enumerate(param_dict):
                subarray_slice = [-1 for _ in indices]
                subarray_slice[index] = slice(None)
                if np.all(converged[tuple(subarray_slice)]):
                    continue
                param_dict[param].append(param_dict[param][-1] + increments[param])
                new_array_shape[index] += 1
                new_array_slice[index] = slice(None, -1)

            new_results = np.empty(new_array_shape)
            new_results[:] = np.nan
            new_results[tuple(new_array_slice)] = results
            results = new_results

def run_singlepoint(workflow_settings, cp_master_calc):

    if workflow_settings['functional'] == 'all':
        # if 'all', create subdirectories and run
        functionals = ['ki', 'pkipz', 'kipz']

        # Make separate directories for KI, pKIPZ, and KIPZ
        for functional in functionals:
            if workflow_settings['from_scratch'] and os.path.isdir(functional):
                utils.system_call(f'rm -r {functional}')
            if not os.path.isdir(functional):
                utils.system_call(f'mkdir {functional}')

        if workflow_settings['alpha_from_file']:
            utils.system_call('cp file_alpharef*.txt ki/')

        for functional in functionals:
            print(f'\n{functional.upper()} CALCULATION')

            # Avoid overwriting defaults
            local_workflow_settings = copy.deepcopy(workflow_settings)
            cp_calc = copy.deepcopy(cp_master_calc)

            local_workflow_settings['functional'] = functional

            # For pKIPZ/KIPZ, use KI as a starting point
            if functional == 'pkipz':
                local_workflow_settings['from_scratch'] = False
                local_workflow_settings['calculate_alpha'] = False
                local_workflow_settings['alpha_from_file'] = False
            elif functional == 'kipz':
                local_workflow_settings['init_density'] = 'ki'
                local_workflow_settings['init_manifold'] = 'ki'
                local_workflow_settings['alpha_from_file'] = True

            # Check that the workflow settings are valid
            check_settings(local_workflow_settings)

            # Change to relevant subdirectory
            os.chdir(functional)

            # Amend the psuedo dir
            if isinstance(cp_calc.pseudo_dir, str) and cp_calc.pseudo_dir[0] != '/':
                cp_calc.pseudo_dir = '../' + cp_calc.pseudo_dir

            # Run workflow for this particular functional
            solved_calc = kc_with_cp.run(cp_calc, local_workflow_settings)

            # Return to the base directory
            os.chdir('..')

            # Provide the pKIPZ and KIPZ calculations with a KI starting point
            if functional == 'ki':
                # pKIPZ
                utils.system_call('cp -r ki/* pkipz/')
                utils.system_call('cp pkipz/final/file_alpharef*.txt pkipz/')

                # KIPZ
                utils.system_call('cp -r ki/final kipz/init')
                utils.system_call(f'cp -r {cp_calc.outdir} kipz/')
                utils.system_call('mv kipz/init/ki_final.cpi kipz/init/ki_init.cpi')
                utils.system_call('mv kipz/init/ki_final.cpo kipz/init/ki_init.cpo')
                utils.system_call('cp -r ki/final/file_alpharef* kipz/')
        return solved_calc

    elif workflow_settings['functional'] in ['ki', 'kipz', 'pkipz']:
        return kc_with_cp.run(cp_master_calc, workflow_settings)

    elif workflow_settings['functional'] == 'pbe':
        return pbe_with_cp.run(cp_master_calc, workflow_settings)

    return
