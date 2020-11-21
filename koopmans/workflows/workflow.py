"""

Generic workflow functions for python_KI

Written by Edward Linscott Jan 2020
KI/pKIPZ/KIPZ moved to kc_with_cp.py Oct 2020

"""

import os
import copy
import numpy as np
import itertools
from koopmans import utils, io
from koopmans.ase import write_json
from koopmans.defaults import load_defaults
from koopmans.workflows import kc_with_cp, pbe_with_cp


def run_convergence(workflow_settings, master_calcs_dct, initial_depth=3):

    from koopmans.config import from_scratch

    if 'cp' in master_calcs_dct:
        cp_master_calc = master_calcs_dct['cp']
    else:
        raise ValueError(
            'run_convergence() has not been generalised beyond cp.x')

    increments = {'cell_size': 0.1, 'ecutwfc': 10, 'empty_states_nbnd': 1}

    if from_scratch:
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
            raise NotImplementedError(
                f'Convergence wrt {parameter} has not yet been implemented')

        if parameter == 'cell_size':
            param_dict[parameter] = [1 + i * increments['cell_size']
                                     for i in range(initial_depth)]
        else:
            attr = getattr(cp_master_calc, parameter, None)
            if attr is None:
                raise AttributeError(
                    f'In order to converge wrt {parameter} specify a baseline value for it')
            param_dict[parameter] = [attr + i * increments[parameter]
                                     for i in range(initial_depth)]

    # Create array for storing calculation results
    results = np.empty([initial_depth for _ in param_dict])
    results[:] = np.nan

    # Continue to increment the convergence parameters until convergence in the convergence observable
    # is achieved. A lot of the following code is quite obtuse, in order to make it work for an arbitrary
    # number of convergence parameters

    # Record the alpha values for the original calculation
    provide_alpha = 'empty_states_nbnd' in param_dict and workflow_settings['functional'] in ['ki', 'kipz', 'pkipz', 'all'] \
        and workflow_settings['alpha_from_file']
    if provide_alpha:
        master_alphas = io.read_alpharef(directory='.')
        if workflow_settings['orbital_groups'] is None:
            workflow_settings['orbital_groups'] = list(
                range(cp_master_calc.nelec // 2 + cp_master_calc.empty_states_nbnd))
        master_orbital_groups = copy.deepcopy(
            workflow_settings['orbital_groups'])

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
                subdir = f'{param}_{value_str}'.replace(
                    ' ', '_').replace('.', 'd')
                if not os.path.isdir(subdir):
                    utils.system_call(f'mkdir {subdir}')
                os.chdir(subdir)

                if param == 'cell_size':
                    cp_calc._ase_calc.atoms.cell *= value
                else:
                    setattr(cp_calc, param, value)
                if param == 'ecutwfc':
                    setattr(cp_calc, 'ecutrho', 4 * value)

            if provide_alpha:
                # Update alpha files and orbitals
                extra_orbitals = cp_calc.empty_states_nbnd - cp_master_calc.empty_states_nbnd
                filling = [True] * (cp_calc.nelec // 2) + \
                    [False] * cp_calc.empty_states_nbnd
                alphas = master_alphas + [master_alphas[-1]
                                          for _ in range(extra_orbitals)]
                workflow_settings['orbital_groups'] = master_orbital_groups + \
                    [master_orbital_groups[-1] for _ in range(extra_orbitals)]
                io.write_alpharef(alphas, directory='.', filling=filling)

            print(header.rstrip(', '))

            # Perform calculation
            calcs_dct = copy.deepcopy(master_calcs_dct)
            calcs_dct['cp'] = cp_calc
            solved_calc = run_singlepoint(workflow_settings, calcs_dct)

            # Store the result
            obs = workflow_settings['convergence_observable']
            obs = obs.replace('total ', '').replace(' ', '_')
            if obs not in solved_calc.results:
                raise ValueError(
                    f'{solved_calc.name} has not returned a value for {obs}')
            result = solved_calc.results[obs]
            results[indices] = result

            # Move back to the base directory:
            for _ in param_dict:
                os.chdir('..')

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
            indices = np.array(np.where(converged[subarray_slice]))[
                tuple(slice_where)]

            # Extract the corresponding parameters
            converged_parameters = {}
            for index, param in zip(indices, param_dict.keys()):
                converged_parameters[param] = param_dict[param][index]

            print('Converged parameters are '
                  ', '.join([f'{k} = {v}' for k, v in converged_parameters.items()]))

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
                    setattr(cp_master_calc, 'ecutrho', 4 * value)

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
                param_dict[param].append(
                    param_dict[param][-1] + increments[param])
                new_array_shape[index] += 1
                new_array_slice[index] = slice(None, -1)

            new_results = np.empty(new_array_shape)
            new_results[:] = np.nan
            new_results[tuple(new_array_slice)] = results
            results = new_results


def run_singlepoint(workflow_settings, calcs_dct):

    from koopmans.config import from_scratch

    if workflow_settings['functional'] == 'all':
        # if 'all', create subdirectories and run
        functionals = ['ki', 'pkipz', 'kipz']

        # Make separate directories for KI, pKIPZ, and KIPZ
        for functional in functionals:
            if from_scratch and os.path.isdir(functional):
                utils.system_call(f'rm -r {functional}')
            if not os.path.isdir(functional):
                utils.system_call(f'mkdir {functional}')

        if workflow_settings['alpha_from_file']:
            utils.system_call('cp file_alpharef*.txt ki/')

        for functional in functionals:
            print(f'\n{functional.upper().replace("PKIPZ", "pKIPZ")} CALCULATION')

            # Avoid overwriting defaults
            local_workflow_settings = copy.deepcopy(workflow_settings)

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

            # Change to relevant subdirectory
            os.chdir(functional)

            # Run workflow for this particular functional
            solved_calc = kc_with_cp.run(local_workflow_settings, calcs_dct)

            # Return to the base directory
            os.chdir('..')

            # Provide the pKIPZ and KIPZ calculations with a KI starting point
            if functional == 'ki':
                # pKIPZ
                utils.system_call('cp -r ki/* pkipz/')
                utils.system_call('cp pkipz/final/file_alpharef*.txt pkipz/')

                # KIPZ
                utils.system_call('cp -r ki/final kipz/init')
                utils.system_call(f'cp -r {solved_calc.outdir} kipz/')
                utils.system_call(
                    'mv kipz/init/ki_final.cpi kipz/init/ki_init.cpi')
                utils.system_call(
                    'mv kipz/init/ki_final.cpo kipz/init/ki_init.cpo')
                utils.system_call('cp -r ki/final/file_alpharef* kipz/')
        return solved_calc

    elif workflow_settings['functional'] in ['ki', 'kipz', 'pkipz']:
        return kc_with_cp.run(workflow_settings, calcs_dct)

    elif workflow_settings['functional'] == 'pbe':
        return pbe_with_cp.run(workflow_settings, calcs_dct)

    return
