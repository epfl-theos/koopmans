"""

convergence workflow object for koopmans

Written by Edward Linscott Oct 2020
Converted to a workflow object Nov 2020

"""

import copy
import itertools
from pathlib import Path
from typing import List, Union

import numpy as np

from koopmans import utils

from ._workflow import Workflow


class ConvergenceWorkflow(Workflow):

    def _run(self, initial_depth: int = 3) -> None:

        # Deferred import to allow for monkeypatchings
        from koopmans import workflows

        if 'kcp' in self.master_calc_params:
            kcp_master_params = self.master_calc_params['kcp']
        else:
            raise NotImplementedError(
                'Convergence.run() has not been generalized beyond kcp.x')

        increments = {'cell_size': 0.1, 'ecutwfc': 10, 'nbnd': 1}

        if self.parameters.from_scratch:
            for key in increments.keys():
                utils.system_call(f'rm -r {key}* 2> /dev/null', False)

        # Populate the param_dict
        param_dict = {}

        # Convert to list
        convergence_parameters = self.parameters.convergence_parameters
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
                attr = getattr(kcp_master_params, parameter, None)
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
        provide_alpha = 'nbnd' in param_dict \
            and self.parameters.functional in ['ki', 'kipz', 'pkipz', 'all'] \
            and self.parameters.alpha_from_file
        if provide_alpha:
            master_alphas = utils.read_alpha_file(directory=Path())
            if self.parameters.orbital_groups is None:
                self.parameters.orbital_groups = list(range(kcp_master_params.nbnd))
            master_orbital_groups = copy.deepcopy(self.parameters.orbital_groups)

        while True:

            # Loop over all possible permutations of convergence parameter settings
            for indices in itertools.product(*[range(len(x)) for x in param_dict.values()]):
                # If we already have results for this permutation, don't recalculate it
                if not np.isnan(results[tuple(indices)]):
                    continue

                # Create duplicate kcp calculator settings and atoms
                kcp_params = copy.deepcopy(self.master_calc_params['kcp'])
                atoms = copy.deepcopy(self.atoms)

                # For each parameter we're converging wrt...
                header = ''
                subdir = ''
                for index, param, values in zip(indices, param_dict.keys(), param_dict.values()):
                    value = values[index]
                    if isinstance(value, int):
                        value_str = str(value)
                    elif isinstance(value, float):
                        value_str = '{:.1f}'.format(value)
                    header += f'{param} = {value_str}, '

                    # Create new working directory
                    subdir += f'{param}_{value_str}/'.replace(
                        ' ', '_').replace('.', 'd')

                    if param == 'cell_size':
                        atoms.cell *= value
                    else:
                        setattr(kcp_params, param, value)
                    if param == 'ecutwfc':
                        setattr(kcp_params, 'ecutrho', 4 * value)

                if provide_alpha:
                    # Update alpha files and orbitals
                    extra_orbitals = kcp_params.nbnd - kcp_master_params.nbnd
                    filling = [True] * (kcp_params.nelec // 2) + \
                        [False] * (kcp_params.nbnd - kcp_params.nelec // 2)
                    alphas = master_alphas + [master_alphas[-1] for _ in range(extra_orbitals)]
                    self.parameters.orbital_groups = master_orbital_groups + \
                        [master_orbital_groups[-1]
                            for _ in range(extra_orbitals)]
                    utils.write_alpha_file(directory=Path(), alphas=alphas, filling=filling)

                self.print(header.rstrip(', '), style='subheading')

                # Perform calculation
                wf_kwargs = self.wf_kwargs
                wf_kwargs['atoms'] = atoms
                wf_kwargs['master_calc_params']['kcp'] = kcp_params
                singlepoint = workflows.SinglepointWorkflow(**wf_kwargs)
                self.run_subworkflow(singlepoint, subdirectory=subdir)
                solved_calc = singlepoint.calculations[-1]

                # Store the result
                obs = self.parameters.convergence_observable
                obs = obs.replace('total ', '').replace(' ', '_')
                if obs not in solved_calc.results:
                    raise ValueError(
                        f'{solved_calc.prefix} has not returned a value for {obs}')
                result = solved_calc.results[obs]
                results[indices] = result

            # Check convergence
            converged = np.empty([len(x) for x in param_dict.values()], dtype=bool)
            converged[:] = False

            converged_result = results[tuple([-1 for _ in param_dict])]
            threshold = self.parameters.convergence_threshold
            subarray_slice: List[Union[int, slice]]
            for indices in itertools.product(*[range(len(x)) for x in param_dict.values()]):
                subarray_slice = [slice(None) for _ in range(len(indices))]
                for i_index, index in enumerate(indices):
                    subarray_slice[i_index] = slice(index, None)
                if np.all(np.abs(results[tuple(subarray_slice)] - converged_result) < threshold):
                    converged[tuple(indices)] = True

            # Check convergence, ignoring the calculations with the most fine parameters because
            # we have nothing to compare them against
            subarray_slice = [slice(0, -1) for _ in param_dict]
            if np.any(converged[tuple(subarray_slice)]):
                # Work out which was the fastest calculation, and propose those parameters

                # First, find the indices of the converged array
                slice_where: List[Union[slice, int]] = [slice(None) for _ in param_dict]
                slice_where[-1] = 0
                converged_indices = np.array(np.where(converged[tuple(subarray_slice)]))[tuple(slice_where)]

                # Extract the corresponding parameters
                converged_parameters = {}
                for index, param in zip(converged_indices, param_dict.keys()):
                    converged_parameters[param] = param_dict[param][index]

                self.print('\n Converged parameters are '
                           + ', '.join([f'{k} = {v}' for k, v in converged_parameters.items()]))

                return
            else:
                # Work out which parameters are yet to converge, and line up more calculations
                # for increased values of those parameters

                new_array_shape = list(np.shape(results))
                new_array_slice: List[Union[int, slice]] = [slice(None) for _ in indices]
                self.print('Progress update', style='heading')
                for index, param in enumerate(param_dict):
                    subarray_slice = [slice(None) for _ in param_dict]
                    subarray_slice[index] = slice(0, -1)
                    if np.any(converged[tuple(subarray_slice)]):
                        self.print(f'{param} appears converged')
                        continue
                    new_val = param_dict[param][-1] + increments[param]
                    param_dict[param].append(new_val)
                    new_array_shape[index] += 1
                    new_array_slice[index] = slice(None, -1)
                    self.print(f'{param} still appears unconverged, will reattempt using a finer value')

                new_results = np.empty(new_array_shape)
                new_results[:] = np.nan
                new_results[tuple(new_array_slice)] = results
                results = new_results
