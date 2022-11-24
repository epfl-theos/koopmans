"""

convergence workflow object for koopmans

Written by Edward Linscott Oct 2020
Converted to a workflow object Nov 2020

"""

import copy
import itertools
import shutil
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from typing import Any, Callable, Generic, List, Optional, TypeVar, Union, cast

import numpy as np

from koopmans import utils

from ._workflow import Workflow

T = TypeVar('T')


@dataclass
class ConvergenceParameter(Generic[T]):
    name: str
    increment: T
    get_value: Callable[[Workflow], T]
    set_value: Callable[[Workflow, T], None]
    values: List[T] = field(default_factory=list)
    converged_value: Optional[T] = None

    def extend(self):
        if isinstance(self.values[-1], list):
            self.values.append([x + self.increment for x in self.values[-1]])
        else:
            self.values.append(self.values[-1] + self.increment)


def get_calc_value(wf: Workflow, key: str) -> None:
    '''
    Gets a particular calculator setting from a workflow
    '''
    values = set([])
    for settings in wf.calculator_parameters.values():
        if settings.is_valid(key):
            val = settings[key]
            if val is None:
                raise AttributeError(f'In order to converge wrt {key}, specify a baseline value for it')
            values.add(val)
    if len(values) > 1:
        raise ValueError(f'{key} has different values for different calculators. This is not compatible with the'
                         'convergence workflow')
    return values.pop()


def set_calc_value(wf: Workflow, value: Union[int, float], key: str) -> None:
    '''
    Sets a particular calculator setting for every calculator in the workflow
    '''
    for settings in wf.calculator_parameters.values():
        if settings.is_valid(key):
            settings[key] = value


def fetch_result_default(wf: Workflow, observable: str) -> float:
    calc = wf.calculations[-1]
    observable = observable.replace('total ', '').replace(' ', '_')
    if observable not in calc.results:
        raise ValueError(
            f'{calc.prefix} has not returned a value for {observable}')
    return calc.results[observable]


def fetch_eps_inf(wf: Workflow) -> float:
    calc = wf.calculations[-1]
    return np.mean(np.diag(calc.results['dielectric tensor']))


class ConvergenceWorkflow(Workflow):

    def _run(self, initial_depth: int = 3) -> None:

        # Deferred import to allow for monkeypatchings
        from koopmans import workflows

        convergence_parameters = []

        for conv_param in self.parameters.convergence_parameters:

            get_value: Callable[[Workflow], Any]
            set_value: Callable[[Workflow, Any], None]
            values: List[Any]

            if conv_param == 'cell_size':
                increment = 0.1

                def get_value(wf: Workflow) -> float:
                    mask = self.atoms.cell != 0
                    cell_ratio = wf.atoms.cell[mask] / self.atoms.cell[mask]
                    return list(set(cell_ratio))[0]

                def set_value(wf: Workflow, value: float) -> None:
                    wf.atoms.cell = self.atoms.cell * value
                values = [get_value(self) + i * increment for i in range(initial_depth)]

            elif conv_param == 'ecutwfc':
                increment = 10

                def set_value(wf: Workflow, value: float) -> None:
                    set_calc_value(wf, value, key='ecutwfc')
                    set_calc_value(wf, 4 * value, key='ecutrho')
                get_value = cast(Callable[[Workflow], float], partial(get_calc_value, key='ecutwfc'))
                values = [get_value(self) + i * increment for i in range(initial_depth)]
            elif conv_param == 'nbnd':
                increment = 1
                set_value = cast(Callable[[Workflow, int], None], partial(set_calc_value, key='nbnd'))
                get_value = cast(Callable[[Workflow], int], partial(get_calc_value, key='nbnd'))
                values = [get_value(self) + i * increment for i in range(initial_depth)]

            elif conv_param == 'kgrid':
                increment = 2

                def get_value(wf: Workflow) -> List[int]:
                    assert wf.kpoints.grid
                    return wf.kpoints.grid

                def set_value(wf: Workflow, value: List[int]) -> None:
                    wf.kpoints.grid = value
                values = [[x + i * increment for x in get_value(self)] for i in range(initial_depth)]

            else:
                raise NotImplementedError(f'Convergence with respect to {conv_param} has not yet been implemented')
            convergence_parameters.append(ConvergenceParameter(conv_param, increment, get_value, set_value, values))

        if self.parameters.convergence_observable == 'eps_inf':
            fetch_result = fetch_eps_inf
        else:
            fetch_result = partial(fetch_result_default, observable=self.parameters.convergence_observable)

        if self.parameters.from_scratch:
            for c in convergence_parameters:
                for path in Path().glob(c.name + '*'):
                    shutil.rmtree(str(path))

        # Create array for storing calculation results
        results = np.empty([initial_depth for _ in convergence_parameters])
        results[:] = np.nan

        # Continue to increment the convergence parameters until convergence in the convergence observable
        # is achieved. A lot of the following code is quite obtuse, in order to make it work for an arbitrary
        # number of convergence parameters

        # Record the alpha values for the original calculation
        provide_alpha = any([p.name == 'nbnd' for p in convergence_parameters]) \
            and self.parameters.functional in ['ki', 'kipz', 'pkipz', 'all'] \
            and self.parameters.alpha_from_file
        if provide_alpha:
            master_alphas = utils.read_alpha_file(directory=Path())
            if self.parameters.orbital_groups is None:
                self.parameters.orbital_groups = list(range(self.calculator_parameters['kcp'].nbnd))
            master_orbital_groups = copy.deepcopy(self.parameters.orbital_groups)

        while True:

            # Loop over all possible permutations of convergence parameter settings
            for indices in itertools.product(*[range(len(x.values)) for x in convergence_parameters]):
                # If we already have results for this permutation, don't recalculate it
                if not np.isnan(results[tuple(indices)]):
                    continue

                # Create a new workflow
                subwf: workflows.Workflow
                if self.parameters.convergence_observable == 'eps_inf':
                    subwf = workflows.DFTPhWorkflow.fromparent(self)
                else:
                    subwf = workflows.SinglepointWorkflow.fromparent(self)

                # For each parameter we're converging wrt...
                header = ''
                subdir = Path()
                for index, parameter in zip(indices, convergence_parameters):
                    value = parameter.values[index]
                    if isinstance(value, float):
                        value_str = f'{value:.1f}'
                    else:
                        value_str = str(value)
                    header += f'{parameter.name} = {value_str}, '

                    if isinstance(value, list):
                        value_str = ''.join([str(x) for x in value])

                    # Create new working directory
                    subdir /= f'{parameter.name}_{value_str}'.replace(' ', '_').replace('.', 'd')

                    # Set the value
                    parameter.set_value(subwf, value)

                if provide_alpha:
                    # Update alpha files and orbitals
                    raise NotImplementedError('TODO')
                    # extra_orbitals = subwf.calculator_parameters.nbnd - kcp_master_params.nbnd
                    # filling = [True] * (kcp_params.nelec // 2) + \
                    #     [False] * (kcp_params.nbnd - kcp_params.nelec // 2)
                    # alphas = master_alphas + [master_alphas[-1] for _ in range(extra_orbitals)]
                    # self.parameters.orbital_groups = master_orbital_groups + \
                    #     [master_orbital_groups[-1]
                    #         for _ in range(extra_orbitals)]
                    # utils.write_alpha_file(directory=Path(), alphas=alphas, filling=filling)

                self.print(header.rstrip(', '), style='subheading')

                # Perform calculation
                subwf.run(subdirectory=subdir)

                # Store the result
                results[indices] = fetch_result(subwf)

            # Check convergence
            converged = np.empty([len(p.values) for p in convergence_parameters], dtype=bool)
            converged[:] = False

            converged_result = results[tuple([-1 for _ in convergence_parameters])]
            threshold = self.parameters.convergence_threshold
            subarray_slice: List[Union[int, slice]]
            for indices in itertools.product(*[range(len(p.values)) for p in convergence_parameters]):
                subarray_slice = [slice(None) for _ in range(len(indices))]
                for i_index, index in enumerate(indices):
                    subarray_slice[i_index] = slice(index, None)
                if np.all(np.abs(results[tuple(subarray_slice)] - converged_result) < threshold):
                    converged[tuple(indices)] = True

            # Check convergence, ignoring the calculations with the most fine parameters because
            # we have nothing to compare them against
            subarray_slice = [slice(0, -1) for _ in convergence_parameters]
            if np.any(converged[tuple(subarray_slice)]):
                # Work out which was the fastest calculation, and propose those parameters

                # First, find the indices of the converged array
                slice_where: List[Union[slice, int]] = [slice(None) for _ in convergence_parameters]
                slice_where[-1] = 0
                converged_indices = np.array(np.where(converged[tuple(subarray_slice)]))[tuple(slice_where)]

                # Extract the corresponding parameters
                for index, parameter in zip(converged_indices, convergence_parameters):
                    parameter.converged_value = parameter.values[index]

                self.print('\n Converged parameters are '
                           + ', '.join([f'{p.name} = {p.converged_value}' for p in convergence_parameters]))

                return
            else:
                # Work out which parameters are yet to converge, and line up more calculations
                # for increased values of those parameters

                new_array_shape = list(np.shape(results))
                new_array_slice: List[Union[int, slice]] = [slice(None) for _ in indices]
                self.print('Progress update', style='heading')
                for index, param in enumerate(convergence_parameters):
                    subarray_slice = [slice(None) for _ in convergence_parameters]
                    subarray_slice[index] = slice(0, -1)

                    if np.any(converged[tuple(subarray_slice)]):
                        self.print(f'{param.name} appears converged')
                        continue

                    # Add one finer value to try
                    param.extend()
                    new_array_shape[index] += 1
                    new_array_slice[index] = slice(None, -1)
                    self.print(f'{param.name} still appears unconverged, will reattempt using a finer value')

                new_results = np.empty(new_array_shape)
                new_results[:] = np.nan
                new_results[tuple(new_array_slice)] = results
                results = new_results
