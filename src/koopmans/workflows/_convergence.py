"""

convergence workflow object for koopmans

Written by Edward Linscott Oct 2020
Converted to a workflow object Nov 2020

"""

from __future__ import annotations

import copy
import itertools
import shutil
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from typing import (Any, Callable, Dict, Generic, List, Optional, Type,
                    TypeVar, Union, cast)

import numpy as np
import numpy.typing as npt

from koopmans import cell, utils

from ._workflow import Workflow

# T is either an int or a float or a list of ints/floats
T = TypeVar('T', int, float, List[int], List[float])


@dataclass
class ConvergenceParameter(Generic[T]):
    name: str
    increment: T
    get_value: Callable[[Workflow], T]
    set_value: Callable[[Workflow, T], None]
    initial_value: Optional[T] = None
    length: int = 3
    converged_value: Optional[T] = None

    def extend(self):
        self.length += 1

    def __len__(self) -> int:
        return self.length

    @property
    def values(self) -> List[T]:
        if self.initial_value is None:
            raise ValueError(
                f'{self.__class__.__name__}.initial_value has not been set. Use {self.__class__.__name__}.get_initial_value() to set it.')

        if isinstance(self.initial_value, list):
            return [[v + i * delta for v, delta in zip(self.initial_value, self.increment)] for i in range(self.length)]
        else:
            return [self.initial_value + i * self.increment for i in range(self.length)]

    def get_initial_value(self, workflow: Workflow):
        if self.initial_value is not None:
            raise ValueError(
                f'Do not call {self.__class__.__name__}.get_initial_value() once {self.__class__.__name__}.initial_value has been initialized')
        self.initial_value = self.get_value(workflow)


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


def fetch_observable_factory(obs_str: str) -> Callable[[Workflow], float]:
    if obs_str == 'eps_inf':
        return fetch_eps_inf
    else:
        return partial(fetch_result_default, observable=obs_str)


def conv_param_celldm1(increment: float = 1.0, **kwargs) -> ConvergenceParameter:

    def get_value(wf: Workflow) -> float:

        params = cell.cell_to_parameters(wf.atoms.cell)
        assert isinstance(params['celldms'], dict)
        assert 1 in params['celldms']
        import ipdb
        ipdb.set_trace()
        return params['celldms'][1]

    def set_value(wf: Workflow, value: float) -> None:
        params = cell.cell_to_parameters(wf.atoms.cell)
        assert isinstance(params['celldms'], dict)
        assert 1 in params['celldms']
        params['celldms'][1] = value
        wf.atoms.cell = cell.parameters_to_cell(**params)

    return ConvergenceParameter('cell_size', increment, get_value, set_value, **kwargs)


def conv_param_ecutwfc(increment: float = 10.0, **kwargs) -> ConvergenceParameter:

    def set_value(wf: Workflow, value: float) -> None:
        set_calc_value(wf, value, key='ecutwfc')
        set_calc_value(wf, 4 * value, key='ecutrho')

    get_value = cast(Callable[[Workflow], float], partial(get_calc_value, key='ecutwfc'))

    return ConvergenceParameter('ecutwfc', increment, get_value, set_value, **kwargs)


def conv_param_nbnd(increment: int = 1, **kwargs) -> ConvergenceParameter:

    set_value = cast(Callable[[Workflow, int], None], partial(set_calc_value, key='nbnd'))

    get_value = cast(Callable[[Workflow], int], partial(get_calc_value, key='nbnd'))

    return ConvergenceParameter('nbnd', increment, get_value, set_value, **kwargs)


def conv_param_kgrid(increment: List[int] = [2, 2, 2], **kwargs) -> ConvergenceParameter:

    def get_value(wf: Workflow) -> List[int]:
        assert wf.kpoints.grid
        return wf.kpoints.grid

    def set_value(wf: Workflow, value: List[int]) -> None:
        wf.kpoints.grid = value

    return ConvergenceParameter('kgrid', increment, get_value, set_value, **kwargs)


def ConvergenceParameterFactory(conv_param, **kwargs) -> ConvergenceParameter:

    if conv_param == 'celldm1':
        return conv_param_celldm1(**kwargs)
    elif conv_param == 'ecutwfc':
        return conv_param_ecutwfc(**kwargs)
    elif conv_param == 'nbnd':
        return conv_param_nbnd(**kwargs)
    elif conv_param == 'kgrid':
        return conv_param_kgrid(**kwargs)
    else:
        raise NotImplementedError(f'Convergence with respect to {conv_param} has not been directly implemented. You '
                                  'can still perform a convergence calculation with respect to this parameter, but '
                                  'you must first create an appropriate ConvergenceParameter object and then '
                                  'construct your ConvergenceWorkflow using the ConvergenceWorkflowFactory')


class ConvergenceWorkflow(Workflow):

    def __init__(self, subworkflow_class: Type[Workflow], convergence_parameters: List[ConvergenceParameter] = [], fetch_observable: Optional[Callable[[Workflow], float]] = None, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._subworkflow_class = subworkflow_class
        self._convergence_parameters = convergence_parameters
        self._fetch_observable = fetch_observable

        # Set the initial value for each of the convergence parameters
        for c in self._convergence_parameters:
            c.get_initial_value(self)

    def _run(self, initial_depth: int = 3) -> None:

        # Deferred import to allow for monkeypatching
        from koopmans import workflows

        # Check that everything has been initialized
        if self._fetch_observable is None:
            raise ValueError(
                f'{self.__class__.__name__} has not been provided with a "convergence_observable" to converge')

        if len(self._convergence_parameters) == 0:
            raise ValueError(
                f'{self.__class__.__name__} has not been provided with any "convergence_parameters" with which to perform convergence')

        if self.parameters.from_scratch:
            for c in self._convergence_parameters:
                for path in Path().glob(c.name + '*'):
                    shutil.rmtree(str(path))

        # Create array for storing calculation results
        results = np.empty([initial_depth for _ in self._convergence_parameters])
        results[:] = np.nan

        # Continue to increment the convergence parameters until convergence in the convergence observable
        # is achieved. A lot of the following code is quite obtuse, in order to make it work for an arbitrary
        # number of convergence parameters

        # Record the alpha values for the original calculation
        provide_alpha = any([p.name == 'nbnd' for p in self._convergence_parameters]) \
            and self.parameters.functional in ['ki', 'kipz', 'pkipz', 'all'] \
            and self.parameters.alpha_from_file
        if provide_alpha:
            master_alphas = utils.read_alpha_file(directory=Path())
            if self.parameters.orbital_groups is None:
                self.parameters.orbital_groups = list(range(self.calculator_parameters['kcp'].nbnd))
            master_orbital_groups = copy.deepcopy(self.parameters.orbital_groups)

        while True:

            # Loop over all possible permutations of convergence parameter settings
            for indices in itertools.product(*[range(len(x.values)) for x in self._convergence_parameters]):
                # If we already have results for this permutation, don't recalculate it
                if not np.isnan(results[tuple(indices)]):
                    continue

                # Create a new subworkflow
                subwf = self._subworkflow_class.fromparent(self)

                # For each parameter we're converging wrt...
                header = ''
                subdir = Path()
                for index, parameter in zip(indices, self._convergence_parameters):
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
                results[indices] = self._fetch_observable(subwf)

            # Check convergence
            converged = np.empty([len(p) for p in self._convergence_parameters], dtype=bool)
            converged[:] = False

            converged_result = results[tuple([-1 for _ in self._convergence_parameters])]
            threshold = self.parameters.convergence_threshold
            subarray_slice: List[Union[int, slice]]
            for indices in itertools.product(*[range(len(p.values)) for p in self._convergence_parameters]):
                subarray_slice = [slice(None) for _ in range(len(indices))]
                for i_index, index in enumerate(indices):
                    subarray_slice[i_index] = slice(index, None)
                if np.all(np.abs(results[tuple(subarray_slice)] - converged_result) < threshold):
                    converged[tuple(indices)] = True

            # Check convergence, ignoring the calculations with the most fine parameters because
            # we have nothing to compare them against
            subarray_slice = [slice(0, -1) for _ in self._convergence_parameters]
            if np.any(converged[tuple(subarray_slice)]):
                # Work out which was the fastest calculation, and propose those parameters

                # First, find the indices of the converged array
                slice_where: List[Union[slice, int]] = [slice(None) for _ in self._convergence_parameters]
                slice_where[-1] = 0
                converged_indices = np.array(np.where(converged[tuple(subarray_slice)]))[tuple(slice_where)]

                # Extract the corresponding parameters
                for index, parameter in zip(converged_indices, self._convergence_parameters):
                    parameter.converged_value = parameter.values[index]

                self.print('\n Converged parameters are '
                           + ', '.join([f'{p.name} = {p.converged_value}' for p in self._convergence_parameters]))

                return
            else:
                # Work out which parameters are yet to converge, and line up more calculations
                # for increased values of those parameters

                new_array_shape = list(np.shape(results))
                new_array_slice: List[Union[int, slice]] = [slice(None) for _ in indices]
                self.print('Progress update', style='heading')
                for index, param in enumerate(self._convergence_parameters):
                    subarray_slice = [slice(None) for _ in self._convergence_parameters]
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


def ConvergenceWorkflowFactory(subworkflow: Workflow, parameters: List[Union[str, ConvergenceParameter]], fetch_observable: Union[str, Callable[[Workflow], float]]) -> Workflow:
    '''
    Generates a ConvergenceWorkflow based on...
        subworkflow_class: the class of the Workflow that we will repeatedly run
        parameters: a list of parameters with respect to which we will perform the convergence
        fetch_observable: a function with which we extract the observable from the child workflow
    '''

    # Ensure parameters are all ConvergenceParameters
    parameters = [p if isinstance(p, ConvergenceParameter) else ConvergenceParameterFactory(p) for p in parameters]

    # Ensure fetch_observable is a Callable
    fetch_observable = fetch_observable_factory(fetch_observable) if isinstance(
        fetch_observable, str) else fetch_observable

    # Co-opt the fromparent method to use the subworkflow settings to initialize a ConvergenceWorkflow...
    wf = ConvergenceWorkflow.fromparent(subworkflow, subworkflow_class=subworkflow.__class__,
                                        convergence_parameters=parameters, fetch_observable=fetch_observable)
    wf.parent = None  # ... making sure we remove wf.parent immediately afterwards

    return wf
