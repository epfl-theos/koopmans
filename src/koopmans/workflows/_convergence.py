"""

convergence workflow object for koopmans

Written by Edward Linscott Oct 2020
Converted to a workflow object Nov 2020

"""

from __future__ import annotations

import copy
import itertools
import shutil
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import (Any, Callable, Dict, Generator, Generic, List, Optional,
                    Type, TypeVar, Union, cast)

import numpy as np

from koopmans import cell, utils
from koopmans.outputs import OutputModel
from koopmans.status import Status

from ._workflow import Workflow

# T is either an int or a float or a list of ints/floats
T = TypeVar('T', int, float, List[int], List[float])


def generate_values_via_addition(initial_value: T, increment: T, length: int) -> List[T]:
    if isinstance(initial_value, list):
        return [[v + i * delta for v, delta in zip(initial_value, increment)] for i in range(length)]
    else:
        return [initial_value + i * increment for i in range(length)]


@dataclass
class ConvergenceVariable(Generic[T]):
    name: str
    increment: T
    get_value: Callable[[Workflow], T]
    set_value: Callable[[Workflow, T], None]
    initial_value: Optional[T] = None
    length: int = 3
    generate_values: Callable[[T, T, int], List[T]] = generate_values_via_addition
    converged_value: Optional[T] = None

    def extend(self):
        self.length += 1

    def __len__(self) -> int:
        return self.length

    @property
    def values(self) -> List[T]:
        if self.initial_value is None:
            raise ValueError(
                f'`{self.__class__.__name__}.initial_value` has not been set. Use `{self.__class__.__name__}.get_initial_value()` to set it.')

        return self.generate_values(self.initial_value, self.increment, self.length)

    def get_initial_value(self, workflow: Workflow):
        if self.initial_value is not None:
            raise ValueError(
                f'Do not call `{self.__class__.__name__}.get_initial_value()` once `{self.__class__.__name__}.initial_value` has been initialized')
        self.initial_value = self.get_value(workflow)

    def todict(self) -> Dict[str, Any]:
        dct = dict(self.__dict__)
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @classmethod
    def fromdict(cls, dct: Dict[str, Any]):
        return cls(**dct)


def get_calculator_parameter(wf: Workflow, key: str) -> None:
    '''
    Gets a particular calculator setting from a workflow
    '''
    values = set([])
    for settings in wf.calculator_parameters.values():
        if settings.is_valid(key):
            val = settings[key]
            if val is None:
                raise AttributeError(
                    f'In order to converge wrt `{key}`, specify a baseline value for it')
            values.add(val)
    if len(values) > 1:
        raise ValueError(f'`{key}` has different values for different calculators. This is not compatible with the'
                         'convergence workflow')
    return values.pop()


def set_calculator_parameter(wf: Workflow, key: str, value: Union[int, float]) -> None:
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
            f'`{calc.prefix}` has not returned a value for `{observable}`')
    return calc.results[observable]


def fetch_eps_inf(wf: Workflow) -> float:
    calc = wf.calculations[-1]
    return np.mean(np.diag(calc.results['dielectric tensor']))


def ObservableFactory(obs_str: str) -> Callable[[Workflow], float]:
    if obs_str == 'eps_inf':
        return fetch_eps_inf
    else:
        return partial(fetch_result_default, observable=obs_str)


def _get_celldm1(wf: Workflow) -> float:

    params = cell.cell_to_parameters(wf.atoms.cell)
    assert 1 in params['celldms']
    return params['celldms'][1]


def _set_celldm1(wf: Workflow, value: float) -> None:
    params = cell.cell_to_parameters(wf.atoms.cell)
    params['celldms'][1] = value
    wf.atoms.cell = cell.parameters_to_cell(**params)


def conv_var_celldm1(increment: float = 1.0, **kwargs) -> ConvergenceVariable:

    return ConvergenceVariable('celldm1', increment, _get_celldm1, _set_celldm1, **kwargs)


def _set_ecutwfc(wf: Workflow, value: float) -> None:
    set_calculator_parameter(wf, 'ecutwfc', value)
    set_calculator_parameter(wf, 'ecutrho', 4 * value)


_get_ecutwfc = cast(Callable[[Workflow], float], partial(get_calculator_parameter, key='ecutwfc'))


def conv_var_ecutwfc(increment: float = 10.0, **kwargs) -> ConvergenceVariable:

    return ConvergenceVariable('ecutwfc', increment, _get_ecutwfc, _set_ecutwfc, **kwargs)


_set_nbnd = cast(Callable[[Workflow, int], None], partial(set_calculator_parameter, key='nbnd'))
_get_nbnd = cast(Callable[[Workflow], int], partial(get_calculator_parameter, key='nbnd'))


def conv_var_nbnd(increment: int = 1, **kwargs) -> ConvergenceVariable:

    return ConvergenceVariable('nbnd', increment, _get_nbnd, _set_nbnd, **kwargs)


def _get_kgrid(wf: Workflow) -> List[int]:
    assert wf.kpoints.grid
    return wf.kpoints.grid


def _set_kgrid(wf: Workflow, value: List[int]) -> None:
    wf.kpoints.grid = value


def _generate_kgrid_values(initial_value: List[int], increment: List[int], length: int) -> List[List[int]]:
    values: List[List[int]] = [initial_value]

    while len(values) < length:
        candidates = []
        for delta in itertools.product([0, 1], repeat=len(initial_value)):
            if 1 not in delta:
                continue
            candidates.append([x + d * i for x, d, i in zip(values[-1], delta, increment)])

        def cos_angle(c):
            return np.dot(c, initial_value) / np.linalg.norm(c) / np.linalg.norm(initial_value)
        best_candidate = sorted(candidates, key=cos_angle)[-1]
        values.append(best_candidate)

    return values


def conv_var_kgrid(increment: List[int] = [1, 1, 1], **kwargs) -> ConvergenceVariable:

    return ConvergenceVariable('kgrid', increment, _get_kgrid, _set_kgrid,
                               generate_values=_generate_kgrid_values, **kwargs)


def _get_num_training_snapshots(wf: Workflow) -> int:
    nts = wf.ml.number_of_training_snapshots
    assert isinstance(nts, int)
    return nts


def _set_num_training_snapshots(wf: Workflow, value: int) -> None:
    wf.ml.number_of_training_snapshots = value


def conv_var_number_of_training_snapshots(increment: int = 1, **kwargs) -> ConvergenceVariable:
    if 'length' not in kwargs:
        kwargs['length'] = 10
    return ConvergenceVariable('number_of_training_snapshots', increment, _get_num_training_snapshots,
                               _set_num_training_snapshots, **kwargs)


def ConvergenceVariableFactory(conv_var, **kwargs) -> ConvergenceVariable:

    if conv_var == 'celldm1':
        return conv_var_celldm1(**kwargs)
    elif conv_var == 'ecutwfc':
        return conv_var_ecutwfc(**kwargs)
    elif conv_var == 'nbnd':
        return conv_var_nbnd(**kwargs)
    elif conv_var == 'kgrid':
        return conv_var_kgrid(**kwargs)
    elif conv_var == 'number_of_training_snapshots':
        return conv_var_number_of_training_snapshots(**kwargs)
    else:
        raise NotImplementedError(f'Convergence with respect to `{conv_var}` has not been directly implemented. You '
                                  'can still perform a convergence calculation with respect to this variable, but '
                                  'you must first create an appropriate `ConvergenceVariable` object and then '
                                  'construct your `ConvergenceWorkflow` using the `ConvergenceWorkflowFactory`')


class ConvergenceOutputs(OutputModel):

    converged_values: Dict[str, Any]


class ConvergenceWorkflow(Workflow):

    output_model = ConvergenceOutputs  # type: ignore
    outputs: ConvergenceOutputs

    '''
    A Workflow class that wraps another workflow in a convergence procedure in order to converge the observable within the specified tolerance with respect to the variables
    '''

    def __init__(self, subworkflow_class: Type[Workflow], observable: Optional[Callable[[Workflow], float]] = None, threshold: Optional[float] = None, variables: List[ConvergenceVariable] = [], *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._subworkflow_class = subworkflow_class
        self.observable = observable
        self.threshold = threshold
        self.variables = variables

    @classmethod
    def fromdict(cls, dct: Dict[str, Any], **kwargs) -> Workflow:
        kwargs.update(subworkflow_class=dct.pop('_subworkflow_class'),
                      observable=dct.pop('observable'),
                      threshold=dct.pop('threshold'),
                      variables=dct.pop('variables'))
        return super(ConvergenceWorkflow, cls).fromdict(dct, **kwargs)

    def _run(self) -> None:
        # Set the initial value for each of the convergence variables
        for c in self.variables:
            if c.initial_value is None:
                c.get_initial_value(self)

        # Check that everything has been initialized
        if self.observable is None:
            raise ValueError(
                f'`{self.__class__.__name__}` has not been provided with an observable to converge')

        if len(self.variables) == 0:
            raise ValueError(
                f'`{self.__class__.__name__}` has not been provided with any variables with which to perform convergence')

        if self.threshold is None:
            raise ValueError(
                f'`{self.__class__.__name__}` has not been provided with a threshold with which to perform convergence')

        # Create array for storing calculation results
        results = np.empty([len(v) for v in self.variables])
        results[:] = np.nan

        # Continue to increment the convergence variables until convergence in the convergence observable
        # is achieved. A lot of the following code is quite obtuse, in order to make it work for an arbitrary
        # number of convergence variables

        # Record the alpha values for the original calculation
        provide_alpha = any([p.name == 'nbnd' for p in self.variables]) \
            and self.parameters.functional in ['ki', 'kipz', 'pkipz', 'all'] \
            and self.parameters.alpha_from_file
        if provide_alpha:
            master_alphas = utils.read_alpha_file(self)
            if self.parameters.orbital_groups is None:
                self.parameters.orbital_groups = list(
                    range(self.calculator_parameters['kcp'].nbnd))
            master_orbital_groups = copy.deepcopy(
                self.parameters.orbital_groups)

        while True:
            # Loop over all possible permutations of convergence variables
            subwfs = []
            indices_to_run = []
            for indices in itertools.product(*[range(len(x)) for x in self.variables]):
                # If we already have results for this permutation, don't recalculate it
                if not np.isnan(results[tuple(indices)]):
                    continue

                # Create a new subworkflow
                subwf = self._subworkflow_class.fromparent(self)

                # For each parameter we're converging wrt...
                label = ''
                for index, variable in zip(indices, self.variables):
                    value = variable.values[index]
                    if isinstance(value, float):
                        value_str = f'{value:.1f}'
                    else:
                        value_str = str(value)

                    if isinstance(value, list):
                        value_str = ''.join([str(x) for x in value])

                    # Create new working directory
                    label += f' {variable.name} {value_str}'.replace('.', '_')

                    # Set the value
                    variable.set_value(subwf, value)

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

                # Perform calculation
                subwf.name += label

                subwfs.append(subwf)
                indices_to_run.append(indices)

            for w in subwfs:
                w.run()

            if any([w.status != Status.COMPLETED for w in subwfs]):
                return

            # Store the result
            for indices, subwf in zip(indices_to_run, subwfs):
                results[indices] = self.observable(subwf)

            # Check convergence
            converged = np.empty([len(p) for p in self.variables], dtype=bool)
            converged[:] = False

            converged_result = results[tuple([-1 for _ in self.variables])]
            subarray_slice: List[Union[int, slice]]
            for indices in itertools.product(*[range(len(p)) for p in self.variables]):
                subarray_slice = [slice(None) for _ in range(len(indices))]
                for i_index, index in enumerate(indices):
                    subarray_slice[i_index] = slice(index, None)
                if np.all(np.abs(results[tuple(subarray_slice)] - converged_result) < self.threshold):
                    converged[tuple(indices)] = True

            # Check convergence, ignoring the calculations with the most fine variables because
            # we have nothing to compare them against
            subarray_slice = [slice(0, -1) for _ in self.variables]

            if np.any(converged[tuple(subarray_slice)]):
                # Work out which was the fastest calculation, and propose those variables

                # First, find the indices of the converged array
                slice_where: List[Union[slice, int]] = [
                    slice(None) for _ in self.variables]
                slice_where[-1] = 0
                converged_indices = np.array(np.where(converged[tuple(subarray_slice)]))[
                    tuple(slice_where)]

                # Extract the corresponding variables
                for index, variable in zip(converged_indices, self.variables):
                    variable.converged_value = variable.values[index]

                self.print('\n Converged variables are '
                           + ', '.join([f'{p.name} = {p.converged_value}' for p in self.variables]))

                self.outputs = self.output_model(converged_values={v.name: v.converged_value for v in self.variables})

                self.status = Status.COMPLETED

                return
            else:
                # Work out which variables are yet to converge, and line up more calculations
                # for increased values of those variables

                new_array_shape = list(np.shape(results))
                new_array_slice: List[Union[int, slice]] = [
                    slice(None) for _ in indices]
                self.print('\nProgress update', style='heading')
                for index, var in enumerate(self.variables):
                    subarray_slice = [slice(None) for _ in self.variables]
                    subarray_slice[index] = slice(0, -1)

                    if np.any(converged[tuple(subarray_slice)]):
                        self.print(f'{var.name} appears converged')
                        continue

                    # Add one finer value to try
                    var.extend()
                    new_array_shape[index] += 1
                    new_array_slice[index] = slice(None, -1)
                    self.print(
                        f'{var.name} still appears unconverged, will reattempt using a finer value')

                # Deal with the edge case where all variables appear converged for the finest value of all the other variables
                if results.shape == tuple([len(v) for v in self.variables]):
                    self.print(
                        '... but the variables are not collectively converged. Cautiously incrementing all variables.')
                    for index, var in enumerate(self.variables):
                        var.extend()
                        new_array_shape[index] += 1
                        new_array_slice[index] = slice(None, -1)
                self.print()

                new_results = np.empty(new_array_shape)
                new_results[:] = np.nan
                new_results[tuple(new_array_slice)] = results
                results = new_results


def ConvergenceWorkflowFactory(subworkflow: Workflow, observable: Union[str, Callable[[Workflow], float]], threshold: float,
                               variables: List[Union[str, ConvergenceVariable]]) -> ConvergenceWorkflow:
    '''
    Generates a ConvergenceWorkflow based on...
        subworkflow_class: the class of the Workflow that we will repeatedly run
        observable: a function with which we extract the observable from the child workflow
        threshold: the threshold within which the observable will be converged
        variables: a list of variables with respect to which we will perform the convergence
    '''

    # Ensure variables are all ConvergenceVariables
    variables = [v if isinstance(
        v, ConvergenceVariable) else ConvergenceVariableFactory(v) for v in variables]

    # Ensure observable is a Callable
    observable = ObservableFactory(observable) if isinstance(
        observable, str) else observable

    # Initialize a ConvergenceWorkflow copying the settings of the subworkflow
    wf = ConvergenceWorkflow.from_other(subworkflow, subworkflow_class=subworkflow.__class__,
                                        variables=variables, observable=observable, threshold=threshold)

    return wf
