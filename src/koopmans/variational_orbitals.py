import itertools
from collections import UserList
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from pydantic import BaseModel, Field

from koopmans.files import File
from koopmans.utils import indented_print, warn


class VariationalOrbital(BaseModel):

    index: int
    spin: int = 0
    filled: bool = True
    group: int
    alpha_history: List[float] = Field(default_factory=list)
    error_history: List[float] = Field(default_factory=list)
    predicted_alpha: Optional[float] = None
    self_hartree: Optional[float] = None
    spread: Optional[float] = None
    center: Optional[np.ndarray] = None
    power_spectrum: Optional[File] = None
    model_config = {'arbitrary_types_allowed': True}

    @classmethod
    def fromdict(cls, dct):
        alpha_history = dct.pop('alpha_history')
        error_history = dct.pop('error_history')
        var_orb = cls(**dct)
        var_orb.alpha_history = alpha_history
        var_orb.error_history = error_history
        return var_orb

    def todict(self) -> dict:
        dct = self.__dict__
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    def __eq__(self, other):
        if not isinstance(other, VariationalOrbital):
            return False
        return self.__dict__ == other.__dict__

    def __repr__(self) -> str:
        info = f'{self.__class__.__name__}(index={self.index}, spin={self.spin}, filled={self.filled}, group={self.group}'
        for attr in ['alpha', 'self_hartree', 'spread', 'center']:
            val = getattr(self, attr, None)
            if val:
                info += f', {attr.replace("_", "-")}={val}'
        return info + ')'

    @property
    def alpha(self) -> float:
        if len(self.alpha_history) == 0:
            raise AttributeError('VariationalOrbital does not have screening parameters')
        return self.alpha_history[-1]

    @alpha.setter
    def alpha(self, value: Optional[float]):
        if value is not None:
            assert isinstance(value, float)
            self.alpha_history.append(value)

    @property
    def error(self):
        if len(self.error_history) == 0:
            raise AttributeError('VariationalOrbital does not have error data')
        return self.error_history[-1]

    @error.setter
    def error(self, value: Optional[float]):
        if value is not None:
            assert isinstance(value, float)
            self.error_history.append(value)


class VariationalOrbitals(BaseModel):

    orbitals: List[VariationalOrbital] = Field(default_factory=list)
    n_spin: int = 1
    spin_polarized: bool = False
    tolerances: Dict[str, float] = Field(default_factory=dict)

    def __iter__(self):
        for o in self.orbitals:
            yield o

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}([\n  ' + '\n  '.join([str(o) for o in self]) + '\n])'

    def __len__(self) -> int:
        return len(self.orbitals)

    @classmethod
    def fromlist(cls, orbitals: List[VariationalOrbital]):
        raise NotImplementedError('TODO')
        # return orbitals

    def get(self, spin: Optional[int] = None, filled: Optional[bool] = None, group: Optional[int] = None,
            to_solve: Optional[bool] = None) -> List[VariationalOrbital]:
        if to_solve:
            selected_orbs = self.to_solve
        else:
            selected_orbs = self
        if filled is not None:
            selected_orbs = [o for o in selected_orbs if o.filled == filled]
        if group is not None:
            selected_orbs = [o for o in selected_orbs if o.group == group]
        if spin is not None:
            selected_orbs = [o for o in selected_orbs if o.spin == spin]

        return selected_orbs

    def __getitem__(self, key):
        return self.orbitals[key]

    def num(self, filled=None, spin=None):
        return len(self.get(filled=filled, spin=spin))

    def index(self, orbital: VariationalOrbital) -> int:
        if orbital not in self.orbitals:
            raise ValueError(f"{orbital} is not in {self}")
        [i_match] = [i for i, b in enumerate(self.orbitals) if b == orbital]
        return i_match

    # def _check_array_shape_match(self, array, array_name):
    #     raise NotImplementedError()
    #     assert len(array) == self.n_spin, f'Bands.{array_name} must be length {self.n_spin} but you provided an ' \
    #         f'array of length {len(array)}'
    #     for i, (subarray, n_bands) in enumerate(zip(array, self.n_bands)):
    #         assert len(subarray) == n_bands, f'Bands.{array_name}[{i}] must be length {n_bands} but you provided an ' \
    #             f'array with length {len(subarray)}. The file_alpharef files should reflect the number of states in the supercell.'

    @property
    def filling(self) -> List[List[bool]]:
        return [[o.filled for o in self if o.spin == i_spin] for i_spin in range(self.n_spin)]

    @filling.setter
    def filling(self, value: List[List[bool]]):
        self._check_array_shape_match(value, 'filling')
        for o, v in zip(self, [v for subarray in value for v in subarray]):
            o.filled = v

    @property
    def indices(self):
        return [[o.index for o in self if o.spin == i_spin] for i_spin in range(self.n_spin)]

    @property
    def groups(self):
        return [[o.group for o in self if o.spin == i_spin] for i_spin in range(self.n_spin)]

    @groups.setter
    def groups(self, value: List[List[int]]):
        self._check_array_shape_match(value, 'groups')
        for o, v in zip(self, [v for subarray in value for v in subarray]):
            o.group = v

    def assign_groups(self, sort_by: str = 'self_hartree', tol: Optional[float] = None, allow_reassignment: bool = False):
        # Basic clustering algorithm for assigning groups

        if self.tolerances == {}:
            # Do not perform clustering
            return

        if sort_by not in self.tolerances:
            return ValueError(f'Cannot sort bands according to {sort_by}; valid choices are' + '/'.join(self.tolerances.keys()))

        # By default use the settings provided when Bands() was initialized
        tol = tol if tol is not None else self.tolerances[sort_by]

        # Separate the orbitals into different subsets, where we don't want any grouping of orbitals belonging to
        # different subsets
        if self.spin_polarized:
            # Separate by both spin and filling
            unassigned_sets = [[o for o in self if o.filled == filled and o.spin == i_spin]
                               for i_spin in range(self.n_spin) for filled in [True, False]]
        else:
            # Separate by filling and focus only on the spin=0 channel
            unassigned_sets = [[b for b in self if b.filled == filled and b.spin == 0] for filled in [True, False]]

        def points_are_close(p0: VariationalOrbital, p1: VariationalOrbital, factor: Union[int, float] = 1) -> bool:
            # Determine if two bands are "close"
            assert tol is not None
            return abs(getattr(p0, sort_by) - getattr(p1, sort_by)) < tol * factor

        group = 0
        for unassigned in unassigned_sets:
            while len(unassigned) > 0:
                # Select one band
                guess = unassigned[0]

                # Find the neighborhood of adjacent orbitals (with 2x larger threshold)
                neighborhood = [b for b in unassigned if points_are_close(guess, b)]

                # Find the center of that neighborhood
                av = np.mean([getattr(b, sort_by) for b in neighborhood])
                center = VariationalOrbital(**{sort_by: av})

                # Find a revised neighborhood close to the center (using a factor of 0.5 because we want points that
                # are on opposite sides of the neighborhood to be within "tol" of each other which means they can be
                # at most 0.5*tol away from the neighborhood center
                neighborhood = [b for b in unassigned if points_are_close(center, b, 0.5)]

                # Check the neighborhood is isolated
                wider_neighborhood = [b for b in unassigned if points_are_close(center, b)]

                if neighborhood != wider_neighborhood:
                    if self.tolerances[sort_by] and tol < 0.01 * self.tolerances[sort_by]:
                        # We have recursed too deeply, abort
                        raise Exception('Clustering algorithm failed')
                    else:
                        self.assign_groups(sort_by, tol=0.9 * tol if tol else None,
                                           allow_reassignment=allow_reassignment)
                        return

                for o in neighborhood:
                    unassigned.remove(o)

                    if allow_reassignment:
                        # Perform the reassignment
                        o.group = group
                    else:
                        # Check previous values exist
                        if o.group is None:
                            o.group = group
                        # Check the new grouping matches the old grouping
                        if o.group != group:
                            raise Exception('Clustering algorithm found different grouping')

                # Move on to next group
                group += 1

        if not self.spin_polarized and self.n_spin == 2:
            for o in self.get(spin=1):
                [match] = [b_op for b_op in self.get(spin=0) if b_op.index == o.index]
                o.group = match.group

        if tol != self.tolerances[sort_by]:
            warn(f'It was not possible to group orbitals with the {sort_by} tolerance of {self.tolerances[sort_by]:.2e} eV. '
                 f'A grouping was found for a tolerance of {tol:.2e} eV.\n'
                 f'Try a larger tolerance to group more orbitals together')

        return

    @property
    def to_solve(self):
        # Dynamically generate a list of orbitals that require solving explicitly

        # If groups have not been assigned...
        if None in [o.group for o in self]:
            if self.spin_polarized:
                # ... and spin-polarized, solve all bands
                return self.get()
            else:
                # ... and not spin-polarized, solve the spin-up bands only
                return self.get(spin=0)

        # If not, work out which bands to solve explicitly
        groups_found = set([])
        to_solve = []

        for orb in [o for i_spin in range(self.n_spin) for o in self.get(spin=i_spin)[::-1] if o.filled] \
                + [o for i_spin in range(self.n_spin) for o in self.get(spin=i_spin) if not o.filled]:
            # Looping through the filled bands from highest to lowest (first high spin then low spin), then empty
            # bands from lowest to highest (but note since these are variational orbitals "highest" and "lowest")
            # is not especially meaningful, unless we happen to be using KS orbitals as variational orbitals...
            if orb.group not in groups_found:
                groups_found.add(orb.group)
                to_solve.append(orb)

        if groups_found != set([b.group for b in self]):
            raise ValueError('Splitting of orbitals into groups failed')

        return sorted(to_solve, key=lambda x: (x.spin, x.index))

    @property
    def self_hartrees(self) -> List[float]:
        return [b.self_hartree for b in self]

    @self_hartrees.setter
    def self_hartrees(self, value: List[List[float]]) -> None:
        self._check_array_shape_match(value, 'self_hartrees')
        for b, v in zip(self, [v for subarray in value for v in subarray]):
            b.self_hartree = v

    @property
    def predicted_alphas(self) -> List[List[float]]:
        return [[b.predicted_alpha for b in self if b.spin == i_spin] for i_spin in range(self.n_spin)]

    @property
    def power_spectrum(self) -> List[List[float]]:
        return [b.power_spectrum for b in self]

    def update_attrib_with_history(self, name: str, value: Union[float, List[List[float]], np.ndarray, pd.DataFrame],
                                   group=None) -> None:
        '''
        Generic function for setting the band's screening parameters/errors to the value provided
         - "value" can be a scalar, a list, or a pandas DataFrame of the alpha_history
         - if "group" is provided then it applies this value to the orbitals belonging to this group only
        '''

        if isinstance(value, pd.DataFrame):
            assert group is None, 'Cannot update only one group via a pandas DataFrame'
            if self.spin_polarized:
                raise NotImplementedError()
            else:
                tmp_arr = np.transpose(value.values)
                array = [tmp_arr for _ in range(self.n_spin)]

            for spin, s_array in enumerate(array):
                for b, history in zip(self.get(spin=spin), s_array):
                    # Make sure to exclude NaNs
                    setattr(b, f'{name}_history', [a for a in history.tolist() if not np.isnan(a)])
            return

        if isinstance(value, float):
            value = [[value for _ in range(n_bands_spin)] for n_bands_spin in self.n_bands]

        self._check_array_shape_match(value, name)
        for b, v in zip(self, [v for subarray in value for v in subarray]):
            if group:
                if b.group != group:
                    continue
            setattr(b, name, v)

    @property
    def alphas(self):
        # This returns the alpha values for the iteration number where we have alpha for all bands
        i = min([len(b.alpha_history) for b in self]) - 1
        if i == -1:
            raise AttributeError()
        return [[b.alpha_history[i] for b in self if b.spin == i_spin] for i_spin in range(self.n_spin)]

    @alphas.setter
    def alphas(self, value):
        self.update_attrib_with_history('alpha', value)

    @property
    def errors(self):
        return [b.error for b in self]

    @errors.setter
    def errors(self, value):
        self.update_errors(value)

    def update_errors(self, value, group=None):
        self.update_attrib_with_history('error', value)

    def _create_dataframe(self, attr, spin=None, only_to_solve=True) -> pd.DataFrame:
        # Generate a dataframe containing the requested attribute, sorting the bands first by index, then by spin
        if self.spin_polarized and spin is None:
            if only_to_solve:
                blist = self.to_solve
            else:
                blist = self

            columns = pd.MultiIndex.from_tuples([(f'spin {b.spin}', b.index) for b in blist])
            band_subset = sorted(blist, key=lambda x: (x.spin, x.index))
        else:
            spin = 0 if spin is None else spin
            columns = [b.index for b in self.get(spin=spin, to_solve=only_to_solve)]
            band_subset = self.get(spin=spin, to_solve=only_to_solve)

        if isinstance(getattr(band_subset[0], attr), list):
            # Create an array of values padded with NaNs
            arr = np.array(list(itertools.zip_longest(*[getattr(b, attr) for b in band_subset], fillvalue=np.nan)))
        else:
            arr = np.array([[getattr(b, attr) for b in band_subset]])
        if arr.size == 0:
            df = pd.DataFrame(columns=columns)
        else:
            df = pd.DataFrame(arr, columns=columns)
        return df

    def alpha_history(self, spin=None) -> pd.DataFrame:
        return self._create_dataframe('alpha_history', spin=spin)

    def error_history(self, spin=None) -> pd.DataFrame:
        return self._create_dataframe('error_history', spin=spin)

    def predicted_alpha_history(self, spin=None) -> pd.DataFrame:
        return self._create_dataframe('predicted_alpha', spin=spin)
