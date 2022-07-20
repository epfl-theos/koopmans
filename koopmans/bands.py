import itertools
from typing import List, Optional, Union

import numpy as np
import pandas as pd

from koopmans.utils import indented_print, warn


class Band(object):
    def __init__(self, index: Optional[int] = None, spin: int = 0, filled: bool = True, group: Optional[int] = None,
                 alpha: Optional[float] = None, error: Optional[float] = None,
                 self_hartree: Optional[float] = None,
                 center: Optional[np.ndarray] = None) -> None:
        self.index = index
        self.spin = spin
        self.filled = filled
        self.group = group
        self.alpha_history: List[float] = []
        self.error_history: List[float] = []
        self.alpha = alpha
        self.error = error
        self.self_hartree = self_hartree
        self.center = center

    @classmethod
    def fromdict(cls, dct):
        alpha_history = dct.pop('alpha_history')
        error_history = dct.pop('error_history')
        band = cls(**dct)
        band.alpha_history = alpha_history
        band.error_history = error_history
        return band

    def todict(self) -> dict:
        dct = self.__dict__
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    def __repr__(self) -> str:
        return f'Band(index={self.index}, spin={self.spin}, filled={self.filled}, group={self.group})'

    @property
    def alpha(self) -> Union[float, None]:
        assert len(self.alpha_history) > 0, 'Band does not have screening parameters'
        return self.alpha_history[-1]

    @alpha.setter
    def alpha(self, value: Optional[float]):
        if value is not None:
            assert isinstance(value, float)
            self.alpha_history.append(value)

    @property
    def error(self):
        assert len(self.error_history) > 0, 'Band does not have error data'
        return self.error_history[-1]

    @error.setter
    def error(self, value: Optional[float]):
        if value is not None:
            assert isinstance(value, float)
            self.error_history.append(value)


class Bands(object):
    def __init__(self, n_bands: Union[int, List[int]], n_spin: int = 1, spin_polarized: bool = False,
                 self_hartree_tol=None, **kwargs):
        if isinstance(n_bands, int):
            self.n_bands = [n_bands for _ in range(n_spin)]
        else:
            if len(n_bands) != n_spin:
                raise ValueError(f'n_bands = {n_bands} should have length matching n_spin = {n_spin}')
            self.n_bands = n_bands
        self.n_spin = n_spin
        self.spin_polarized = spin_polarized
        if self.spin_polarized:
            # Assign every single band a distinct group
            self._bands: List[Band] = []
            for i_spin, n_bands_spin in enumerate(self.n_bands):
                for i_band in range(n_bands_spin):
                    self._bands.append(Band(i_band + 1, i_spin, group=len(self._bands)))
        else:
            # Assign bands with the same index but opposite spin the same group
            self._bands = [Band(i_band + 1, i_spin, group=i_band) for i_spin in range(n_spin)
                           for i_band in range(self.n_bands[i_spin])]

        self.self_hartree_tol = self_hartree_tol
        for k, v in kwargs.items():
            assert hasattr(self, k)
            if v:
                setattr(self, k, v)

    def __iter__(self):
        for b in self._bands:
            yield b

    @classmethod
    def fromdict(cls, dct):
        bands = dct.pop('_bands')
        obj = cls(**dct)
        obj._bands = bands
        return obj

    @classmethod
    def fromlist(cls, bands: List[Band]):
        raise NotImplementedError('TODO')
        # return bands

    def todict(self):
        dct = self.__dict__
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    def get(self, spin: Optional[int] = None, filled: Optional[bool] = None, group: Optional[int] = None,
            to_solve: Optional[bool] = None) -> List[Band]:
        if to_solve:
            selected_bands = self.to_solve
        else:
            selected_bands = self
        if filled is not None:
            selected_bands = [b for b in selected_bands if b.filled == filled]
        if group is not None:
            selected_bands = [b for b in selected_bands if b.group == group]
        if spin is not None:
            selected_bands = [b for b in selected_bands if b.spin == spin]

        return selected_bands

    def __getitem__(self, key):
        return self._bands[key]

    def num(self, filled=None, spin=None):
        return len(self.get(filled=filled, spin=spin))

    def index(self, band: Band) -> int:
        if band not in self._bands:
            raise ValueError(f"{band} is not in Bands object")
        [i_match] = [i for i, b in enumerate(self._bands) if b == band]
        return i_match

    def _check_array_shape_match(self, array, array_name):
        assert len(array) == self.n_spin, f'Bands.{array_name} must be length {self.n_spin} but you provided an ' \
            f'array of length {len(array)}'
        for i, (subarray, n_bands) in enumerate(zip(array, self.n_bands)):
            assert len(subarray) == n_bands, f'Bands.{array_name}[{i}] must be length {n_bands} but you provided an ' \
                f'array with length {len(subarray)}'

    @property
    def filling(self) -> List[List[bool]]:
        return [[b.filled for b in self if b.spin == i_spin] for i_spin in range(self.n_spin)]

    @filling.setter
    def filling(self, value: List[List[bool]]):
        self._check_array_shape_match(value, 'filling')
        for b, v in zip(self, [v for subarray in value for v in subarray]):
            b.filled = v

    @property
    def indices(self):
        return [[b.index for b in self if b.spin == i_spin] for i_spin in range(self.n_spin)]

    @property
    def groups(self):
        return [[b.group for b in self if b.spin == i_spin] for i_spin in range(self.n_spin)]

    @groups.setter
    def groups(self, value: List[List[int]]):
        self._check_array_shape_match(value, 'groups')
        for b, v in zip(self, [v for subarray in value for v in subarray]):
            b.group = v

    def assign_groups(self, sh_tol: Optional[float] = None, allow_reassignment: bool = False):
        # Basic clustering algorithm for assigning groups

        if self.self_hartree_tol is None:
            # Do not perform clustering
            return

        # By default use the settings provided when Bands() was initialized
        sh_tol = sh_tol if sh_tol is not None else self.self_hartree_tol

        # Separate the orbitals into different subsets, where we don't want any grouping of orbitals belonging to
        # different subsets
        if self.spin_polarized:
            # Separate by both spin and filling
            unassigned_sets = [[b for b in self if b.filled == filled and b.spin == i_spin]
                               for i_spin in range(self.n_spin) for filled in [True, False]]
        else:
            # Separate by filling and focus only on the spin=0 channel
            unassigned_sets = [[b for b in self if b.filled == filled and b.spin == 0] for filled in [True, False]]

        def points_are_close(p0: Band, p1: Band, factor: Union[int, float] = 1) -> bool:
            # Determine if two bands are "close"
            for obs, tol in (('self_hartree', sh_tol),):
                if tol is not None and abs(getattr(p0, obs) - getattr(p1, obs)) > tol * factor:
                    return False
            return True

        group = 0
        for unassigned in unassigned_sets:
            while len(unassigned) > 0:
                # Select one band
                guess = unassigned[0]

                # Find the neighborhood of adjacent bands (with 2x larger threshold)
                neighborhood = [b for b in unassigned if points_are_close(guess, b)]

                # Find the center of that neighborhood
                av_sh = np.mean([b.self_hartree for b in neighborhood])
                center = Band(self_hartree=av_sh)

                # Find a revised neighborhood close to the center (using a factor of 0.5 because we want points that
                # are on opposite sides of the neighborhood to be within "tol" of each other which means they can be
                # at most 0.5*tol away from the neighborhood center
                neighborhood = [b for b in unassigned if points_are_close(center, b, 0.5)]

                # Check the neighborhood is isolated
                wider_neighborhood = [b for b in unassigned if points_are_close(center, b)]

                if neighborhood != wider_neighborhood:
                    if self.self_hartree_tol and sh_tol < 0.01 * self.self_hartree_tol:
                        # We have recursed too deeply, abort
                        raise Exception('Clustering algorithm failed')
                    else:
                        self.assign_groups(sh_tol=0.9 * sh_tol if sh_tol else None,
                                           allow_reassignment=allow_reassignment)
                        return

                for b in neighborhood:
                    unassigned.remove(b)

                    if allow_reassignment:
                        # Perform the reassignment
                        b.group = group
                    else:
                        # Check previous values exist
                        if b.group is None:
                            b.group = group
                        # Check the new grouping matches the old grouping
                        if b.group != group:
                            raise Exception('Clustering algorithm found different grouping')

                # Move on to next group
                group += 1

        if not self.spin_polarized and self.n_spin == 2:
            for b in self.get(spin=1):
                [match] = [b_op for b_op in self.get(spin=0) if b_op.index == b.index]
                b.group = match.group

        if sh_tol != self.self_hartree_tol:
            warn(f'It was not possible to group orbitals with the self-Hartree tolerance of {self.self_hartree_tol:.2e} eV. '
                 f'A grouping was found for a self-Hartree tolerance of {sh_tol:.2e} eV.\n'
                 f'Try a larger self-Hartree tolerance to group more orbitals together')

        return

    @property
    def to_solve(self):
        # Dynamically generate a list of bands that require solving explicitly

        # If groups have not been assigned...
        if None in [b.group for b in self]:
            if self.spin_polarized:
                # ... and spin-polarized, solve all bands
                return self.get()
            else:
                # ... and not spin-polarized, solve the spin-up bands only
                return self.get(spin=0)

        # If not, work out which bands to solve explicitly
        groups_found = set([])
        to_solve = []

        for band in [b for i_spin in range(self.n_spin) for b in self.get(spin=i_spin)[::-1] if b.filled] \
                + [b for i_spin in range(self.n_spin) for b in self.get(spin=i_spin) if not b.filled]:
            # Looping through the filled bands from highest to lowest (first high spin then low spin), then empty
            # bands from lowest to highest (but note since these are variational orbitals "highest" and "lowest")
            # is not especially meaningful, unless we happen to be using KS orbitals as variational orbitals...
            if band.group not in groups_found:
                groups_found.add(band.group)
                to_solve.append(band)

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

    def _create_dataframe(self, attr) -> pd.DataFrame:
        # Generate a dataframe containing the requested attribute, sorting the bands first by index, then by spin
        if self.spin_polarized:
            columns = pd.MultiIndex.from_tuples([(f'spin {b.spin}', b.index) for b in self])
            band_subset = sorted(self, key=lambda x: (x.spin, x.index))
        else:
            columns = [b.index for b in self.get(spin=0)]
            band_subset = self.get(spin=0)

        # Create an array of values padded with NaNs
        arr = np.array(list(itertools.zip_longest(*[getattr(b, attr) for b in band_subset], fillvalue=np.nan)))
        if arr.size == 0:
            df = pd.DataFrame(columns=columns)
        else:
            df = pd.DataFrame(arr, columns=columns)
        return df

    @property
    def alpha_history(self) -> List[pd.DataFrame]:
        return self._create_dataframe('alpha_history')

    @property
    def error_history(self):
        return self._create_dataframe('error_history')

    def print_history(self, indent: int = 0):
        # Printing out a progress summary
        indented_print(f'\nalpha', indent=indent)
        indented_print(str(self.alpha_history), indent=indent)
        if not self.error_history.empty:
            indented_print(f'\nDelta E_i - epsilon_i (eV)', indent=indent)
            indented_print(str(self.error_history), indent=indent)
        indented_print('')
