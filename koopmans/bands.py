import itertools
from typing import Optional, List, Union, Any, Dict
from pathlib import Path
import numpy as np
import pandas as pd
from ase import Atoms
from ase.io.wannier90 import num_wann_from_projections, proj_string_to_dict, proj_dict_to_string
from koopmans.utils import indented_print, list_to_formatted_str


class Band(object):
    def __init__(self, index: Optional[int] = None, spin: int = 0, filled: bool = True, group: Optional[int] = None,
                 alpha: Optional[float] = None, error: Optional[float] = None,
                 self_hartree: Optional[float] = None,
                 centre: Optional[np.ndarray] = None) -> None:
        self.index = index
        self.spin = spin
        self.filled = filled
        self.group = group
        self.alpha_history: List[float] = []
        self.error_history: List[float] = []
        self.alpha = alpha
        self.error = error
        self.self_hartree = self_hartree
        self.centre = centre

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
    def __init__(self, n_bands: Union[int, List[int]], n_spin: int = 1, spin_polarised: bool = False, self_hartree_tol=None, **kwargs):
        if isinstance(n_bands, int):
            self.n_bands = [n_bands for _ in range(n_spin)]
        else:
            if len(n_bands) != n_spin:
                raise ValueError(f'n_bands = {n_bands} should have length matching n_spin = {n_spin}')
            self.n_bands = n_bands
        self.n_spin = n_spin
        self.spin_polarised = spin_polarised
        if self.spin_polarised:
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
        bands = dct['_bands']
        return cls(bands=bands, **dct)

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
        assert len(
            array) == self.n_spin, f'Bands.{array_name} must be length {self.n_spin} but you provided an array of length {len(array)}'
        for i, (subarray, n_bands) in enumerate(zip(array, self.n_bands)):
            assert len(
                subarray) == n_bands, f'Bands.{array_name}[{i}] must be length {n_bands} but you provided an array with length {len(subarray)}'

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

        # By default use the settings provided when Bands() was initialised
        sh_tol = sh_tol if sh_tol is not None else self.self_hartree_tol

        # Separate the orbitals into different subsets, where we don't want any grouping of orbitals belonging to
        # different subsets
        if self.spin_polarised:
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

                # Find the neighbourhood of adjacent bands (with 2x larger threshold)
                neighbourhood = [b for b in unassigned if points_are_close(guess, b)]

                # Find the centre of that neighbourhood
                av_sh = np.mean([b.self_hartree for b in neighbourhood])
                centre = Band(self_hartree=av_sh)

                # Find a revised neighbourhood close to the centre (using a factor of 0.5 because we want points that
                # are on opposite sides of the neighbourhood to be within "tol" of each other which means they can be
                # at most 0.5*tol away from the neighbourhood centre
                neighbourhood = [b for b in unassigned if points_are_close(centre, b, 0.5)]

                # Check the neighbourhood is isolated
                wider_neighbourhood = [b for b in unassigned if points_are_close(centre, b)]

                if neighbourhood != wider_neighbourhood:
                    if self.self_hartree_tol and sh_tol < 0.01 * self.self_hartree_tol:
                        # We have recursed too deeply, abort
                        raise Exception('Clustering algorithm failed')
                    else:
                        self.assign_groups(sh_tol=0.9 * sh_tol if sh_tol else None,
                                           allow_reassignment=allow_reassignment)
                        return

                for b in neighbourhood:
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

        if not self.spin_polarised and self.n_spin == 2:
            for b in self.get(spin=1):
                [match] = [b_op for b_op in self.get(spin=1) if b_op.index == b.index]
                b.group = match.group

        return

    @property
    def to_solve(self):
        # Dynamically generate a list of bands that require solving explicitly

        # If groups have not been assigned...
        if None in [b.group for b in self]:
            if self.spin_polarised:
                # ... and spin-polarised, solve all bands
                return self.get()
            else:
                # ... and not spin-polarised, solve the spin-up bands only
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
            if self.spin_polarised:
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
        if self.spin_polarised:
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


class WannierBandBlock(object):
    # This simple object contains th projections, filling, and spin corresponding to a block of bands
    def __init__(self,
                 projections: List[Union[str, Dict[str, Any]]],
                 filled: bool,
                 spin: Optional[int]):

        self.projections = []
        for proj in projections:
            if isinstance(proj, str):
                proj = proj_string_to_dict(proj)
            self.projections.append(proj)
        self.filled = filled
        self.spin = spin

    def __repr__(self) -> str:
        out = f'WannierBandBlock(projections={[self.projections]}, filled={self.filled}'
        if self.spin is not None:
            out += ', spin={self.spin}'
        return out + ')'

    def todict(self) -> dict:
        dct = {k: v for k, v in self.__dict__.items()}
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @property
    def w90_kwargs(self) -> Dict[str, Any]:
        # Returns the keywords to provide when constructing a new calculator corresponding to this block
        kwargs = {}
        for key in ['projections', 'num_wann', 'num_bands', 'exclude_bands']:
            val = getattr(self, key, None)
            if val is None:
                raise AttributeError(f'You must define {self.__class__.__name__}.{key} before requesting w90_kwargs')
            kwargs[key] = val
        if self.spin is not None:
            kwargs['spin_component'] = self.spin
        return kwargs

    @classmethod
    def fromdict(cls, dct):
        return cls(**dct)


class WannierBandBlocks(object):
    """
    This object is a collection of blocks of bands. In addition to the band blocks themselves, it also stores
    system-wide properties such as how many extra conduction bands we have.

    It allow easy iteration over all the blocks with the syntax

    > for block in bandblocks:
    >     ...

    where in this case "block" is a dictionary that is populated with all the relevant keywords pertaining
    to that band block, such as spin (which could alternatively be extracted from the individual band block) but also
    keywords such as exclude_bands (which require knowledge of the total number of bands and wannier functions to
    construct). See self.__iter__() for more details.
    """

    def __init__(self, blocks: List[WannierBandBlock], atoms: Atoms):
        self._blocks = blocks
        self._atoms = atoms
        # This BandBlocks object must keep track of how many bands we have not belonging to any block
        self._n_bands_below = {spin: 0 for spin in set([b.spin for b in blocks])}
        self._n_bands_above = {spin: 0 for spin in set([b.spin for b in blocks])}

    def __repr__(self):
        out = 'WannierBandBlocks('
        for b in self._blocks:
            out += f'{b}\n                  '
        out = out.strip()
        out += ')'
        return out

    @property
    def blocks(self):
        # Before returning all the blocks, add more global information such as exclude_bands
        for spin in [None, 'up', 'down']:
            subset = self.get_subset(spin=spin)
            if len(subset) == 0:
                continue
            is_last = [False for _ in subset]
            is_last[-1] = True
            wann_counter = self._n_bands_below[spin] + 1
            for b, include_above in zip(subset, is_last):
                # Construct num_wann
                b.num_wann = num_wann_from_projections(b.projections, self._atoms)

                # Construct num_bands
                b.num_bands = b.num_wann
                if include_above:
                    b.num_bands += self._n_bands_above[spin]

                # Construct exclude_bands
                band_indices = range(wann_counter, wann_counter + b.num_wann)
                wann_counter += b.num_wann
                if include_above:
                    # For the uppermost block we don't want to exclude any extra bands from the wannierisation
                    upper_bound = max(band_indices)
                else:
                    upper_bound = self.num_bands(spin=spin)
                to_exclude = [i for i in range(1, upper_bound + 1) if i not in band_indices]
                b.exclude_bands = list_to_formatted_str(to_exclude)

                # Construct the calc_type
                if b.filled:
                    label = 'occ'
                else:
                    label = 'emp'
                if b.spin is not None:
                    label += '_' + spin
                b.calc_type = 'w90_' + label

                # Construct directory and info for merging
                b.directory = label

                subset = self.get_subset(occ=b.filled, spin=b.spin)
                if len(subset) > 1:
                    b.to_merge = True
                    b.merge_directory = b.directory
                    b.directory = label + f'_block{subset.index(b) + 1}'
                else:
                    b.to_merge = False
        return self._blocks

    def __iter__(self):
        for b in self.blocks:
            yield b

    @classmethod
    def fromprojections(cls,
                        list_of_projections: List[List[Union[str, Dict[str, Any]]]],
                        fillings: List[bool],
                        spins: List[Union[str, None]],
                        atoms: Atoms):

        # Make sure to store all filled blocks before any empty blocks
        blocks: List[WannierBandBlock] = []
        for filled in [True, False]:
            blocks += [WannierBandBlock(p, f, s)
                       for p, f, s in zip(list_of_projections, fillings, spins) if f is filled]

        return cls(blocks, atoms)

    def add_bands(self, num: int, above: bool = False, spin: Optional[str] = None):
        if above:
            self._n_bands_above[spin] += num
        else:
            self._n_bands_below[spin] += num

    def get_subset(self, occ: Optional[bool] = None, spin: Optional[str] = None):
        return [b for b in self._blocks if (occ is None or b.filled == occ) and (spin is None or b.spin == spin)]

    def num_wann(self, occ: Optional[bool] = None, spin: Optional[str] = None):
        return sum([num_wann_from_projections(b.projections, self._atoms) for b in self.get_subset(occ, spin)])

    def num_bands(self, occ: Optional[bool] = None, spin: Optional[str] = None):
        nbands = self.num_wann(occ, spin)
        try:
            if (occ is True or occ is None):
                nbands += self._n_bands_below[spin]
            if (occ is False or occ is None):
                nbands += self._n_bands_above[spin]
        except KeyError:
            raise ValueError('num_bands does not support summing over all spins')
        return nbands

    def todict(self) -> dict:
        dct: Dict[str, Any] = {k: v for k, v in self.__dict__.items()}
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @ classmethod
    def fromdict(cls, dct):
        new_bandblock = cls(dct.pop('_blocks'), dct.pop('_atoms'))
        for k, v in dct.items():
            setattr(new_bandblock, k, v)
        return new_bandblock
