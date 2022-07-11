from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from ase import Atoms
from ase.io.wannier90 import (list_to_formatted_str, num_wann_from_projections,
                              proj_string_to_dict)


class ProjectionBlock(object):
    # This simple object contains the projections, filling, and spin corresponding to a block of bands
    def __init__(self,
                 projections: List[Union[str, Dict[str, Any]]],
                 filled: bool,
                 spin: Optional[str] = None,
                 directory: Optional[Path] = None,
                 merge_directory: Optional[Path] = None,
                 num_wann: Optional[int] = None,
                 num_bands: Optional[int] = None,
                 exclude_bands: Optional[str] = None,
                 calc_type: Optional[str] = None,
                 to_merge: Optional[bool] = None):

        self.projections = []
        for proj in projections:
            if isinstance(proj, str):
                proj = proj_string_to_dict(proj)
            self.projections.append(proj)
        self.filled = filled
        self.spin = spin
        self.directory = directory
        self.merge_directory = merge_directory
        self.num_wann = num_wann
        self.num_bands = num_bands
        self.exclude_bands = exclude_bands
        self.calc_type = calc_type
        self.to_merge = to_merge

    def __repr__(self) -> str:
        out = f'ProjectionBlock({[self.projections]}, filled={self.filled}'
        if self.spin is not None:
            out += f', spin={self.spin}'
        return out + ')'

    def __len__(self):
        # Count the number of projections
        return len(self.projections)

    def __bool__(self):
        # Return true if this ProjectionBlock is non-empty
        return len(self) > 0

    def __eq__(self, other):
        if isinstance(other, ProjectionBlock):
            return self.__dict__ == other.__dict__
        return False

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
            if val is None and key != 'exclude_bands':
                raise AttributeError(f'You must define {self.__class__.__name__}.{key} before requesting w90_kwargs')
            kwargs[key] = val
        if self.spin is not None:
            kwargs['spin'] = self.spin
        return kwargs

    @classmethod
    def fromdict(cls, dct):
        return cls(**dct)


class ProjectionBlocks(object):
    """
    This object is a collection of blocks of projections. In addition to the projections blocks themselves, it also
    stores system-wide properties such as how many extra conduction bands we have.

    Whenever a user queries self.blocks (e.g. when they iterate over this object) it will first propagate these
    system-wide properties down to the individual ProjectionBlock objects. See self.blocks() for more details.
    """

    def __init__(self, blocks: List[ProjectionBlock], atoms: Atoms):
        self._blocks = blocks
        self._atoms = atoms
        # This BandBlocks object must keep track of how many bands we have not belonging to any block
        if len(blocks) > 0:
            self._n_bands_below = {spin: 0 for spin in set([b.spin for b in blocks])}
            self._n_bands_above = {spin: 0 for spin in set([b.spin for b in blocks])}
        else:
            # By default, assume spin-unpolarized
            self._n_bands_below = {None: 0}
            self._n_bands_above = {None: 0}

    def __repr__(self):
        out = 'ProjectionBlocks('
        for b in self._blocks:
            out += f'{b}\n                  '
        out = out.strip()
        out += ')'
        return out

    def __len__(self):
        # Count the number of non-empty "ProjectionBlock"s
        return sum([1 for b in self._blocks if b])

    def __bool__(self):
        # Return true if this contains any non-empty "ProjectionBlock"s
        return len(self) > 0

    def __eq__(self, other):
        if isinstance(other, ProjectionBlocks):
            return self.__dict__ == other.__dict__
        return False

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
                    # For the uppermost block we don't want to exclude any extra bands from the wannierization
                    upper_bound = max(band_indices)
                else:
                    upper_bound = self.num_bands(spin=spin)
                to_exclude = [i for i in range(1, upper_bound + 1) if i not in band_indices]
                if len(to_exclude) > 0:
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
                b.directory = Path(label)

                subset = self.get_subset(occ=b.filled, spin=b.spin)
                if len(subset) > 1:
                    b.to_merge = True
                    b.merge_directory = b.directory
                    b.directory = Path(label + f'_block{subset.index(b) + 1}')
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
        blocks: List[ProjectionBlock] = []
        for filled in [True, False]:
            blocks += [ProjectionBlock(p, f, s)
                       for p, f, s in zip(list_of_projections, fillings, spins) if f is filled and len(p) > 0]

        return cls(blocks, atoms)

    def add_bands(self, num: int, above: bool = False, spin: Optional[str] = None):
        if num < 0:
            raise ValueError('num must be > 0')
        if above:
            self._n_bands_above[spin] += num
        else:
            self._n_bands_below[spin] += num

    def to_merge(self):
        for occ in [True, False]:
            for spin in [None, 'up', 'down']:
                subset = self.get_subset(occ=occ, spin=spin)
                if len(subset) > 1:
                    yield subset

    def get_subset(self, occ: Optional[bool] = None, spin: Optional[str] = 'both'):
        return [b for b in self._blocks if (occ is None or b.filled == occ) and (spin == 'both' or b.spin == spin)]

    def num_wann(self, occ: Optional[bool] = None, spin: Optional[str] = 'both'):
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

    @classmethod
    def fromdict(cls, dct):
        new_bandblock = cls(dct.pop('_blocks'), dct.pop('_atoms'))
        for k, v in dct.items():
            if not hasattr(new_bandblock, k):
                raise AttributeError(k)
            setattr(new_bandblock, k, v)
        return new_bandblock
