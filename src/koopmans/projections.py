"""Classes that define blocks of bands."""

from typing import Any, Dict, List, Optional, Union

from ase_koopmans import Atoms
from ase_koopmans.io.wannier90 import (list_to_formatted_str,
                                       num_wann_from_projections,
                                       proj_string_to_dict)
from pydantic import BaseModel, ConfigDict, model_validator

from koopmans import calculators
from koopmans.utils import SpinOptions, SpinType


class BlockID(BaseModel):
    """The ID of a block of bands, which includes a label, filling, and spin."""

    label: Optional[str] = None
    filled: Optional[bool] = None
    spin: SpinType = None
    model_config = ConfigDict(frozen=True)

    def __str__(self):
        if self.label is None:
            return ''
        return f'{self.label}_{self.spin}' if self.spin else self.label

    # Set label to 'occ' or 'emp' if filled is set but label is not
    @model_validator(mode='before')
    @classmethod
    def set_label_from_filled(cls, data: Any):
        """Construct a label (if it is not provided) from the filling and spin arguments."""
        if isinstance(data, dict):
            if data.get('filled', None) is not None and data.get('label', None) is None:
                spin = data.get('spin', None)
                spin_label = '_spin_' + spin if spin is not None else ''
                data['label'] = 'occ' + spin_label if data['filled'] else 'emp' + spin_label
        return data

    def filling(self):
        """Return a string indicating whether the block is occupied or empty."""
        if self.filled:
            return 'occ'
        elif self.filled is False:
            return 'emp'
        else:
            raise ValueError(f'The filling of block {str(self)} is not known')


class ProjectionBlock(object):
    """A class that contains the projections, filling, and spin corresponding to a block of bands."""

    def __init__(self,
                 projections: List[Union[str, Dict[str, Any]]],
                 spin: SpinType = None,
                 name: Optional[str] = None,
                 num_wann: Optional[int] = None,
                 num_bands: Optional[int] = None,
                 include_bands: Optional[List[int]] = None,
                 exclude_bands: Optional[str] = None):

        self.projections = []
        for proj in projections:
            if isinstance(proj, str):
                proj = proj_string_to_dict(proj)
            self.projections.append(proj)
        self.id = BlockID(label=name, filled=None, spin=spin)
        self.num_wann = num_wann
        self.num_bands = num_bands
        self.include_bands = include_bands
        self.exclude_bands = exclude_bands
        self.w90_calc: calculators.Wannier90Calculator | None = None

    @property
    def spin(self):
        """Return the spin of the block."""
        return self.id.spin

    @property
    def name(self):
        """Return the name of the block."""
        return self.id.label

    def __repr__(self) -> str:
        out = f'ProjectionBlock({[self.projections]}'
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
        """Convert a ProjectionBlock object to a dictionary."""
        dct = {k: v for k, v in self.__dict__.items()}
        dct.pop('w90_calc')
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @property
    def w90_kwargs(self) -> Dict[str, Any]:
        """Return the `Wannier90` keywords to provide when constructing a new calculator corresponding to this block."""
        kwargs = {}
        for key in ['projections', 'num_wann', 'num_bands', 'exclude_bands']:
            val = getattr(self, key, None)
            if val is None and key != 'exclude_bands':
                raise AttributeError(
                    f'You must define `{self.__class__.__name__}.{key}` before requesting `w90_kwargs`')
            kwargs[key] = val
        if self.spin is not None:
            kwargs['spin'] = self.spin
        return kwargs

    @classmethod
    def fromdict(cls, dct):
        """Construct a ProjectionBlock object from a dictionary."""
        return cls(**dct)


class ProjectionBlocks(object):
    """A collection of blocks of projections.

    In addition to the projections blocks themselves, it also stores system-wide properties such as how many extra
    conduction bands we have.

    Whenever a user queries self.blocks (e.g. when they iterate over this object) it will first propagate these
    system-wide properties down to the individual ProjectionBlock objects. See self.blocks() for more details.
    """

    def __init__(self, blocks: List[ProjectionBlock], atoms: Atoms):
        self._blocks = blocks
        self._atoms = atoms
        # This BandBlocks object must keep track of how many bands we have not belonging to any block
        self.exclude_bands: Dict[Optional[str], List[int]] = {k: [] for k in SpinOptions}
        self.num_extra_bands: Dict[Optional[str], int] = {k: 0 for k in SpinOptions}
        self._num_occ_bands: Dict[Optional[str], int] = {}

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

    def __getitem__(self, key):
        return self.blocks[key]

    def __setitem__(self, key, value):
        self._blocks[key] = value

    def divisions(self, spin: Optional[str]) -> List[int]:
        """Work out the size of individual "blocks" in the set of bands."""
        divs: List[int] = []
        excl_bands = set(self.exclude_bands[spin])
        for block in self.get_subset(spin):
            block_size = num_wann_from_projections(block.projections, self._atoms)
            for excl_band in sorted(excl_bands):
                if sum(divs) == excl_band - 1:
                    excl_bands.remove(excl_band)
                    divs.append(1)
                elif excl_band > sum(divs) and excl_band <= sum(divs) + block_size:
                    raise ValueError('The `Wannier90` excluded bands are mixed with a block of projections. '
                                     'Please redefine excluded_bands such that the remaining included bands are '
                                     'commensurate with the provided projections')
            divs.append(block_size)
        assert len(excl_bands) == 0
        return divs

    @property
    def blocks(self):
        """Return all blocks of projections, with additional global information added to each block."""
        for spin in SpinOptions:
            subset = self.get_subset(spin=spin)
            if len(subset) == 0:
                continue
            is_last = [False for _ in subset]
            is_last[-1] = True
            wann_counter = 1
            for iblock, (b, include_above) in enumerate(zip(subset, is_last)):
                # Construct num_wann
                b.num_wann = num_wann_from_projections(b.projections, self._atoms)

                # Construct num_bands
                b.num_bands = b.num_wann
                if include_above:
                    b.num_bands += self.num_extra_bands[spin]

                # Construct exclude_bands
                while wann_counter in self.exclude_bands[spin]:
                    wann_counter += 1
                band_indices = range(wann_counter, wann_counter + b.num_wann)
                wann_counter += b.num_wann
                if include_above:
                    # For the uppermost block we don't want to exclude any extra bands from the wannierization
                    upper_bound = max(band_indices)
                else:
                    upper_bound = self.num_bands(spin=spin)
                to_exclude = [i for i in range(1, upper_bound + 1) if i not in band_indices]
                b.include_bands = list(band_indices)
                if len(to_exclude) > 0:
                    b.exclude_bands = list_to_formatted_str(to_exclude)

                # Construct the ID
                spin_str = f'_spin_{spin}' if spin is not None else ''
                b.id = BlockID(label=f'block_{iblock + 1}{spin_str}', spin=spin)

        return self._blocks

    def __iter__(self):
        for b in self.blocks:
            yield b

    @classmethod
    def fromlist(cls,
                 list_of_projections: List[List[Union[str, Dict[str, Any]]]],
                 spins: List[SpinType],
                 atoms: Atoms):
        """Construct a ProjectionBlocks object from a list of projections."""
        if not all([isinstance(p, list) for p in list_of_projections]):
            raise ValueError('`list_of_projections` must be a list of lists')
        blocks: List[ProjectionBlock] = []
        blocks += [ProjectionBlock(p, s) for p, s in zip(list_of_projections, spins) if len(p) > 0]

        return cls(blocks, atoms)

    def get_subset(self, spin: Optional[str] = 'both') -> List[ProjectionBlock]:
        """Return all blocks with a particular spin."""
        return [b for b in self._blocks if (spin == 'both' or b.spin == spin)]

    def num_wann(self, spin: Optional[str] = 'both') -> int:
        """Return the number of Wannier functions in the entire system."""
        return sum([num_wann_from_projections(b.projections, self._atoms) for b in self.get_subset(spin)])

    def num_bands(self, spin: Optional[str] = None) -> int:
        """Return the number of bands in the entire system."""
        nbands = self.num_wann(spin)
        nbands += len(self.exclude_bands[spin])
        nbands += self.num_extra_bands[spin]
        return nbands

    @property
    def num_occ_bands(self):
        """The number of occupied bands in the entire system."""
        return self._num_occ_bands

    @num_occ_bands.setter
    def num_occ_bands(self, value: Dict[Optional[str], int]):
        self._num_occ_bands = value

    def todict(self) -> dict:
        """Convert a ProjectionBlocks object to a dictionary."""
        dct: Dict[str, Any] = {k: v for k, v in self.__dict__.items()}
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @classmethod
    def fromdict(cls, dct):
        """Convert a dictionary to a ProjectionBlocks object."""
        new_bandblock = cls(dct.pop('_blocks'), dct.pop('_atoms'))
        for k, v in dct.items():
            if not hasattr(new_bandblock, k):
                raise AttributeError(k)
            setattr(new_bandblock, k, v)
        return new_bandblock

    @property
    def to_merge(self) -> Dict[BlockID, List[ProjectionBlock]]:
        """Determine the sets of blocks that should be merged with one another.

        Group the blocks by their correspondence to occupied/empty bands, and by their spin
        """
        dct: Dict[BlockID, List[ProjectionBlock]] = {}
        for block in self.blocks:
            try:
                n_occ_bands = self.num_occ_bands[block.spin]
            except KeyError:
                raise AssertionError(
                    'Initialize `ProjectionBlocks.num_occ_bands` before calling `ProjectionBlocks.to_merge()`')
            if max(block.include_bands) <= n_occ_bands:
                filled = True
            elif min(block.include_bands) > n_occ_bands:
                filled = False
            else:
                raise ValueError('Block spans both occupied and empty manifolds')
            block_id = BlockID(filled=filled, spin=block.spin)
            if block_id in dct:
                dct[block_id].append(block)
            else:
                dct[block_id] = [block]
        return dct
