"""Classes that define blocks of bands."""
from __future__ import annotations

import string
from abc import ABC
from typing import Any, Dict, Generic, List, Optional, TypeVar, Union

from ase_koopmans import Atoms
from ase_koopmans.io.wannier90 import (list_to_formatted_str,
                                       num_wann_from_projections)
from pydantic import BaseModel, ConfigDict, Field, model_validator
from wannier90_input.models.parameters import Projection

from koopmans import calculators
from koopmans.utils import Spin


class BlockID(BaseModel):
    """The ID of a block of bands, which includes a label, filling, and spin."""

    label: Optional[str] = None
    filled: Optional[bool] = None
    spin: Spin = Spin.NONE
    model_config = ConfigDict(frozen=True)

    def __str__(self):
        if self.label is None:
            return ''
        return f'{self.label}_{self.spin}' if self.spin != Spin.NONE else self.label

    # Set label to 'occ' or 'emp' if filled is set but label is not
    @model_validator(mode='before')
    def set_label_from_filled(cls, data: Any):
        """Construct a label (if it is not provided) from the filling and spin arguments."""
        if isinstance(data, dict):
            if data.get('filled', None) is not None and data.get('label', None) is None:
                spin = data.get('spin', Spin.NONE)
                spin_label = f'_spin_{spin}' if spin != Spin.NONE else ''
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


class ProjectionsBlock(BaseModel):
    """A class that contains the projections, filling, and spin corresponding to a block of bands."""

    num_wann: int
    filled: Optional[bool] = None
    spin: Spin = Spin.NONE
    label: Optional[str] = None
    num_bands: Optional[int] = None
    include_bands: Optional[List[int]] = None
    exclude_bands: Optional[List[int]] = None
    w90_calc: calculators.Wannier90Calculator | None = None

    model_config = ConfigDict(frozen=False, arbitrary_types_allowed=True)

    @property
    def id(self):
        """Return the ID of this block."""
        return BlockID(label=self.label, filled=self.filled, spin=self.spin)

    def __len__(self) -> int:
        return self.num_wann

    def __bool__(self):
        # Return true if this ProjectionBlock is non-empty
        return len(self) > 0

    def __eq__(self, other):
        if isinstance(other, ProjectionsBlock):
            return self.__dict__ == other.__dict__
        return False

    @property
    def w90_kwargs(self) -> Dict[str, Any]:
        """Return the `Wannier90` keywords to provide when constructing a new calculator corresponding to this block."""
        kwargs = {}
        for key in ['num_wann', 'num_bands', 'exclude_bands']:
            val = getattr(self, key, None)
            if val is None and key != 'exclude_bands':
                raise AttributeError(
                    f'You must define `{self.__class__.__name__}.{key}` before requesting `w90_kwargs`')
            kwargs[key] = val
        if self.spin is not Spin.NONE:
            kwargs['spin'] = self.spin.value
        return kwargs


class ImplicitProjectionsBlock(ProjectionsBlock):
    """This class implements ProjectionsBlock with automated projections."""

    def split(self, groups: List[List[int]]) -> List[ImplicitProjectionsBlock]:
        """Split the block into sub-blocks according to their groupings.

        Note that these blocks will have `num_wann` == `num_bands` and will not have any disentanglement.
        It is assumed that any disentanglement has already been performed on the full block.
        """
        # Sanity checking
        band_indices = [i for group in groups for i in group]
        if not len(band_indices) == len(set(band_indices)):
            raise ValueError('The provided groups contain duplicate band indices')

        assert self.include_bands is not None
        if not set(band_indices) == set(self.include_bands):
            raise ValueError('The provided groups do not span the same bands as this block of bands')

        # Construct the sub-groups
        blocks = []
        for letter, include_bands in zip(string.ascii_lowercase, groups):
            exclude_bands = [i for i in self.include_bands if i not in include_bands]
            new_block = ImplicitProjectionsBlock(num_wann=len(include_bands),
                                                 num_bands=len(include_bands),
                                                 spin=self.spin,
                                                 filled=self.filled,
                                                 label=self.label + letter if self.label else None,
                                                 include_bands=include_bands,
                                                 exclude_bands=exclude_bands)
            blocks.append(new_block)

        return blocks

    @property
    def w90_kwargs(self) -> Dict[str, Any]:
        """Return the `Wannier90` keywords to provide when constructing a new calculator corresponding to this block."""
        kwargs = super().w90_kwargs
        kwargs['auto_projections'] = True
        return kwargs


class ExplicitProjectionsBlock(ProjectionsBlock):
    """This class extends ProjectionsBlock to have explicit projections instead of simply num_wann."""

    projections: list[Projection]

    @property
    def w90_kwargs(self) -> Dict[str, Any]:
        """Return the `Wannier90` keywords to provide when constructing a new calculator corresponding to this block."""
        kwargs = super().w90_kwargs
        kwargs['auto_projections'] = False
        kwargs['projections'] = [p.dict() for p in self.projections]
        return kwargs


ProjectionsBlockType = TypeVar('ProjectionsBlockType', bound=ProjectionsBlock)


class Projections(BaseModel, ABC, Generic[ProjectionsBlockType]):
    """A collection of blocks of projections.

    In addition to the projections blocks themselves, it also
    stores system-wide properties such as how many extra conduction bands we have.

    Whenever a user iterates through this object it will first propagate these system-wide properties down
    to the individual ProjectionBlock objects. See self._populate_blocks() for more details.
    """

    blocks: list[ProjectionsBlockType]
    atoms: Atoms
    exclude_bands: Dict[Spin, List[int]] = Field(default_factory=lambda: {k: [] for k in Spin})
    num_extra_bands: Dict[Spin, int | None] = Field(default_factory=lambda: {k: 0 for k in Spin})
    num_occ_bands: Dict[Spin, int | None] = Field(default_factory=lambda: {k: None for k in Spin})

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def __repr__(self) -> str:
        out = self.__class__.__name__ + '('
        for b in self.blocks:
            out += f'{b}\n                  '
        out = out.strip()
        out += ')'
        return out

    def __len__(self):
        # Count the number of non-empty "ProjectionBlock"s
        self._populate_blocks()
        return sum([1 for b in self.blocks if b])

    def __bool__(self):
        # Return true if this contains any non-empty "ProjectionBlock"s
        return len(self) > 0

    def __eq__(self, other):
        if isinstance(other, Projections):
            return self.__dict__ == other.__dict__
        return False

    def __getitem__(self, key):
        self._populate_blocks()
        return self.blocks[key]

    def __setitem__(self, key, value):
        self.blocks[key] = value

    def divisions(self, spin: Spin) -> List[int]:
        """Work out the size of individual "blocks" in the set of bands."""
        divs: List[int] = []
        excl_bands = set(self.exclude_bands[spin])
        for block in self.get_subset(spin):
            for excl_band in sorted(excl_bands):
                if sum(divs) == excl_band - 1:
                    excl_bands.remove(excl_band)
                    divs.append(1)
                elif excl_band > sum(divs) and excl_band <= sum(divs) + block.num_wann:
                    raise ValueError('The `Wannier90` excluded bands are mixed with a block of projections. '
                                     'Please redefine excluded_bands such that the remaining included bands are '
                                     'commensurate with the provided projections')
            divs.append(block.num_wann)
        assert len(excl_bands) == 0
        return divs

    def _populate_blocks(self):
        """Populate all blocks of projections with additional global information."""
        for spin in Spin:
            subset = self.get_subset(spin=spin)
            if len(subset) == 0:
                continue
            is_last = [False for _ in subset]
            is_last[-1] = True
            wann_counter = 1
            for iblock, (b, include_above) in enumerate(zip(subset, is_last)):
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
                spin_str = '' if spin == Spin.NONE else f'_spin_{spin.value}'
                b.label = f'block_{iblock + 1}{spin_str}'

    def __iter__(self):
        # Before returning all the blocks, add more global information such as exclude_bands
        self._populate_blocks()

        for b in self.blocks:
            yield b

    def get_subset(self, spin: Spin | str = 'both') -> List[ProjectionsBlockType]:
        """Return all blocks with a particular spin."""
        return [b for b in self.blocks if (spin == 'both' or b.spin == spin)]

    def num_wann(self, spin: Spin | str = 'both') -> int:
        """Return the number of Wannier functions in the entire system."""
        return sum([b.num_wann for b in self.get_subset(spin)])

    def num_bands(self, spin: Spin) -> int:
        """Return the number of bands in the entire system."""
        nbands = self.num_wann(spin)
        nbands += len(self.exclude_bands[spin])
        num_extra_bands = self.num_extra_bands[spin]
        if num_extra_bands is None:
            raise ValueError(f"`num_extra_bands` is not set for spin = {spin}")
        nbands += num_extra_bands
        return nbands

    @property
    def to_merge(self) -> Dict[BlockID, List[ProjectionsBlockType]]:
        """Determine the sets of blocks that should be merged with one another.

        Group the blocks by their correspondence to occupied/empty bands, and by their spin
        """
        dct: Dict[BlockID, List[ProjectionsBlockType]] = {}
        for block in self:
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


class ExplicitProjections(Projections[ExplicitProjectionsBlock]):
    """A set of projections with explicitly specified projections."""

    @classmethod
    def fromlist(cls,
                 list_of_projections: List[List[Union[str, Dict[str, Any]]]],
                 spins: List[Spin],
                 atoms: Atoms):
        """Create a set of projections from a list of projections."""
        if not all([isinstance(p, list) for p in list_of_projections]):
            raise ValueError('`list_of_projections` must be a list of lists')
        blocks = [ExplicitProjectionsBlock(projections=p, spin=s, num_wann=num_wann_from_projections(
            p, atoms)) for p, s in zip(list_of_projections, spins) if len(p) > 0]
        # TODO achieve the above via a model validator to populate the projections
        return cls(blocks=blocks, atoms=atoms)


class ImplicitProjections(Projections[ImplicitProjectionsBlock]):
    """A set of projections with the projections specified only via num_wann."""

    @classmethod
    def from_block_lengths(cls, lengths: List[int], spins: List[Spin], atoms: Atoms):
        """Construct a set of implicit projections purely from the number of Wannier orbitals."""
        blocks = [ImplicitProjectionsBlock(num_wann=nw, spin=s) for nw, s in zip(lengths, spins) if nw > 0]
        return cls(blocks=blocks, atoms=atoms)
