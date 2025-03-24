from pathlib import Path
from typing import Any, Dict, List, Optional, TypeVar, Union

from ase_koopmans import Atoms
from ase_koopmans.io.wannier90 import (list_to_formatted_str,
                                       num_wann_from_projections,
                                       proj_string_to_dict)
from pydantic import BaseModel, ConfigDict, model_validator

from koopmans import calculators
from koopmans.utils import Spin


class BlockID(BaseModel):

    label: Optional[str] = None
    filled: Optional[bool] = None
    spin: Spin = Spin.NONE
    model_config = ConfigDict(frozen=True)

    def __str__(self):
        if self.label is None:
            return ''
        return f'{self.label}_{self.spin}' if self.spin else self.label

    # Set label to 'occ' or 'emp' if filled is set but label is not
    @model_validator(mode='before')
    def set_label_from_filled(cls, data: Any):
        if isinstance(data, dict):
            if data.get('filled', None) is not None and data.get('label', None) is None:
                spin = data.get('spin', None)
                spin_label = '_spin_' + spin if spin is not None else ''
                data['label'] = 'occ' + spin_label if data['filled'] else 'emp' + spin_label
        return data

    def filling(self):
        if self.filled:
            return 'occ'
        elif self.filled is False:
            return 'emp'
        else:
            raise ValueError(f'The filling of block {str(self)} is not known')


TProjBlock = TypeVar("TProjBlock", bound="ProjectionBlock")


class ProjectionBlock(object):
    # List of keys to be provided by self.w90_kwargs
    _w90_keys = ['num_wann', 'num_bands', 'exclude_bands']

    # This simple object contains the projections, filling, and spin corresponding to a block of bands
    def __init__(self,
                 num_wann: int,
                 num_bands: Optional[int] = None,
                 spin: Spin = Spin.NONE,
                 name: Optional[str] = None,
                 include_bands: Optional[List[int]] = None,
                 exclude_bands: Optional[str] = None):

        self.id = BlockID(label=name, filled=None, spin=spin)
        self.num_wann = num_wann
        self.num_bands = num_bands
        self.include_bands = include_bands
        self.exclude_bands = exclude_bands
        self.w90_calc: calculators.Wannier90Calculator | None = None
        raise NotImplementedError('spin is read-only')
        # self.spin = spin

    @property
    def spin(self):
        return self.id.spin

    @property
    def name(self):
        return self.id.label

    def __repr__(self) -> str:
        out = f'{self.__class__.__name__}(num_wann={[self.num_wann]}'
        if self.spin is not None:
            out += f', spin={self.spin}'
        return out + ')'

    def __len__(self):
        # The number of Projection functions
        return self.num_wann

    def __bool__(self):
        # Return true if this ProjectionBlock is non-empty
        return len(self) > 0

    def __eq__(self, other):
        if isinstance(other, ProjectionBlock):
            return self.__dict__ == other.__dict__
        return False

    def todict(self) -> dict:
        dct = {k: v for k, v in self.__dict__.items()}
        dct.pop('w90_calc')
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @property
    def w90_kwargs(self) -> Dict[str, Any]:
        # Returns the keywords to provide when constructing a new calculator corresponding to this block
        kwargs = {}
        for key in self._w90_keys:
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
        return cls(**dct)

    def split(self: TProjBlock, groups: List[List[int]], atoms: Optional[Atoms] = None) -> List[TProjBlock]:
        # Split this block into several blocks, as indicated by "groups", a list of lists of orbital indices, with each
        # sublist corresponding to a new block to be created. For the base class "atoms" is not required, but it is
        # required for subclasses

        assert self.num_wann is not None
        assert self.num_bands is not None
        assert self.include_bands is not None

        mismatched_bands = set([i for subblock in groups for i in subblock]).difference(self.include_bands)
        if len(mismatched_bands) > 0:
            raise ValueError('The list of oribtals "groups" does not perfectly overlap with this block')

        # Create the two blocks
        dct = {k: v for k, v in self.__dict__.items() if k in ['num_wann', 'projections']}
        blocks = [self.__class__(**dct) for _ in groups]
        for block, group in zip(blocks, groups):
            # num_wann
            block.num_wann = len(group)
            # num_bands (equal to num_wann, since after splitting we don't have extra bands)
            block.num_bands = len(group)

        return blocks


class ExplicitProjectionBlock(ProjectionBlock):
    # This class extends ProjectionBlock to have explicit projections instead of simply num_wann
    def __init__(self, projections: List[str | Dict[str, Any]], *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.projections = []
        for proj in projections:
            if isinstance(proj, str):
                proj = proj_string_to_dict(proj)
            self.projections.append(proj)

        self._w90_keys.append('projections')

    def __repr__(self) -> str:
        out = f'{self.__class__.__name__}({[self.projections]}'
        if self.spin is not None:
            out += f', spin={self.spin}'
        return out + ')'

    def __len__(self):
        # Count the number of projections
        return len(self.projections)

    def split(self, groups: List[List[int]], atoms: Optional[Atoms] = None) -> List[ExplicitProjectionBlock]:
        blocks = super().split(groups)

        if atoms is None:
            raise ValueError('Please provide an Atoms object when splitting a block with explicit projections')

        # Split self.projections
        for block in blocks:
            block.projections = []
        i_block = 0
        for i, proj in enumerate(self.projections):
            current_num_wann = num_wann_from_projections(blocks[i_block].projections, atoms)
            target_num_wann = len(groups[i_block])
            additional_num_wann = num_wann_from_projections([proj], atoms)

            if current_num_wann == target_num_wann:
                i_block += 1
            if current_num_wann + additional_num_wann <= target_num_wann:
                blocks[i_block].projections.append(proj)
            else:
                raise ValueError(
                    f'The projections \n{self.projections}\n cannot be split into blocks of length ' + ', '.join([str(len(x)) for x in groups]))

        return blocks


class ProjectionBlocks(object):
    """
    This object is a collection of blocks of bands that will be separately wannierized. In addition to the blocks
    of bands themselves, it also stores system-wide properties such as how many extra conduction bands we have.

    Whenever a user queries self.blocks (e.g. when they iterate over this object) it will first propagate these
    system-wide properties down to the individual ProjectionBlock objects. See self.blocks for more details.
    """

    def __init__(self, blocks: List[ProjectionBlock]):
        self._blocks = blocks
        # This BandBlocks object must keep track of how many bands we have not belonging to any block
        self.exclude_bands: Dict[Spin, List[int]] = {k: [] for k in Spin}
        self.num_extra_bands: Dict[Spin, int] = {k: 0 for k in Spin}
        self._num_occ_bands: Dict[Spin, int] = {}

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

    def divisions(self, spin: Spin) -> List[int]:
        # This algorithm works out the size of individual "blocks" in the set of bands
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

    @property
    def blocks(self):
        # Before returning all the blocks, add more global information such as exclude_bands
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
                spin_str = f'_spin_{spin}' if spin is not None else ''
                b.id = BlockID(label=f'block_{iblock + 1}{spin_str}', spin=spin)

        return self._blocks

    def __iter__(self):
        for b in self.blocks:
            yield b

    @classmethod
    def fromlist(cls,
                 list_of_projections: List[List[Union[str, Dict[str, Any]]]],
                 spins: List[Spin],
                 atoms: Atoms):

        if not all([isinstance(p, list) for p in list_of_projections]):
            raise ValueError('`list_of_projections` must be a list of lists')
        blocks: List[ProjectionBlock] = []
        for projs, spin in zip(list_of_projections, spins):
            if isinstance(projs, int):
                blocks.append(ProjectionBlock(projs, spin=spin))
            else:
                assert isinstance(projs, list)
                if len(projs) == 0:
                    continue
                if atoms is None:
                    raise ValueError('To construct a ProjectionBlocks object from a list of list of '
                                     'strings/dictionaries, you must provide an "atoms" argument')
                for i, proj in enumerate(projs):
                    if isinstance(proj, str):
                        projs[i] = proj_string_to_dict(proj)
                num_wann = num_wann_from_projections(projs, atoms)
                blocks.append(ExplicitProjectionBlock(projs, num_wann, spin=spin))

        return cls(blocks)

    def get_subset(self, spin: Spin | str = 'both') -> List[ProjectionBlock]:
        return [b for b in self._blocks if (spin == 'both' or b.spin == spin)]

    def num_wann(self, spin: Spin | str = 'both') -> int:
        return sum([b.num_wann for b in self.get_subset(spin)])

    def num_bands(self, spin: Spin = Spin.NONE) -> int:
        nbands = self.num_wann(spin)
        nbands += len(self.exclude_bands[spin])
        nbands += self.num_extra_bands[spin]
        return nbands

    @property
    def num_occ_bands(self):
        return self._num_occ_bands

    @num_occ_bands.setter
    def num_occ_bands(self, value: Dict[Spin, int]):
        self._num_occ_bands = value

    def todict(self) -> dict:
        dct: Dict[str, Any] = {k: v for k, v in self.__dict__.items()}
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @classmethod
    def fromdict(cls, dct):
        new_bandblock = cls(dct.pop('_blocks'))
        for k, v in dct.items():
            if not hasattr(new_bandblock, k):
                raise AttributeError(k)
            setattr(new_bandblock, k, v)
        return new_bandblock

    def block_spans_occ_and_emp(self, block: ProjectionBlock) -> bool:
        # Works out if the provided block spans both occupied and empty
        try:
            n_occ_bands = self.num_occ_bands[block.spin]
        except KeyError:
            raise AssertionError(
                'Initialize ProjectionBlocks.num_occ_bands before calling ProjectionBlocks.to_merge()')
        assert block.include_bands is not None
        return max(block.include_bands) > n_occ_bands and min(block.include_bands) <= n_occ_bands

    def block_is_occupied(self, block: ProjectionBlock) -> bool:
        # Works out if the provided block is occupied
        try:
            n_occ_bands = self.num_occ_bands[block.spin]
        except KeyError:
            raise AssertionError(
                'Initialize ProjectionBlocks.num_occ_bands before calling ProjectionBlocks.to_merge()')
        assert block.include_bands is not None
        if max(block.include_bands) <= n_occ_bands:
            return True
        elif min(block.include_bands) > n_occ_bands:
            return False
        else:
            raise ValueError('Block spans both occupied and empty manifolds')

    @property
    def to_merge(self) -> Dict[BlockID, List[ProjectionBlock]]:
        # Group the blocks by their correspondence to occupied/empty bands, and by their spin
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

    def split(self, block: ProjectionBlock, groups: Optional[List[List[int]]] = None, atoms: Optional[Atoms] = None):
        # Split block
        assert block.include_bands is not None

        if groups is None:
            # Default to splitting valence and conduction
            groups = [[i for i in block.include_bands if i <= self.num_occ_bands[block.spin]],
                      [i for i in block.include_bands if i > self.num_occ_bands[block.spin]]]

        blocks = block.split(groups, atoms)

        # Record the current self.num_wann
        num_wann_orig = self.num_wann(block.spin)

        # Replace block with blocks in self._blocks
        i = self._blocks.index(block)
        self._blocks = self._blocks[:i] + blocks + self._blocks[i+1:]

        # Check self.num_wann has not changed
        assert self.num_wann(block.spin) == num_wann_orig

        # After splitting we do away with extra bands
        self.num_extra_bands[block.spin] = 0

        # Call self.blocks to re-populate the global information assigned to each block
        self.blocks

        return

    def to_wjl_config(self, fd):
        dct = {}
        dct['splitvc'] = {'indices': [p.include_bands for p in self],
                          'outdirs': [str(p.directory) for p in self]}
        toml.dump(dct, fd)
