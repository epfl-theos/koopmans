"""Settings module for the UI calculator."""

from pathlib import Path
from typing import Any, Literal

from ase_koopmans.dft.kpoints import BandPath
from pydantic import Field, field_validator

from koopmans.pydantic import BaseModel


class UnfoldAndInterpolateConfig(BaseModel):
    """Model for the UI calculator settings"""
    kc_ham_file: Path | None = Field(
        default=None,
        description='the name of the Hamiltonian file to read in'
    )
    wannier90_seedname: Path = Field(
        default=Path('wannier90'),
        description=('wannier90_seedname must be equal to the seedname used in the previous Wannier90 calculation. The code '
                     'will look for a file called wannier90_seedname.wout')
    )
    wannier90_calc: Literal['pc', 'sc'] = Field(
        default='pc',
        description=('Specifies the type of PW/Wannier90 calculation preceding the koopmans calculation. If the latter '
                     'is done in a supercell at Gamma then wannier90_calc must be equal to \'sc\', otherwise if it comes from '
                     'a calculation with k-points it must be equal to \'pc\'.\n')
    )
    do_map: bool = Field(
        default=False,
        description=('if True, it realizes the map |m> --> |Rn>, that connects the Wannier functions in the supercell to '
                     'those in the primitive cell. This is basically the unfolding procedure. It can be activated only '
                     'if wannier90_calc=\'sc\'')
    )
    use_ws_distance: bool = Field(
        default=True,
        description=('if True, the real Wigner-Seitz distance between the Wannier functions centers is considered as in '
                     'the Wannier90 code. In particular, this accounts for the periodic boundary conditions and it is '
                     'crucial for a good interpolation when using coarse MP meshes or, equivalently, small supercells')
    )
    kpath: str | BandPath | None = Field(
        default=None,
        description='path in the Brillouin zone for generating the band structure, specified by a string e.g. "GXG"'
    )
    smooth_int_factor: tuple[int, int, int] = Field(
        default=(1, 1, 1),
        description=('if this is > 1 (or is a 3-element list with at least one entry > 1), the smooth interpolation '
                     'method is used. This consists of removing the DFT part of the Hamiltonian from the full Koopmans '
                     'Hamiltonian and adding the DFT Hamiltonian from a calculation with a denser k-points mesh, where '
                     'this keyword defines how many times denser to make the mesh. (If this is set to a scalar a, the '
                     'new k-grid will be [a*kx_old, a*ky_old, a*kz_old]. If it is a list [a, b, c], the dense k-grid '
                     'will be [a*kx_old, b*ky_old, c*kz_old].) This works only for a non self-consistent Koopmans '
                     'calculation using Wannier since, to be consistent, all the Hamiltonians must be in the same '
                     'gauge, i.e. the Wannier gauge')
    )
    dft_ham_file: Path | None = Field(
        default=None,
        description=''
    )
    dft_smooth_ham_file: Path | None = Field(
        default=None,
        description=''
    )
    do_dos: bool = Field(
        default=True,
        description=('if True, the density-of-states is interpolated along the input kpath. The DOS is written to a '
                     'file called "dos_interpolated.dat"')
    )
    num_wann: int | None = Field(
        default=None,
        description=''
    )
    num_wann_sc: int | None = Field(
        default=None,
        description=''
    )
    wannier90_input_sc: bool = Field(
        default=False,
        description='')

    @field_validator('smooth_int_factor', mode='before')
    @classmethod
    def ensure_smooth_int_factor_is_tuple(cls, v: Any) -> Any:
        if isinstance(v, int):
            v = (v, v, v)
        elif isinstance(v, list):
            v = tuple(v)
        return v

    @property
    def do_smooth_interpolation(self):
        """Return True if the smooth interpolation is used."""
        return any([f > 1 for f in self.smooth_int_factor])
