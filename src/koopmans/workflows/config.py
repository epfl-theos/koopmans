"""Settings for a `Workflow` object."""

from enum import Enum
from pathlib import Path
from typing import Any

from pydantic import Field, field_validator, model_validator

from koopmans.pydantic import BaseModel


class Task(Enum):
    SINGLEPOINT = 'singlepoint'
    CONVERGENCE = 'convergence'
    WANNIERIZE = 'wannierize'
    UI = 'ui'
    DFT_BANDS = 'dft_bands'
    DFT_EPS = 'dft_eps'
    TRAJECTORY = 'trajectory'


class Correction(Enum):
    KI = 'ki'
    KIPZ = 'kipz'
    PKIPZ = 'pkipz'
    NONE = 'none'
    ALL = 'all'


class CalculateScreeningMethod(Enum):
    DSCF = 'dscf'
    DFPT = 'dfpt'


class VariationalOrbital(Enum):
    PZ = 'pz'
    KOHN_SHAM = 'kohn-sham'
    MLWFS = 'mlwfs'
    PROJWFS = 'projwfs'


class WorkflowConfig(BaseModel):
    """Model for the configuration of a `Workflow`."""

    task: Task = Field(default=Task.SINGLEPOINT, description="Task to perform")
    correction: Correction = Field(
        default=Correction.KI, description="orbital-density-dependent-functional/density-functional to use")
    calculate_alpha: bool = Field(
        default=True, description="whether or not to calculate the screening parameters ab-initio")
    pseudo_library: str = Field(
        description="the pseudopotential library to use (for valid options, run `koopmans pseudos list`)")
    screening_method: CalculateScreeningMethod = Field(
        default=CalculateScreeningMethod.DSCF, description="the method to calculate the screening parameters: either with ΔSCF or DFPT")
    init_orbitals: VariationalOrbital = Field(
        default=VariationalOrbital.PZ, description="which orbitals to use as an initial guess for the variational orbitals")
    init_empty_orbitals: VariationalOrbital = Field(
        default=None, description="which orbitals to use as an initial guess for the empty variational orbitals")
    frozen_orbitals: bool = Field(
        default=False, description="if True, freeze the variational orbitals for the duration of the calculation once they've been initialized")
    calculate_bands: bool = Field(default=False, description="Calculate the band structure of the system (if relevant)")
    spin_polarized: bool = Field(
        default=False, description="if True, the system will be allowed to break spin symmetry i.e. n^{up}(r) != n^{down}(r)")
    initialize_with_smearing: bool = Field(
        default=False, description="if True, the first step of the workflow will use smearing. This can help convergence in some difficult cases.")
    fix_spin_contamination: bool = Field(
        default=False, description="if True, steps will be taken to try and avoid spin contamination. This is only sensible when performing a non-spin-polarized calculation, and is turned on by default for such calculations")
    npool: int | None = Field(
        default=None, description="Number of pools for parallelizing over kpoints (should be commensurate with the k-point grid)")
    gb_correction: bool | None = Field(
        default=None, description="if True, apply the Gygi-Baldereschi scheme to deal with the q->0 divergence of the Coulomb interation for periodic systems")
    mp_correction: bool | None = Field(
        default=None, description="if True, apply the Makov-Payne correction for charged periodic systems")
    mt_correction: bool | None = Field(
        default=None, description="if True, apply the Martyna-Tuckerman correction for charged aperiodic systems")
    eps_inf: float | str | None = Field(
        default=None, description='dielectric constant of the system used by the Gygi-Baldereschi and Makov-Payne corrections; either provide an explicit value or set to "auto" to calculate it ab initio')
    alpha_numsteps: int = Field(default=10, description="Number of steps for alpha calculation")
    alpha_conv_thr: float = Field(
        default=1e-3, description="convergence threshold for |Delta E_i - epsilon_i|; if below this threshold, the corresponding alpha value is not updated")
    alpha_guess: float | list[float] = Field(
        default=0.6, description="starting guess for alpha (overridden if alpha_from_file is true)")
    alpha_mixing: float = Field(default=1.0, description="mixing parameter for updating alpha")
    alpha_from_file: bool = Field(
        default=False, description="if True, uses the file_alpharef.txt from the base directory as a starting guess")
    orbital_groups: list[list[int]] | None = Field(
        default=None, description="a list of integers the same length as the total number of bands, denoting which bands to assign the same screening parameter to")
    orbital_groups_self_hartree_tol: float | None = Field(
        default=None, description="when calculating alpha parameters, the code will group orbitals together only if their self-Hartree energy is within this threshold")
    orbital_groups_spread_tol: float | None = Field(
        default=None, description="when calculating alpha parameters, the code will group orbitals together only if their spread is within this threshold")
    converge: bool = Field(
        default=False, description="If True, repeat the workflow increasing the convergence_parameters until the convergence_observable converges within the convergence_threshold")
    dfpt_coarse_grid: tuple[int, int, int] | None = Field(
        default=None, description="The coarse k-point grid on which to perform the DFPT calculations")
    block_wannierization_threshold: float | None = Field(
        default=None, description="blocks of bands separated by this threshold will be Wannierized separately")
    max_time: int | float | None = Field(
        default=None, description="maximum time in seconds to wait for the workflow to complete; if None, no timeout is applied")
    wait_time: int | float = Field(
        default=5, description="time in seconds to wait between checking the status of in-progress calculations")

    @field_validator('task', 'correction', 'screening_method', 'init_orbitals', 'init_empty_orbitals', mode='before')
    @classmethod
    def make_lowercase(cls, v: Any) -> Any:
        if isinstance(v, str):
            return v.lower()
        return v

    @field_validator('orbital_groups', mode='before')
    @classmethod
    def ensure_orbital_groups_is_list_of_lists(cls, v: Any) -> Any:
        if v is not None:
            if len(v) == 0 or not isinstance(v[0], list):
                v = [v]
        return v

    @model_validator(mode="before")
    def empty_variational_orbitals_default_to_same_as_filled(cls, values: dict[str, Any]) -> dict[str, Any]:
        """If init_empty_orbitals is not specified, set it to the same value as init_orbitals."""
        if values.get("init_empty_orbitals", None) is None:
            values["init_empty_orbitals"] = values.get("init_orbitals", VariationalOrbital.PZ)
        return values
