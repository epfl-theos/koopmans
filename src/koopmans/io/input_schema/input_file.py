from __future__ import annotations

from collections import OrderedDict
from json import load
from pathlib import Path
from typing import Any, ClassVar, Literal, Type, Union

import numpy as np
from ase import Atoms
from ase.cell import Cell
from pydantic import Field, field_validator, model_validator
from pydantic_core import ErrorDetails
from pydash import get
from upf_tools.pydantic import UPFModel
from wannier90_input.models.parameters import Projection
from yaml import safe_load

from koopmans.commands import CommandConfigs
from koopmans.engines.ase import LocalASEEngine
from koopmans.engines.engine import Engine
from koopmans.kpoints import Kpoints, kpoints_from_path_as_string
from koopmans.parameters import CalculatorsInputParameters
from koopmans.parameters.kcp import KCPInputParameters
from koopmans.parameters.pw import PWInputParameters
from koopmans.parameters.pw2wannier90 import PW2Wannier90InputParameters
from koopmans.parameters.ui import UnfoldAndInterpolateConfig
from koopmans.parameters.wannier90 import RestrictedWannier90InputParameters
from koopmans.plotting import PlottingConfig
from koopmans.projections import ExplicitProjections, ImplicitProjections
from koopmans.pseudopotentials import nwfcs_from_projectors, nwfcs_from_pseudos
from koopmans.pydantic import BaseModel
from koopmans.utils import Spin
from koopmans.workflows import (DFTBandsWorkflow, DFTPhWorkflow,
                                SinglepointWorkflow, TrajectoryWorkflow,
                                WannierizeWorkflow, Workflow, WorkflowInput)
from koopmans.workflows.config import Task, WorkflowConfig
from koopmans.workflows.utils import get_pseudopotentials

from .atomic_positions import AtomicAlatPositionsInput, AtomicPositionsInput
from .cell_parameters import (CellParametersViaAlat, CellParametersViaIbrav,
                              CellParametersViaVectors)


class AtomsInput(BaseModel):
    cell_parameters: CellParametersViaIbrav | CellParametersViaVectors | CellParametersViaAlat
    atomic_positions: AtomicPositionsInput | AtomicAlatPositionsInput

    def to_atoms_obj(self) -> Atoms:
        cell = self.cell_parameters.to_cell()
        symbols = self.atomic_positions.symbols

        position_kwarg: dict[str, Any] = {}
        if self.atomic_positions.units == 'crystal':
            position_kwarg['scaled_positions'] = self.atomic_positions.positions_array()
        else:
            position_kwarg['positions'] = self.atomic_positions.in_angstrom()

        atoms = Atoms(symbols=symbols, cell=cell, pbc=self.cell_parameters.periodic, **position_kwarg)

        atoms.set_tags(self.atomic_positions.tags)

        return atoms

    def to_snapshots(self) -> list[Atoms]:
        """Generate a list of Atoms objects for each snapshot, if applicable."""
        if self.atomic_positions.snapshots is None:
            return [self.to_atoms_obj()]

        raise NotImplementedError('Need to implement this')


class KpointsInput(BaseModel):
    grid: tuple[int, int, int]
    offset: tuple[int, int, int] = (0, 0, 0)
    path: str | None = None
    gamma_only: bool = False
    density: float = 10.0

    def to_kpoints_obj(self, cell: Cell) -> Kpoints:
        if self.gamma_only:
            return Kpoints(gamma_only=True)
        elif self.path is None:
            # By default, use ASE's default bandpath for this cell (see
            # https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#brillouin-zone-data)
            default_path = cell.bandpath(eps=1e-10).path
            return Kpoints(path=default_path, grid=self.grid, offset=self.offset)
        else:
            return kpoints_from_path_as_string(path=self.path, cell=cell, density=self.density,
                                               grid=self.grid, offset=self.offset)


class SpinSpecificWannierInput(BaseModel):
    dis_froz_max: float | None = None
    dis_froz_min: float | None = None
    dis_win_max: float | None = None
    dis_win_min: float | None = None
    auto_projections: bool = False
    projections: list[list[Projection]] = Field(default_factory=list)

    @field_validator('projections', mode='before')
    @classmethod
    def convert_strings_to_projections(cls, value: Any) -> list[list[Projection]]:
        return [[Projection.from_string(proj) if isinstance(proj, str) else proj for proj in block] for block in value]


class Wannier90InputParametersWithUpDown(RestrictedWannier90InputParameters):
    up: SpinSpecificWannierInput | None = None
    down: SpinSpecificWannierInput | None = None
    # In the input file, projections are specified as a list of lists to separate each block
    projections: list[list[Projection]] = Field(default_factory=list)

    @model_validator(mode='after')
    def check_up_down_exclusivity(self) -> Wannier90InputParametersWithUpDown:
        if (self.up is None) != (self.down is None):
            raise ValueError("Both 'up' and 'down' must be specified together.")
        return self


class CalculatorParametersInput(BaseModel):
    ecutwfc: float | None = None
    nbnd: float | None = None
    tot_magnetization: float | None = None
    pw: PWInputParameters = Field(default_factory=lambda: PWInputParameters())
    pw2wannier90: PW2Wannier90InputParameters = Field(default_factory=lambda: PW2Wannier90InputParameters())
    wannier90: Wannier90InputParametersWithUpDown = Field(default_factory=lambda: Wannier90InputParametersWithUpDown())
    ui: UnfoldAndInterpolateConfig = Field(default_factory=lambda: UnfoldAndInterpolateConfig())
    kcp: KCPInputParameters = Field(default_factory=lambda: KCPInputParameters())

    def to_obj(self, spin_polarized: bool = False) -> dict[str, Any]:
        dct = self.model_dump()
        ecutwfc = dct.pop('ecutwfc')
        nbnd = dct.pop('nbnd')
        tot_magnetization = dct.pop('tot_magnetization')
        if ecutwfc is not None:
            dct['pw']['system']['ecutwfc'] = ecutwfc
            dct['kcp']['system']['ecutwfc'] = ecutwfc
        if nbnd is not None:
            dct['pw']['system']['nbnd'] = nbnd
        if tot_magnetization is not None:
            dct['pw']['system']['tot_magnetization'] = tot_magnetization
            dct['kcp']['system']['tot_magnetization'] = tot_magnetization

        # Split wannier90 up and down into their own blocks
        up_dct = dct['wannier90'].pop('up')
        down_dct = dct['wannier90'].pop('down')
        if spin_polarized:
            wannier90 = dct.pop('wannier90')
            dct['wannier90_up'] = wannier90.copy()
            dct['wannier90_down'] = wannier90.copy()
            if up_dct:
                dct['wannier90_up'].update(up_dct)
            if down_dct:
                dct['wannier90_down'].update(down_dct)

        for key in ['wannier90_up', 'wannier90_down', 'wannier90']:
            # Projections are excluded from the CalculatorsInputParameters object
            if key in dct:
                dct[key].pop('projections')

        return dct


class KoopmansInput(BaseModel):
    workflow: WorkflowConfig
    atoms: AtomsInput
    kpoints: KpointsInput = Field(default_factory=lambda: KpointsInput(gamma_only=True))
    calculator_parameters: CalculatorParametersInput
    plotting: PlottingConfig = Field(default_factory=lambda: PlottingConfig())

    @classmethod
    def from_file(cls, filename: str | Path) -> KoopmansInput:
        """Load an input file and return a KoopmansInput object."""
        filename = Path(filename)
        if filename.suffix in {'.yaml', '.yml'}:
            with open(filename, 'r') as f:
                input_dict = safe_load(f)
        elif filename.suffix == '.json':
            with open(filename, 'r') as f:
                input_dict = load(f)
        else:
            raise ValueError(f'Unrecognized file type for `{filename}`')

        return cls.model_validate(input_dict)

    def get_projections(self, atoms: Atoms, pseudopotentials: OrderedDict[str, UPFModel]) -> ExplicitProjections | ImplicitProjections:
        spins = [Spin.UP, Spin.DOWN] if self.workflow.spin_polarized else [Spin.NONE]
        if self.workflow.automated_wannierization:
            return self.get_implicit_projections(atoms, pseudopotentials, spins)
        else:
            return self.get_explicit_projections(atoms, pseudopotentials, spins)

    def get_implicit_projections(self, atoms: Atoms, pseudopotentials: OrderedDict[str, UPFModel], spins: list[Spin]) -> ImplicitProjections:
        if self.calculator_parameters.pw2wannier90.atom_proj_ext:
            proj_dir = self.calculator_parameters.pw2wannier90.atom_proj_dir
            if proj_dir is None:
                raise ValueError("When `atom_proj_ext` is True, `atom_proj_dir` must be specified.")
            projectors = OrderedDict([(k, proj_dir / p.filename.with_suffix('.dat').name)
                                      for k, p in pseudopotentials.items()])
            num_wann = [nwfcs_from_projectors(atoms, projectors) for _ in spins]
        else:
            num_wann = [nwfcs_from_pseudos(atoms, pseudopotentials) for _ in spins]
        return ImplicitProjections.from_block_lengths(num_wann, spins, atoms)

    def get_explicit_projections(self, atoms: Atoms, pseudopotentials: OrderedDict[str, UPFModel], spins: list[Spin]) -> ExplicitProjections:
        if self.workflow.spin_polarized:
            projs_by_spin = [self.calculator_parameters.wannier90.up.projections,
                             self.calculator_parameters.wannier90.down.projections]
        else:
            projs_by_spin = [self.calculator_parameters.wannier90.projections]
        blocks = []
        for spin, projs in zip(spins, projs_by_spin):
            blocks += [{'spin': spin, 'projections': proj_block} for proj_block in projs]
        return ExplicitProjections(blocks=blocks, atoms=atoms)

    def to_workflow(self, engine: Engine | None = None) -> Workflow[Any]:
        """Generate a Workflow object from the input data."""

        workflow_cls: Type[Workflow[Any]]
        if self.workflow.task == Task.SINGLEPOINT:
            workflow_cls = SinglepointWorkflow
        elif self.workflow.task == Task.WANNIERIZE:
            workflow_cls = WannierizeWorkflow
        elif self.workflow.task == Task.DFT_BANDS:
            workflow_cls = DFTBandsWorkflow
        elif self.workflow.task == Task.DFT_EPS:
            workflow_cls = DFTPhWorkflow
        elif self.workflow.task == Task.TRAJECTORY:
            workflow_cls = TrajectoryWorkflow
        else:
            raise NotImplementedError(f'Need to implement `{self.workflow.task.value}`')

        if engine is None:
            engine = LocalASEEngine()

        # Here need to define the logic constructing the workflow properly
        atoms = self.atoms.to_atoms_obj()
        kpoints = self.kpoints.to_kpoints_obj(cell=atoms.cell)
        pseudopotentials = get_pseudopotentials(atoms, engine, self.workflow.pseudo_library)
        projections = self.get_projections(atoms, pseudopotentials)
        calculator_parameters = self.calculator_parameters.to_obj(spin_polarized=self.workflow.spin_polarized)

        inputs = WorkflowInput(
            config=self.workflow,
            atoms=atoms,
            kpoints=kpoints,
            pseudopotentials=pseudopotentials,
            projections=projections,
            calculator_parameters=calculator_parameters,
            plotting_config=self.plotting,
        )

        return workflow_cls(
            name=None,
            inputs=inputs,
            engine=engine,
        )


CUSTOM_MESSAGES = {
    'type': 'is the wrong type (should be "{expected_type}", not "{given_type}")',
    'extra_forbidden': 'is not a valid keyword.',
    'missing': 'was not provided.',
}


def convert_errors(e: ValidationError) -> list[ErrorDetails]:
    """Make the validation errors more user-friendly"""
    new_errors: list[ErrorDetails] = []
    for error in e.errors():
        custom_message = CUSTOM_MESSAGES.get(error['type'], None)
        if custom_message:
            ctx = error.get('ctx')
            error['msg'] = custom_message.format(**ctx) if ctx else custom_message
        new_errors.append(error)
    return new_errors


def prettify_errors(e: ValidationError) -> str:
    """Return a prettified string of validation errors."""
    errors = convert_errors(e)
    error_lines = []
    for error in errors:
        loc = '.'.join(str(l) for l in error['loc'])
        msg = error['msg']
        error_lines.append(f' `{loc}` {msg}')
    return '\n'.join(error_lines)


def read_input_file(filename: str | Path, engine: Engine | None = None) -> Workflow[Any]:
    try:
        koopmans_input = KoopmansInput.from_file(filename)
    except ValidationError as e:
        raise ValueError(
            f"Errors found in the input file: \n\n{prettify_errors(e)}\n\n For more information, see URL_HERE.")
    return koopmans_input.to_workflow(engine=engine)
