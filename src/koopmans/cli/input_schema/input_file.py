from __future__ import annotations

from collections import OrderedDict
from json import load
from pathlib import Path
from typing import Any, ClassVar, Literal, Type, Union

import numpy as np
from ase_koopmans import Atoms
from ase_koopmans.cell import Cell
from pydantic import Field, field_validator, model_validator
from pydash import get
from upf_tools import UPFDict
from wannier90_input.models.parameters import (Coordinate,
                                               FractionalCoordinate,
                                               Projection)
from yaml import safe_load

from koopmans.calculators.parameters import CalculatorsInputParameters
from koopmans.calculators.parameters.pw import (ControlNamelist,
                                                ElectronsNamelist,
                                                SystemNamelist)
from koopmans.calculators.parameters.pw2wannier90 import \
    PW2Wannier90InputParameters
from koopmans.calculators.parameters.wannier90 import Wannier90InputParameters
from koopmans.commands import CommandConfigs
from koopmans.engines import Engine, LocalhostEngine
from koopmans.kpoints import GammaOnlyKpoints, Kpoints
from koopmans.plotting import PlottingConfig
from koopmans.processes.ui.config import UnfoldAndInterpolateConfig
from koopmans.projections import ExplicitProjections, ImplicitProjections
from koopmans.pydantic import BaseModel
from koopmans.utils import Spin
from koopmans.workflows import (DFTBandsWorkflow, DFTPhWorkflow,
                                SinglepointWorkflow, TrajectoryWorkflow,
                                WannierizeWorkflow, Workflow, WorkflowInput)
from koopmans.workflows.config import Task, WorkflowConfig

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

    def to_kpoints_obj(self, cell: Cell) -> Kpoints | GammaOnlyKpoints:
        if self.gamma_only:
            return GammaOnlyKpoints()
        elif self.path is None:
            path = cell.get_bravais_lattice().bandpath()
            return Kpoints(path=path, grid=self.grid, offset=self.offset)
        else:
            return Kpoints.from_path_as_string(path=self.path, cell=cell, density=self.density,
                                               grid=self.grid, offset=self.offset)


class PWInput(BaseModel):
    control: ControlNamelist = Field(default_factory=lambda: ControlNamelist())
    system: SystemNamelist = Field(default_factory=lambda: SystemNamelist())
    electrons: ElectronsNamelist = Field(default_factory=lambda: ElectronsNamelist())


class SpinSpecificWannierInput(BaseModel):
    projections: list[list[Projection]] = []
    dis_froz_max: float | None = None
    dis_froz_min: float | None = None
    dis_win_max: float | None = None
    dis_win_min: float | None = None
    auto_projections: bool = False


class Wannier90Input(Wannier90InputParameters):
    # Hide various fields that should not be set in the w90 block of the Koopmans input file
    num_wann: ClassVar[int] = -1
    num_bands: ClassVar[int | None] = -1
    unit_cell_cart: ClassVar[list[Coordinate]] = []
    mp_grid: ClassVar[tuple[int, int, int]] = (1, 1, 1)
    kpoints: ClassVar[list[FractionalCoordinate]] = []
    # In the Koopmans input file, projections are further split into blocks
    projections: list[list[Projection]] = []
    up: SpinSpecificWannierInput = Field(default_factory=lambda: SpinSpecificWannierInput())
    down: SpinSpecificWannierInput = Field(default_factory=lambda: SpinSpecificWannierInput())

    @model_validator(mode="before")
    @classmethod
    def set_default_num_bands(cls, values: dict[str, Any]) -> dict[str, Any]:
        """Skip default num_bands setting."""
        return values

    @model_validator(mode="before")
    @classmethod
    def atoms_frac_xor_cart(cls, values: dict[str, Any]) -> dict[str, Any]:
        """Skip atoms_frac/atoms_cart check."""
        return values


class CalculatorParametersInput(BaseModel):
    ecutwfc: float | None = None
    nbnd: float | None = None
    pw: PWInput = Field(default_factory=lambda: PWInput())
    pw2wannier90: PW2Wannier90InputParameters = Field(default_factory=lambda: PW2Wannier90InputParameters())
    wannier90: Wannier90Input = Field(default_factory=lambda: Wannier90Input())
    ui: UnfoldAndInterpolateConfig = Field(default_factory=lambda: UnfoldAndInterpolateConfig())

    def to_obj(self) -> CalculatorsInputParameters:
        dct = self.model_dump()
        ecutwfc = dct.pop('ecutwfc')
        nbnd = dct.pop('nbnd')
        if ecutwfc is not None:
            self.pw.system.ecutwfc = ecutwfc
        if nbnd is not None:
            self.pw.system.nbnd = nbnd
        return CalculatorsInputParameters(dct)


class EngineInput(BaseModel):
    from_scratch: bool = True
    keep_tmpdirs: bool = True
    commands: CommandConfigs = Field(default_factory=lambda: CommandConfigs)


class KoopmansInput(BaseModel):
    workflow: WorkflowConfig
    engine: EngineInput = Field(default_factory=lambda: EngineInput)
    atoms: AtomsInput
    kpoints: KpointsInput
    calculator_parameters: CalculatorParametersInput
    plotting: PlottingConfig

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

    def get_pseudopotentials(self, atoms: Atoms, engine: Engine) -> OrderedDict[str, UPFDict]:
        """Get the pseudopotentials as UPFDict objects."""
        upfs = OrderedDict()

        # Work out the unique atom labels. This needs to be done like this rather than with set() to ensure
        # the order is consistent with how the atoms appear in the input file
        unique_labels: list[tuple[str, str]] = []
        for a in atoms:
            label_tuple: tuple[str, str] = (a.symbol, a.tag)
            if label_tuple not in unique_labels:
                unique_labels.append(label_tuple)

        for element, tag in unique_labels:
            pseudo = engine.get_pseudopotential(library=self.workflow.pseudo_library, element=element)

            pseudo_type = get(pseudo, 'header.pseudo_type', None)
            if pseudo_type != 'NC':
                raise ValueError('Koopmans functionals only currently support norm-conserving pseudopotentials; '
                                 f'{pseudo.filename} is {pseudo_type}')

            symbol = element + str(tag) if tag > 0 else element
            upfs[symbol] = pseudo

        return upfs

    def get_projections(self, atoms: Atoms, pseudopotentials: OrderedDict[str, UPFDict]) -> ExplicitProjections | ImplicitProjections:
        if self.workflow.spin_polarized:
            projs_by_spin = [self.calculator_parameters.wannier90.up.projections,
                             self.calculator_parameters.wannier90.down.projections]
            spins = [Spin.UP, Spin.DOWN]
        else:
            projs_by_spin = [self.calculator_parameters.wannier90.projections]
            spins = [Spin.NONE]

        auto = any([self.calculator_parameters.wannier90.auto_projections,
                    self.calculator_parameters.wannier90.up.auto_projections,
                    self.calculator_parameters.wannier90.down.auto_projections])
        if auto:
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
        else:
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
            engine = LocalhostEngine.model_validate(self.engine.model_dump())

        # Here need to define the logic constructing the workflow properly
        atoms = self.atoms.to_atoms_obj()
        kpoints = self.kpoints.to_kpoints_obj(cell=atoms.cell)
        pseudopotentials = self.get_pseudopotentials(atoms, engine)
        projections = self.get_projections(atoms, pseudopotentials)
        calculator_parameters = self.calculator_parameters.to_obj()

        inputs = WorkflowInput(
            atoms=atoms,
            kpoints=kpoints,
            pseudopotentials=pseudopotentials,
            projections=projections
        )

        return workflow_cls(
            inputs=inputs,
            config=self.workflow,
            atoms=atoms,
            kpoints=kpoints,
        )
