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
from wannier90_input.models.parameters import Projection
from yaml import safe_load

from koopmans.commands import CommandConfigs
from koopmans.engines import Engine, LocalhostEngine
from koopmans.kpoints import GammaOnlyKpoints, Kpoints
from koopmans.parameters import CalculatorsInputParameters
from koopmans.parameters.kcp import KCPInputParameters
from koopmans.parameters.pw import PWInputParameters
from koopmans.parameters.pw2wannier90 import PW2Wannier90InputParameters
from koopmans.parameters.ui import UnfoldAndInterpolateConfig
from koopmans.parameters.wannier90 import Wannier90InputParameters
from koopmans.plotting import PlottingConfig
from koopmans.projections import ExplicitProjections, ImplicitProjections
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

    def to_kpoints_obj(self, cell: Cell) -> Kpoints | GammaOnlyKpoints:
        if self.gamma_only:
            return GammaOnlyKpoints()
        elif self.path is None:
            # By default, use ASE's default bandpath for this cell (see
            # https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#brillouin-zone-data)
            default_path = atoms.cell.bandpath(eps=1e-10).path
            return Kpoints(path=default_path, grid=self.grid, offset=self.offset)
        else:
            return Kpoints.from_path_as_string(path=self.path, cell=cell, density=self.density,
                                               grid=self.grid, offset=self.offset)


class SpinSpecificWannierInput(BaseModel):
    dis_froz_max: float | None = None
    dis_froz_min: float | None = None
    dis_win_max: float | None = None
    dis_win_min: float | None = None
    auto_projections: bool = False
    projections: list[list[Projection]] = Field(default_factory=list)


class Wannier90InputParametersWithUpDown(Wannier90InputParameters):
    up: SpinSpecificWannierInput = Field(default_factory=lambda: SpinSpecificWannierInput())
    down: SpinSpecificWannierInput = Field(default_factory=lambda: SpinSpecificWannierInput())
    # In the input file, projections are specified as a list of lists to separate each block
    projections: list[list[Projection]] = Field(default_factory=list)


class CalculatorParametersInput(BaseModel):
    ecutwfc: float | None = None
    nbnd: float | None = None
    pw: PWInputParameters = Field(default_factory=lambda: PWInputParameters())
    pw2wannier90: PW2Wannier90InputParameters = Field(default_factory=lambda: PW2Wannier90InputParameters())
    wannier90: Wannier90InputParametersWithUpDown = Field(default_factory=lambda: Wannier90InputParametersWithUpDown())
    ui: UnfoldAndInterpolateConfig = Field(default_factory=lambda: UnfoldAndInterpolateConfig())
    kcp: KCPInputParameters = Field(default_factory=lambda: KCPInputParameters())

    def to_obj(self) -> dict[str, Any]:
        dct = self.model_dump()
        ecutwfc = dct.pop('ecutwfc')
        nbnd = dct.pop('nbnd')
        if ecutwfc is not None:
            self.pw.system.ecutwfc = ecutwfc
        if nbnd is not None:
            self.pw.system.nbnd = nbnd

        # Split wannier90 up and down into their own blocks
        up_dct = dct['wannier90'].pop('up')
        down_dct = dct['wannier90'].pop('down')
        dct['wannier90_up'] = dct['wannier90'].copy()
        dct['wannier90_down'] = dct['wannier90'].copy()
        dct['wannier90_up'].update(up_dct)
        dct['wannier90_down'].update(down_dct)

        for key in ['wannier90_up', 'wannier90_down', 'wannier90']:
            # Projections are excluded from the CalculatorsInputParameters object
            dct[key].pop('projections')

        return dct


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
        pseudopotentials = get_pseudopotentials(atoms, engine, self.workflow.pseudo_library)
        projections = self.get_projections(atoms, pseudopotentials)
        calculator_parameters = self.calculator_parameters.to_obj()

        inputs = WorkflowInput(
            config=self.workflow,
            atoms=atoms,
            kpoints=kpoints,
            pseudopotentials=pseudopotentials,
            projections=projections,
            calculator_parameters=calculator_parameters,
        )

        return workflow_cls(
            name=None,
            inputs=inputs,
            engine=engine,
        )
