"""Generic workflow object for koopmans."""

from __future__ import annotations

import copy
import json as json_ext
import logging
import os
import re
import sys
import time
from abc import ABC, abstractmethod
from collections import OrderedDict
from contextlib import contextmanager
from pathlib import Path
from typing import (Any, Callable, Dict, Generator, Generic, List, Optional,
                    Sequence, Tuple, Type, TypeVar, Union)

import dill
import numpy as np
from numpy import typing as npt
from pybtex.database import BibliographyData
from pydantic import Field, field_validator, model_validator
from pydash import get, set_
from typing_extensions import Self

# isort: off
import koopmans.mpl_config  # noqa: F401
import matplotlib.pyplot as plt
# isort: on

from ase import Atoms
from ase.build.supercells import make_supercell
from ase.cell import Cell
from ase.dft.dos import DOS
from ase.dft.kpoints import BandPath
from ase.spacegroup import symmetrize
from ase.spectrum.band_structure import BandStructure
from ase.spectrum.doscollection import GridDOSCollection
from ase.spectrum.dosdata import GridDOSData
from ase_koopmans.io import read as ase_read
from ase_koopmans.io import write as ase_write
from ase_koopmans.io.espresso import \
    contruct_kcp_namelist as construct_namelist
from upf_tools.pydantic import UPFModel

from koopmans import settings, utils
from koopmans.calculators.calculator import CalcType, CalculatorProcess
from koopmans.calculators.pw import PWCalculator
from koopmans.calculators.wannier90 import Wannier90Calculator
from koopmans.engines.ase import LocalASEEngine
from koopmans.engines.engine import Engine
from koopmans.files import File
from koopmans.kpoints import Kpoints
from koopmans.ml import AbstractMLModel, MLModel, OccEmpMLModels
from koopmans.parameters import CalculatorsInputParameters
from koopmans.parameters.kcp import KCPInputParameters
from koopmans.parameters.ml import MLConfig
from koopmans.parameters.pw import PWInputParameters
from koopmans.plotting import PlottingConfig
from koopmans.processes.koopmans_cp import (ConvertFilesFromSpin1To2,
                                            ConvertFilesFromSpin2To1)
from koopmans.processes.process import Process
from koopmans.projections import (ExplicitProjections, ImplicitProjections,
                                  Projections)
from koopmans.protocol import ProcessProtocol
from koopmans.pseudopotentials import (nelec_from_pseudos,
                                       nwfcs_from_projectors,
                                       nwfcs_from_pseudos,
                                       pseudopotential_library_citations)
from koopmans.pydantic import BaseModel
from koopmans.references import bib_data
from koopmans.status import Status
from koopmans.utils import Spin
from koopmans.utils.warnings import ImageCorrectionWarning, warn
from koopmans.variational_orbitals import VariationalOrbitals
from koopmans.workflows.config import (CalculateScreeningMethod, Correction,
                                       Task, VariationalOrbital,
                                       WorkflowConfig)
from koopmans.workflows.utils import get_pseudopotentials


class WorkflowInput(BaseModel):
    config: WorkflowConfig
    atoms: Atoms = Field(description='an ASE ``Atoms`` object defining the atomic positions, cell, etc')
    pseudopotentials: OrderedDict[str, UPFModel] = Field(
        description='a dictionary mapping atom labels to pseudopotential objects')
    kpoints: Kpoints
    projections: ExplicitProjections | ImplicitProjections | None = None
    variational_orbitals: VariationalOrbitals | None = None
    calculator_parameters: CalculatorsInputParameters = Field(default_factory=CalculatorsInputParameters)
    plotting_config: PlottingConfig = Field(default_factory=PlottingConfig)
    ml_config: MLConfig = Field(default_factory=lambda: MLConfig())
    ml_model: AbstractMLModel | None = None
    snapshots: List[Atoms]

    @model_validator(mode='before')
    @classmethod
    def setup_implicit_projections(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Set up the projections if they are implicit."""
        calculator_parameters = get(data, 'calculator_parameters')
        auto = [get(calculator_parameters, f'wannier90{key}.auto_projections')
                for key in ['', '_up', '_down']]
        if any(auto):
            raise NotImplementedError()
            # Also override various pw2wannier and wannier90 settings
            calculator_parameters['pw2wannier'].atom_proj = True
            for spin in Spin:
                label = 'wannier90' if spin == Spin.NONE else f'wannier90_{spin.value}'
                calculator_parameters[label].guiding_centres = False
                calculator_parameters[label].dis_froz_proj = True
                if calculator_parameters[label].dis_proj_max is None:
                    calculator_parameters[label].dis_proj_max = 0.95
            data['calculator_parameters'] = calculator_parameters
        return data

    @model_validator(mode='after')
    def ml_checks(self) -> Self:
        """Perform checks on the ML configuration."""
        if self.ml_config.predict or self.ml_config.train or self.ml_config.test:
            warn("Predicting screening parameters with machine-learning is an experimental feature; proceed with "
                 "caution")
            if self.ml_model is None:
                raise ValueError("You have requested to train or predict with a machine-learning model, but no model "
                                 "is attached to this workflow. Either set `ml:train` or `predict` to `True` when "
                                 "initializing the workflow, or directly add a model to the workflow's `ml_model` "
                                 "attribute")
            if self.ml_model.estimator_type != self.ml_config.estimator:
                raise ValueError(f'The estimator type of the loaded ML model (`{self.ml_model.estimator_type}`) does not '
                                 f'match the estimator type specified in the Workflow settings (`{self.ml_config.estimator}`).')
            if self.ml_model.descriptor_type != self.ml_config.descriptor:
                raise ValueError(f'The descriptor type of the loaded ML model (`{self.ml_model.descriptor_type}`) does not '
                                 f'match the descriptor type specified in the Workflow settings (`{self.ml_config.descriptor}`).')
            if self.config.task not in [Task.SINGLEPOINT, Task.TRAJECTORY]:
                raise NotImplementedError(
                    f'Using the ML-prediction for the task `{self.config.task.value}` has not yet been implemented.')
            if self.config.screening_method != CalculateScreeningMethod.DSCF:
                raise NotImplementedError(
                    f"Using the ML-prediction for the method `{self.config.screening_method.value}` has not yet been implemented")
            if self.config.correction != Correction.KI:
                raise NotImplementedError(
                    f'Using the ML-prediction for the `{self.config.correction.value}` functional has not yet been '
                    'implemented.')
            if self.config.init_orbitals != VariationalOrbital.MLWFS:
                raise NotImplementedError(
                    f'Using the ML-prediction for `{self.config.init_orbitals.value}` initial orbitals has not yet been '
                    'implemented.')
            if self.config.init_empty_orbitals != self.config.init_orbitals.value:
                raise NotImplementedError(
                    'Using the ML-prediction for using different initial orbitals for empty states than for occupied '
                    'states has not yet been implemented.')
            if self.config.spin_polarized:
                warn('Using the ML-prediction for spin-polarised systems has not yet been extensively tested.')
            if not all(self.atoms.pbc):
                warn('Using the ML-prediction for non-periodic systems has not yet been extensively tested.')
            if self.config.orbital_groups:
                warn('Using orbital_groups has not yet been extensively tested.')
            if not np.all(self.atoms.cell.angles() == 90.0):
                raise ValueError("The ML-workflow has only been implemented for simulation cells that have 90° angles")

        return self

    @model_validator(mode='before')
    @classmethod
    def populate_snapshots(cls, values: dict[str, Any]) -> dict[str, Any]:
        if get(values, 'snapshots', None) is None:
            set_(values, 'snapshots', [get(values, 'atoms')])
        return values

    @field_validator('atoms', mode="before")
    @classmethod
    def wrap_atoms(cls, value: Atoms) -> Atoms:
        if all(value.pbc):
            value.wrap(pbc=True)
        return value

    @model_validator(mode='before')
    @classmethod
    def add_excluded_bands_info_to_projections(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Add excluded_bands info to the Projections"""
        projections = get(data, 'projections', None)
        calculator_parameters = get(data, 'calculator_parameters')
        if projections:
            exclude_bands = get(projections, 'exclude_bands')
            for spin in [Spin.UP, Spin.DOWN, Spin.NONE]:
                label = 'wannier90'
                if spin != Spin.NONE:
                    label += f'_{spin.value}'
                exclude_bands[spin] = get(calculator_parameters, f'{label}.exclude_bands', [])
        return data

    @field_validator('atoms', mode='before')
    @classmethod
    def detach_atoms_calculator(cls, value: Atoms) -> Atoms:
        if value.calc is not None:
            warn(f'You have initialized a `{cls.__name__}` object with an atoms object that '
                 'possesses a calculator. This calculator will be ignored.')
            value.calc = None
        return value

    @model_validator(mode='after')
    def ensure_wannier90_config_and_spin_polarization_are_consistent(self) -> Self:
        if self.config.spin_polarized != (self.calculator_parameters.wannier90 is None):
            import ipdb
            ipdb.set_trace()
            raise ValueError("Inconsistent spin polarization settings: for spin-polarized calculations, "
                             "both 'wannier90_up' and 'wannier90_down' must be specified; for non-spin-polarized calculations, "
                             "'wannier90' must be specified.")
        return self

    @model_validator(mode='before')
    @classmethod
    def extract_projections(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Adding the projections to the workflow kwargs.

        This is unusual in that this is an attribute of the workflow but can be specified within the wannier90 block
        """
        calc_configs = get(data, 'calculator_parameters')
        atoms = get(data, 'atoms')
        projections = get(data, 'projections', None)
        spin_polarized = get(data, 'config.spin_polarized')

        if projections is None:
            raise NotImplementedError('Need to check the following code')
            wannier90_block_projs: list[dict[str, Any]] = []

            spins = [Spin.UP, Spin.DOWN] if spin_polarized else [Spin.NONE]
            suffixes = ['up', 'down'] if spin_polarized else ['']
            for spin, suffix in zip(spins, suffixes):
                projs = get(calc_configs, f'wannier90{suffix}.projections', [[]])

                # Add the spin information
                wannier90_block_projs += [{'spin': spin, 'projections': proj_block} for proj_block in projs]

            # Suppress type-checking because mypy does not understand that pydantic will cast blocks appropriately
            # to the desired type
            set_(data, 'projections', ExplicitProjections(blocks=wannier90_block_projs, atoms=atoms))
        return data

    @model_validator(mode='before')
    @classmethod
    def default_fix_spin_contamination(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Set `fix_spin_contamination` to the opposite of `spin_polarized`, if not explicitly set."""
        config = get(data, 'config')
        if get(config, 'fix_spin_contamination') is None:
            set_(config, 'fix_spin_contamination', not get(config, 'spin_polarized'))
        return data

    @model_validator(mode='after')
    def check_if_fix_spin_contamination_then_not_spin_polarized(self) -> Self:
        """Check that if `fix_spin_contamination` is True, then `spin_polarized` is False."""
        if self.config.fix_spin_contamination and self.config.spin_polarized:
            raise ValueError('`fix_spin_contamination = True` is incompatible with `spin_polarized = True`')
        return self

    @model_validator(mode='before')
    @classmethod
    def set_default_frozen_orbitals(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Set `frozen_orbitals` to True if `method` is `dfpt`, False otherwise, if not explicitly set."""
        config = get(data, 'config')
        if get(config, 'frozen_orbitals') is None:
            if get(config, 'screening_method') == CalculateScreeningMethod.DFPT:
                set_(config, 'frozen_orbitals', True)
            else:
                set_(config, 'frozen_orbitals', False)
        return data

    @model_validator(mode='after')
    def check_frozen_orbitals_and_method_consistency(self) -> Self:
        """Check that `frozen_orbitals` is True if `method` is `dfpt`."""
        if self.config.screening_method == CalculateScreeningMethod.DFPT and not self.config.frozen_orbitals:
            raise ValueError('`frozen_orbitals` must be equal to `True` when `method` is `dfpt`')

        if self.config.screening_method == CalculateScreeningMethod.DSCF and self.config.frozen_orbitals:
            warn('You have requested a ΔSCF calculation with frozen orbitals. This is unusual; proceed '
                 'only if you know what you are doing')
        return self

    @model_validator(mode='before')
    @classmethod
    def do_not_calculate_alpha_for_dft(cls, data: dict[str, Any]) -> dict[str, Any]:
        config = get(data, 'config')
        if get(config, 'calculate_alpha') is None and get(config, 'correction') == Correction.NONE:
            set_(config, 'calculate_alpha', False)
        return data

    @model_validator(mode='before')
    @classmethod
    def default_periodic_correction_schemes(cls, data: dict[str, Any]) -> dict[str, Any]:
        atoms = get(data, 'atoms')
        config = get(data, 'config')

        if not get(config, 'calculate_alpha'):
            # If we are not calculating alpha, we do not consider charged systems and therefore we don't need image
            # corrections, so we skip the following checks
            return data
        if all(atoms.pbc):

            if get(config, 'screening_method') == CalculateScreeningMethod.DFPT:
                # For DPFT, we use gb_correction
                if get(config, 'gb_correction') is None:
                    set_(config, 'gb_correction', True)

            elif get(config, 'screening_method') == CalculateScreeningMethod.DSCF:
                # For DSCF, we use mp_correction
                if get(config, 'mp_correction') is None:
                    set_(config, 'mp_correction', True)

                if get(config, 'eps_inf') is None:
                    warn('`eps_inf` missing in input; it will default to 1.0. Proceed with caution for periodic '
                         'systems; consider setting `eps_inf == "auto"` to calculate it automatically.')
                    set_(config, 'eps_inf', 1.0)

            if get(config, 'mt_correction') is None:
                set_(config, 'mt_correction', False)

        else:
            if get(config, 'gb_correction') is None:
                set_(config, 'gb_correction', False)

            if get(config, 'mp_correction') is None:
                set_(config, 'mp_correction', False)

            if get(config, 'mt_correction') is None:
                set_(config, 'mt_correction', True)
        return data

    @model_validator(mode='after')
    def check_periodic_correction_schemes(self) -> Self:
        """Perform sanity checks on the periodic image correction schemes."""
        if not self.config.calculate_alpha:
            # If we are not calculating alpha, we do not consider charged systems and therefore we don't need image
            # corrections, so we skip the following checks
            return self

        if all(self.atoms.pbc):
            if self.config.screening_method == CalculateScreeningMethod.DFPT:
                if not self.config.gb_correction:
                    warn('Gygi-Baldereschi corrections are not being used; do this with '
                         'caution for periodic systems', ImageCorrectionWarning)

            elif self.config.screening_method == CalculateScreeningMethod.DSCF:
                if not self.config.mp_correction:
                    warn('Makov-Payne corrections are not being used; do this with '
                         'caution for periodic systems', ImageCorrectionWarning)

            if self.config.mt_correction:
                raise ValueError('Do not use Martyna-Tuckerman corrections for periodic systems')

            # Check symmetry of the system
            dataset = symmetrize.check_symmetry(self.atoms, 1e-6, verbose=False)
            if dataset is None or dataset.number not in range(195, 231):
                warn('This system is not cubic and will therefore not have a uniform dielectric tensor. However, '
                     'the image-correction schemes that are currently implemented assume a uniform dielectric. '
                     'Proceed with caution', ImageCorrectionWarning)
        else:
            if self.config.gb_correction:
                raise ValueError('Do not use Gygi-Baldereschi corrections for aperiodic systems')

            if self.config.mp_correction:
                raise ValueError('Do not use Makov-Payne corrections for aperiodic systems')

            if not self.config.mt_correction:
                warn('Martyna-Tuckerman corrections not applied for an aperiodic calculation; do this with '
                     'caution', ImageCorrectionWarning)

        return self

    @model_validator(mode='after')
    def check_initial_orbitals(self) -> Self:
        if self.config.init_orbitals in [VariationalOrbital.MLWFS, VariationalOrbital.PROJWFS] \
                and self.config.correction != Correction.NONE:
            if self.projections is None or len(self.projections) == 0:
                raise ValueError(f'In order to use `init_orbitals={self.config.init_orbitals}`, projections must '
                                 'be provided')
        return self

    @model_validator(mode='after')
    def check_projections_spin_consistency(self) -> Self:
        if self.projections is None:
            return self

        spin_set = set([p.spin for p in self.projections])
        if self.config.spin_polarized:
            if spin_set != {Spin.UP, Spin.DOWN}:
                raise ValueError('This calculation is spin-polarized; please provide spin-up and spin-down '
                                 'projections')
        else:
            if spin_set != {Spin.NONE}:
                raise ValueError('This calculation is not spin-polarized; please do not provide spin-indexed '
                                 'projections')
        return self

    @model_validator(mode='before')
    @classmethod
    def check_kponts_vs_do_wf_cmplx(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Check the consistency between self.kpoints.gamma_only and KCP's do_wf_cmplx."""
        if not get(data, 'kpoints.gamma_only') and not get(data, 'calculator_parameters.kcp.system.do_wf_cmplx'):
            warn('`do_wf_cmplx = False` is not consistent with `gamma_only = False`. '
                 'Changing `do_wf_cmplx` to `True`')
            set_(data, 'calculator_parameters.kcp.system.do_wf_cmplx', True)
        return data

    @model_validator(mode='before')
    @classmethod
    def set_smearing_if_initializing_with_smearing(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Set smearing parameters if initializing with smearing."""
        if get(data, 'config.initialize_with_smearing'):
            if get(data, 'calculator_parameters.pw.occupations') != 'smearing':
                set_(data, 'calculator_parameters.pw.occupations', 'smearing')
            if get(data, 'calculator_parameters.pw.smearing') is None:
                set_(data, 'calculator_parameters.pw.smearing', 'cold')
            if get(data, 'calculator_parameters.pw.degauss') is None:
                set_(data, 'calculator_parameters.pw.degauss', 2e-3)
        return data


OutputModel = TypeVar('OutputModel', bound=BaseModel)
InputModel = TypeVar('InputModel', bound=WorkflowInput)
W = TypeVar('W', bound='Workflow')


class Workflow(Process[InputModel, OutputModel], ABC, Generic[InputModel, OutputModel]):
    """Abstract base class that defines a Koopmans workflow."""

    version: str = Field(
        'Version number used to create this workflow (becomes relevant when loading workflows from .kwf files)')
    step_counter: int = 0
    print_indent: int = 0
    status: Status = Status.NOT_STARTED
    _variational_orbitals: VariationalOrbitals | None

    @model_validator(mode='before')
    @classmethod
    def create_ml_model(cls, data: dict[str, Any]) -> dict[str, Any]:
        ml_config = data['inputs'].ml_config
        ml_model = data['inputs'].ml_model
        engine = data['engine']
        if ml_model is None:
            if ml_config.train:
                if ml_config.occ_and_emp_together:
                    ml_model = MLModel(ml_config.estimator, ml_config.descriptor, engine=engine)
                else:
                    ml_model = OccEmpMLModels(ml_config.estimator, ml_config.descriptor, engine=engine)
            elif ml_config.predict or ml_config.test:
                if ml_config.model_file is None:
                    raise ValueError('Cannot initialize a Workflow with `ml.predict = True` without providing a model via '
                                     'the `ml.model_file` setting or the `ml_model` argument')
                with open(ml_config.model_file, 'rb') as f:
                    ml_model = dill.load(f)
                assert isinstance(ml_model, AbstractMLModel)
            data['inputs'].ml_model = ml_model
        return data

    @field_validator('name', mode='before')
    @classmethod
    def set_name(cls, name: Optional[str]) -> str:
        if name:
            return name
        name_with_split_acroynms = re.sub(r'([a-z])([A-Z])', r'\1 \2',
                                          cls.__name__.replace('Workflow', ''))
        return re.sub(r'([A-Z])([A-Z][a-z])', r'\1 \2', name_with_split_acroynms)

    @field_validator('version', mode='before')
    @classmethod
    def set_version(cls, version: Optional[str]) -> str:
        from koopmans import __version__
        return version if version is not None else __version__

    @field_validator('version', mode='after')
    @classmethod
    def check_version(cls, version: str) -> str:
        from koopmans import __version__
        if version != __version__:
            warn(f'You are using version {__version__} of `koopmans`, but this workflow was generated with '
                 f'version {version}. Proceed with caution.')
        return version

    @model_validator(mode="before")
    @classmethod
    def populate_calculator_parameters(cls, data: dict[str, Any]) -> dict[str, Any]:

        if get(data, 'config.task') == Task.UI or get(data, 'autogenerate_settings', False):
            return data

        # Automatically calculate nelec/nelup/neldw/etc using information contained in the pseudopotential files
        # and the kcp settings
        calculator_parameters = get(data, 'inputs.calculator_parameters')
        kcp_params = get(calculator_parameters, 'kcp', {})
        atoms = get(data, 'inputs.atoms')
        pseudopotentials = get(data, 'inputs.pseudopotentials')
        spin_polarized = get(data, 'config.spin_polarized', False)
        nelec = _number_of_electrons(atoms, pseudopotentials, spin_polarized, kcp_params, Spin.NONE)
        nelup = _number_of_electrons(atoms, pseudopotentials, spin_polarized, kcp_params, Spin.UP)
        neldw = _number_of_electrons(atoms, pseudopotentials, spin_polarized, kcp_params, Spin.DOWN)
        tot_charge = nelec_from_pseudos(atoms, pseudopotentials) - nelec
        tot_mag = nelup - neldw

        # Setting up the magnetic moments
        starting_mags = calculator_parameters.pw.system.starting_magnetization
        if not starting_mags:
            starting_mags_list = calculator_parameters.kcp.system.starting_magnetization

        if starting_mags:
            raise NotImplementedError('Yet to test the following')
            labels = [a.symbol + str(a.tag) if a.tag > 0 else a.symbol for a in atoms]
            starting_magmoms: Dict[str, float] = {}
            for i, (l, p) in enumerate(pseudopotentials.items()):
                mag = starting_mags.pop(f'starting_magnetization({i + 1})', 0.0)
                if abs(mag) < 1.0:
                    # If |mag| < 1, QE interprets this as site magnetization *per valence electron*, whereas ASE
                    # expects simply the site magnetization
                    valence = p['header']['z_valence']
                    starting_magmoms[l] = mag * valence
                else:
                    # If |mag| >= 1, QE interprets this as site magnetization
                    starting_magmoms[l] = mag
            atoms.set_initial_magnetic_moments([starting_magmoms[label] for label in labels])
        elif tot_mag != 0:
            atoms.set_initial_magnetic_moments([tot_mag / len(atoms) for _ in atoms])

        # Work out the number of bands
        nbnd = get(calculator_parameters, 'pw.system.nbnd', get(
            calculator_parameters, 'kcp.system.nbnd', nelec // 2 + nelec % 2))
        generated_keywords = {'nelec': nelec, 'tot_charge': tot_charge, 'tot_magnetization': tot_mag,
                              'nelup': nelup, 'neldw': neldw, 'nbnd': nbnd}

        for block in [calculator_parameters.pw, calculator_parameters.kcp]:
            # Apply auto-generated keywords
            for k, v in generated_keywords.items():
                # Skipping nbnd for kcp; if necessary it will be generated later
                if isinstance(block, KCPInputParameters) and k == 'nbnd':
                    continue

                if k in type(block.system).model_fields:
                    setattr(block.system, k, v)

        return data

    @model_validator(mode='before')
    @classmethod
    def get_pseudopotentials(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Get the pseudopotentials, ensuring that they are coerced into UPFModel objects.

        Note this is done inside `Workflow` and not `WorkflowInput` because it requires access
        to the engine and config."""
        engine = data['engine']
        atoms = get(data, 'inputs.atoms')
        pseudopotentials = get(data, 'inputs.pseudopotentials', None)
        pseudo_library = get(data, 'config.pseudo_library')

        if pseudopotentials is None:
            upfs = get_pseudopotentials(atoms, engine, pseudo_library)
        else:
            upfs = OrderedDict()
            for symbol, pseudo_filename in pseudopotentials.items():
                if isinstance(pseudo_filename, str):
                    pseudopotentials['symbol'] = engine.get_pseudopotential(pseudo_library, filename=pseudo_filename)
        return data

    def model_post_init(self, context: Any) -> None:
        self.engine.workflows[self.uid] = self
        return super().model_post_init(context)

    def __eq__(self, other: Any):
        if self.__class__ != other.__class__:
            return False
        else:
            for key in self.__slots__:
                if key == 'engine':
                    continue
                if getattr(self, key) != getattr(other, key):
                    return False
        return True

    def run(self):
        """Run the workflow."""
        assert self.engine is not None

        self.print_preamble()

        if not self.parent_process_uid:
            bf = '**' if sys.stdout.isatty() else ''
            self.print(bf + self.name + bf)
            self.print(bf + '-' * len(self.name) + bf)
            if self.engine.from_scratch:
                self.engine.remove_tmpdirs()

        start_time = time.time()

        while self.status != Status.COMPLETED:
            logger = logging.getLogger(__name__)
            logger.info(f'Proceeding with workflow {self.name} ({self.__class__.__name__})')

            self.proceed()

            while self.engine.processes_are_running():
                logger.info('Processes are running; waiting for them to complete...')
                self.wait()
                if self.parameters.max_time is not None and time.time() - start_time > self.parameters.max_time:
                    self.status = Status.FAILED
                    raise ValueError('Workflow failed to complete within the specified time')

        if self.parent_process_uid is None and self.status == Status.COMPLETED:
            self._teardown()

        self.status = Status.COMPLETED

    def wait(self):
        """Wait while calculations are running.

        While waiting, the engine might perform other actions such as fetching remote files etc.
        """
        self.engine.wait(seconds=self.parameters.wait_time)

    def _pre_run(self):
        """Run checks and actions before running the workflow."""
        logger = logging.getLogger(__name__)
        logger.info(f'Proceeding to the next step in the workflow {self.name} ({self.__class__.__name__})')
        self.step_counter = 0

        if self.status == Status.NOT_STARTED:
            self.status = Status.RUNNING

    def proceed(self):
        """Proceed through the workflow until a step has status `RUNNING`.

        Typically, this will involve proceeding to the next step in the workflow which has not been run,
        set it `RUNNING`, and then returning

        Note, however, for some engines this function will execute the _entire_ workflow e.g. `localhost` workflows,
        which, when requested to run a step in the workflow, will run that step and block the interpreter until the
        step is completed. For such a workflow, the status of individual processes proceed immediately from `NOT_STARTED`
        to `COMPLETE` and `proceed()` will not exit until all the processes in the workflow are `COMPLETED`.
        """

        self._pre_run()

        self._run()

        self._post_run()

    @abstractmethod
    def _run(self) -> None:
        """Run the workflow; to be implemented by subclasses."""

    def _post_run(self):
        """Run checks and actions after running the workflow."""
        logger = logging.getLogger(__name__)
        logger.info(f'Exiting workflow {self.name} with status {self.status.name} ({self.__class__.__name__})')
        assert self.engine is not None
        self.engine.update_statuses()

    @property
    def outputs(self) -> OutputModel:
        """The outputs of the workflow."""
        if self._outputs is None:
            raise ValueError('Process has no outputs because it has not been run yet')
        return self._outputs

    @outputs.setter
    def outputs(self, value: OutputModel):
        self._outputs = value

    @classmethod
    def from_other(cls: Type[W], other_wf: Workflow[Any, Any], **kwargs) -> W:
        """Create a new workflow from another workflow, copying all settings and parameters.

        e.g.
        >>> new_wf = Workflow.fromother(other_wf)
        """
        inputs_dct = {k: getattr(other_wf.inputs, k)
                      for k in other_wf.inputs.model_fields if k in cls.input_model.model_fields}
        inputs_dct.update(kwargs)
        inputs = cls.input_model.model_validate(inputs_dct)
        return cls(inputs=inputs, engine=other_wf.engine)

    @classmethod
    def fromparent(cls: Type[W], parent_wf: Workflow[Any, Any], **kwargs) -> W:
        """Create a subworkflow with the same configuration as the parent workflow.

        e.g.
        >>> sub_wf = Workflow.fromparent(self)
        """
        parent_wf.step_counter += 1

        wf = cls.from_other(
            other_wf=parent_wf,
            **kwargs
        )
        wf.print_indent = parent_wf.print_indent + 2
        wf.parent_process_uid = parent_wf.uid

        subdirectory_str = wf.name.replace(' ', '-').lower()
        subdirectory_path = Path(f'{parent_wf.step_counter:02}-{subdirectory_str}')
        wf.directory = parent_wf.directory / subdirectory_path if parent_wf.directory else subdirectory_path

        return wf

    def new_calculator(self,
                       calc_class: type[CalcType],
                       name: str,
                       spin: Spin = Spin.NONE,
                       override: dict[str, Any] = {},
                       kpts: Optional[Union[List[int], BandPath]] = None,
                       **kwargs) -> CalcType:
        """Create a new calculator of a given type.

        Additional calculator keywords can be provided via `override`
        and additional workflow input parameters as keyword arguments."""
        calc_id = calc_class.id

        calc_type = f'{calc_id}_{spin.value}' if spin != Spin.NONE else calc_id

        if calc_class.input_model is None:
            raise ValueError(
                f'Calculator class {calc_class.__name__} does not have an input model defined. Please report this bug.')

        parameters = getattr(self.inputs.calculator_parameters, calc_type).model_copy(deep=True)
        parameters_class = calc_class.input_model.model_fields['parameters'].annotation

        if not type(parameters) is parameters_class:
            dct = parameters.model_dump()
            dct.update(override)
            parameters = parameters_class.model_validate(dct)
        else:
            for k, v in override.items():
                set_(parameters, k, v)

        # For the k-points, the Workflow has two options: self.kpoints.grid and self.kpoints.path. A calculator should
        # only ever have one of these two. By default, use the grid.

        # # Add further information to the calculator as required
        # for kw in ['gamma_only', 'kgrid', 'kpath', 'koffset', 'plotting']:
        #     if kw not in override and kw in parameters.model_fields:
        #         val: Any
        #         if kw == 'kgrid':
        #             val = self.inputs.kpoints.grid
        #         elif kw == 'kpath':
        #             val = self.inputs.kpoints.path
        #         elif kw == 'koffset':
        #             if (calc_class == PWCalculator and parameters.control.calculation == 'nscf'
        #                     or calc_class == Wannier90Calculator) and self.inputs.kpoints.offset_nscf is not None:
        #                 val = self.inputs.kpoints.offset_nscf
        #             else:
        #                 val = self.inputs.kpoints.offset
        #         elif kw == 'gamma_only':
        #             val = self.inputs.kpoints.gamma_only
        #         else:
        #             val = getattr(self, kw)
        #         kwargs[kw] = val

        if 'pseudopotentials' in calc_class.input_model.model_fields:
            kwargs['pseudopotentials'] = OrderedDict(
                **{k: v.filename.name for k, v in self.inputs.pseudopotentials.items() if v.filename is not None})
        if 'pseudo_family' in calc_class.input_model.model_fields:
            kwargs['pseudo_family'] = self.inputs.config.pseudo_library
        if 'kpoints' in calc_class.input_model.model_fields:
            kwargs['kpoints'] = copy.deepcopy(self.inputs.kpoints)

        # Create the calculator
        calc = calc_class(name=name, parent_process_uid=self.uid, inputs={
                          'atoms': copy.deepcopy(self.inputs.atoms), 'parameters': parameters, **kwargs}, engine=self.engine)

        return calc

    def primitive_to_supercell(self, matrix: Optional[npt.NDArray[np.int_]] = None, **kwargs):
        """Convert from a primitive cell to a supercell, as given by a 3x3 transformation matrix."""
        if matrix is None:
            if self.inputs.kpoints.gamma_only:
                matrix = np.identity(3, dtype=float)
            else:
                assert self.inputs.kpoints.grid is not None
                matrix = np.diag(self.inputs.kpoints.grid)
        assert np.shape(matrix) == (3, 3)
        self.atoms = make_supercell(self.inputs.atoms, matrix, **kwargs)

    def supercell_to_primitive(self, matrix: Optional[npt.NDArray[np.int_]] = None):
        """Convert from a supercell to a primitive cell, as given by a 3x3 transformation matrix.

        This is the inverse of self.primitive_to_supercell()
        """
        if matrix is None:
            assert self.kpoints.grid is not None
            matrix = np.diag(self.kpoints.grid)
        assert np.shape(matrix) == (3, 3)

        # # Work out the atoms belonging to the primitive cell
        inv_matrix = np.linalg.inv(matrix)
        self.atoms.cell = np.dot(inv_matrix, self.atoms.cell)

        # Wrap the atomic positions
        wrapped_atoms = copy.deepcopy(self.atoms)
        wrapped_atoms.wrap(pbc=True)

        # Find all atoms whose positions have not moved
        mask = [np.linalg.norm(a.position - wrapped_a.position) < 1e-7 for a,
                wrapped_a in zip(self.atoms, wrapped_atoms)]

        self.atoms = self.atoms[mask]

    def run_processes(self, processes: Process[Any, Any] | Sequence[Process[Any, Any]], **kwargs) -> Status:
        """Run a sequence of processes in the workflow.

        Run either a process or sequence (i.e. list) of processes. If a list is provided, the processes must be able to be run
        in parallel.

        :param processes: The process or processes to run
        :type processes: ProcessProtocol | Sequence[ProcessProtocol]
        :return: The status of this set of processes
        :rtype: Status
        """
        assert self.engine is not None

        logger = logging.getLogger(__name__)

        # Convert processes to an immutable, covariant tuple
        if isinstance(processes, Process):
            processes = (processes,)
        if not isinstance(processes, tuple):
            processes = tuple(processes)

        for process in processes:
            process.parent_process_uid = self.uid
            process.engine = self.engine
            self.step_counter += 1

            # Check that calc.directory was not set
            if process.directory is not None:
                raise ValueError(
                    f'`{process.name}.directory` should not be set manually, but it was set to `{process.directory}`')

            if self.directory is None:
                process.directory = f'{self.step_counter:02d}-{process.name}'
            else:
                process.directory = self.directory / f'{self.step_counter:02d}-{process.name}'

            status = self.engine.get_status(process)
            if status == Status.NOT_STARTED:
                logger.info(f'Submitting {process.directory} ({process.__class__.__name__})')
                self.engine.run(process, **kwargs)
            else:
                logger.info(f'Step {process.uid} is {status.name}', stacklevel=2)

        for process in processes:
            if self.engine.get_status(process) == Status.COMPLETED:
                # Load the results of the step
                logger.info(f'Loading the results for {process.directory} ({process.__class__.__name__})')
                self.engine.load_results(process)

        # If any of the processes are running, return RUNNING
        for process in processes:
            if self.engine.get_status(process) == Status.RUNNING:
                return Status.RUNNING

        # Otherwise, return COMPLETED
        return Status.COMPLETED

    def processes_are_running(self) -> bool:
        """Return True if any processes have status RUNNING."""
        assert self.engine is not None
        if self.engine.processes_are_running():
            return True
        return False

    def print(self, *args, **kwargs) -> None:
        """Print information in a tidy way."""
        utils.indented_print(*args, **kwargs)

    def todict(self):
        """Convert the workflow to a dictionary."""
        # Shallow copy
        dct = dict(self.__dict__)

        # Adding information required by the json decoder
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @classmethod
    def fromdict(cls, dct: Dict[str, Any], **kwargs) -> Workflow:
        """Create a workflow from a dictionary."""
        # Remove __koopmans_name/module__ if present (won't happen if the encoder was used, but will happen if
        # todict and fromdict are used directly)
        dct.pop('__koopmans_name__', None)
        dct.pop('__koopmans_module__', None)

        wf = cls(atoms=dct.pop('atoms'),
                 parameters=dct.pop('parameters'),
                 calculator_parameters=dct.pop('calculator_parameters'),
                 pseudopotentials=dct.pop('pseudopotentials'),
                 kpoints=dct.pop('kpoints'),
                 projections=dct.pop('projections'),
                 autogenerate_settings=False,
                 **kwargs)

        for k, v in dct.items():
            setattr(wf, k, v)

        return wf

    @property
    def variational_orbitals(self) -> VariationalOrbitals | None:
        """The variational_orbitals object associated with this workflow."""
        return self._variational_orbitals

    @variational_orbitals.setter
    def variational_orbitals(self, value: VariationalOrbitals):
        assert isinstance(value, VariationalOrbitals)
        self._variational_orbitals = value

    @classmethod
    def fromjson(cls, fname: str, override: Dict[str, Any] = {}, **kwargs):
        """Generate a workflow from a json file."""
        with open(fname, 'r') as fd:
            bigdct = json_ext.loads(fd.read())
        wf = cls._fromjsondct(bigdct, override, **kwargs)

        # Define the name of the workflow using the name of the json file
        wf.name = Path(fname).stem

        return wf

    @classmethod
    def _fromjsondct(cls, bigdct: Dict[str, Any], override: Dict[str, Any] = {}, **kwargs):
        """Generate a workflow from a dictionary loaded from a json file."""
        # Override all keywords provided explicitly
        utils.update_nested_dict(bigdct, override)

        # Loading atoms object
        atoms_dict = bigdct.pop('atoms', None)
        if atoms_dict:
            atoms, snapshots = read_atoms_dict(utils.parse_dict(atoms_dict))
        else:
            raise ValueError('Please provide an `atoms` block in the json input file')

        # Loading kpoints
        kpts = Kpoints(**utils.parse_dict(bigdct.pop('kpoints', {})), cell=atoms.cell)

        # Engine
        if 'engine' not in kwargs:
            engine = LocalASEEngine(**bigdct.pop('engine', {}))
        else:
            engine = kwargs.pop('engine')
            engine_parameters = bigdct.pop('engine', {})
            for k, v in engine_parameters.items():
                setattr(engine, k, v)

        # Create the workflow. Note that any keywords provided in the calculator_parameters (i.e. whatever is left in
        # calcdict) are provided as kwargs

        wf = cls(inputs={'atoms': atoms, 'snapshots': snapshots, 'kpoints': kpts, **bigdct}, engine=engine, **kwargs)

        if wf.config.converge:
            # Wrap the workflow in a convergence outer loop
            from koopmans.workflows import ConvergenceWorkflowFactory
            return ConvergenceWorkflowFactory(wf, **conv_block)
        else:
            return wf

    def print_bib(self):
        """Generate a .bib file containing the references relevant to this workflow."""
        relevant_references = BibliographyData()

        def add_ref(bibkey: str, note: str):
            if bibkey not in bib_data.entries:
                raise ValueError(f'Could not find bibliography entry for `{bibkey}`')
            else:
                entry = bib_data.entries[bibkey]
                entry.fields['note'] = note
                relevant_references.add_entry(bibkey, entry)

        add_ref('Linscott2023', 'Introduces the koopmans code')
        if self.inputs.config.correction in [Correction.KI, Correction.KIPZ, Correction.PKIPZ, Correction.ALL]:
            add_ref('Dabo2010', 'One of the founding Koopmans functionals papers')
            add_ref('Borghi2014', 'One of the founding Koopmans functionals papers')

            if all(self.inputs.atoms.pbc):
                add_ref('Nguyen2018', 'Describes Koopmans functionals in periodic systems')
                if self.inputs.config.calculate_alpha:
                    if self.inputs.config.screening_method == 'dfpt':
                        add_ref('Colonna2019', 'Introduces the DFPT method for calculating screening parameters')
                        add_ref('Colonna2022', 'Describes the algorithms underpinning the kcw.x code')
                    else:
                        add_ref('DeGennaro2022', 'Describes how to extract band structures from Koopmans functional '
                                'calculations')
                        add_ref('Borghi2015', 'Describes the algorithms underpinning the kcp.x code')
            else:
                add_ref('Borghi2015', 'Describes the algorithms underpinning the kcp.x code')

            for citation in pseudopotential_library_citations(self.inputs.config.pseudo_library):
                add_ref(citation,
                        f'Citation for the {self.inputs.config.pseudo_library.replace("_", " ")} pseudopotential library')

        utils.print_alert(
            'note', f'Please cite the papers listed in `{self.name}.bib` in work involving this calculation')
        relevant_references.to_file(self.name + '.bib')

    def print_preamble(self):
        """Print the preamble of the workflow."""
        if self.parent_process_uid:
            return

        self.print_bib()

    def toinputjson(self) -> Dict[str, Dict[str, Any]]:
        """Convert the workflow to a dictionary that can be written to a json file."""
        bigdct: Dict[str, Dict[str, Any]] = {}

        bigdct['workflow'] = {}

        # "workflow" block (not printing any values that match the defaults except for core keywords)
        for k, v in self.parameters.items():
            if v is None:
                continue
            if isinstance(v, Path):
                v = str(v)
            if k == 'pseudo_directory' and self.parameters.pseudo_library is not None:
                continue
            default = self.parameters.defaults.get(k, None)
            if v != default or k in ['task', 'functional']:
                bigdct['workflow'][k] = v

        # "atoms" block
        # Working out ibrav
        ibrav = self.calculator_parameters['kcp'].get('ibrav', self.calculator_parameters['pw'].get('ibrav', 0))

        bigdct['atoms'] = {}

        # cell parameters
        if ibrav == 0:
            bigdct['atoms']['cell_parameters'] = utils.construct_cell_parameters_block(self.atoms)

        # atomic positions
        if len(self.snapshots) > 1:
            snapshots_file = 'snapshots.xyz'
            ase_write(snapshots_file, self.snapshots)
            bigdct['atoms']['atomic_positions'] = {'snapshots': snapshots_file}
        else:
            if len(set(self.atoms.get_tags())) > 1:
                labels = [s + str(t) if t > 0 else s for s, t in zip(self.atoms.symbols, self.atoms.get_tags())]
            else:
                labels = self.atoms.symbols
            if ibrav == 0:
                bigdct['atoms']['atomic_positions'] = {'positions': [
                    [label] + [str(x) for x in pos] for label, pos in zip(labels, self.atoms.get_positions())],
                    'units': 'angstrom'}
            else:
                bigdct['atoms']['atomic_positions'] = {'positions': [
                    [label] + [str(x) for x in pos] for label, pos in zip(labels, self.atoms.get_scaled_positions())],
                    'units': 'crystal'}

        # k-points
        bigdct['kpoints'] = self.kpoints.tojson()

        # Populating the calculator subdictionary
        bigdct['calculator_parameters'] = {}
        calcdct = bigdct['calculator_parameters']
        calcdct['wannier90'] = {}
        calcdct['ui'] = {}
        for code, params in self.calculator_parameters.items():
            # Remove default settings (ensuring we switch to using relative paths to check this)
            tmp, params.use_relative_paths = params.use_relative_paths, True
            params_dict = {k: v for k, v in params.items() if params.defaults.get(k, None) != v}
            params.use_relative_paths = tmp

            # convert Paths to strings
            for k in params_dict:
                if isinstance(params_dict[k], Path):
                    params_dict[k] = str(params_dict[k])

            # If the params_dict is empty, don't add a block for this calculator
            if not params_dict and not code.startswith('wannier90'):
                continue

            if code in ['pw', 'kcp']:
                calcdct[code] = {}

                # Populate calcdct with the settings
                input_data = construct_namelist(params_dict)
                for key, block in input_data.items():

                    if len(block) > 0:
                        calcdct[code][key] = {k: v for k, v in dict(
                            block).items() if v is not None}

            elif code in ['pw2wannier', 'kcw_wannier', 'kcw_screen', 'kcw_ham', 'projwfc', 'wann2kcp', 'ph']:
                calcdct[code] = params_dict
            elif code.startswith('ui_'):
                calcdct['ui'][code.split('_')[-1]] = params_dict
            elif code == 'ui':
                calcdct['ui'].update(**params_dict)
            elif code.startswith('wannier90'):
                if '_' in code:
                    [spin] = [s for s in Spin if s.value == code.split('_')[1]]
                    calcdct['wannier90'][spin.value] = params_dict
                else:
                    spin = Spin.NONE
                    calcdct['wannier90'] = params_dict
                projections = self.projections.get_subset(spin) if self.projections else None
                if projections:
                    proj_kwarg = {'projections': [p.projections for p in projections]}
                else:
                    proj_kwarg = {}

                if spin != Spin.NONE:
                    calcdct['wannier90'][spin.value].update(**proj_kwarg)
                else:
                    calcdct['wannier90'].update(**proj_kwarg)
            else:
                raise NotImplementedError(
                    f'Writing of `{params.__class__.__name__}` with `write_json` is not yet implemented')

        other_blocks: Dict[str, Any] = {'plotting': self.plotting, 'ml': self.ml}
        for key, params in other_blocks.items():
            dct: Dict[str, Any] = {k: v for k, v in params.items() if params.defaults.get(k, None) != v}
            if dct:
                bigdct[key] = dct

        return bigdct

    def plot_bandstructure(self,
                           bs: Union[BandStructure, List[BandStructure]],
                           dos: Optional[Union[GridDOSCollection, DOS]] = None,
                           filename: Optional[str] = None,
                           bsplot_kwargs: Union[Dict[str, Any], List[Dict[str, Any]]] = {},
                           dosplot_kwargs: Dict[str, Any] = {}) -> None:
        """Plot the provided band structure (and optionally also a provided DOS).

        Arguments
        ---------
        bs -- a bandstructure/list of band structures to be plotted
        dos -- a density of states object to be plotted
        filename -- the name of the file to which the figure will be saved
        bsplot_kwargs -- keyword arguments for when plotting the band structure(s). N.B. if bs is a list of band
                         structures, bsplot_kwargs must be a list of equal length, with kwarg dicts for each entry
        dosplot_kwargs -- keyword arguments for when plotting the DOS
        """
        # Sanitize input
        if isinstance(bs, BandStructure):
            bs = [bs]
        if isinstance(bsplot_kwargs, dict):
            bsplot_kwargs = [bsplot_kwargs for _ in bs]
        if len(bs) != len(bsplot_kwargs):
            raise ValueError('The `bs` and `bsplot_kwargs` arguments to `plot_bandstructure()` should be the same '
                             'length')
        spins: List[Optional[str]]
        if isinstance(dos, DOS):
            if self.parameters.spin_polarized:
                spins = ['up', 'down']
            else:
                spins = [None]
            dos = GridDOSCollection([GridDOSData(dos.get_energies(), dos.get_dos(ispin), info={'spin': spin})
                                     for ispin, spin in enumerate(spins)])

        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        if dos is not None:
            # Construct the axes
            _, axes = plt.subplots(1, 2, sharey=True, gridspec_kw={'width_ratios': [3, 1]})
            ax_bs = axes[0]
            ax_dos = axes[1]
        else:
            ax_bs = None
            ax_dos = None

        # Plot the band structure
        defaults = {'colors': colors, 'emin': self.inputs.plotting_config.emin,
                    'emax': self.inputs.plotting_config.emax}
        for b, kwargs in zip(bs, bsplot_kwargs):
            for k, v in defaults.items():
                if k not in kwargs:
                    kwargs[k] = v
            ax_bs = b.plot(ax=ax_bs, **kwargs)

        # Move the legend (if there is one)
        _, labels = ax_bs.get_legend_handles_labels()
        if len(labels) > 0:
            ax_bs.legend(loc='lower left', bbox_to_anchor=(0, 1), frameon=False, ncol=min((2, len(bs))))

        if dos is None:
            axes = [ax_bs]
        else:
            # Assemble the densities of state
            dos_summed = dos.sum_by('symbol', 'n', 'l', 'spin')
            if self.inputs.config.spin_polarized:
                dos_up = dos_summed.select(spin='up')
                dos_down = dos_summed.select(spin='down')
                dos_down._weights *= -1
                doss = [dos_up, dos_down]
            else:
                doss = [dos_summed]

            # Plot the DOSs
            spdf_order = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
            [xmin, xmax] = ax_bs.get_ylim()
            for dos in doss:
                # Before iterating through the DOSs, sort them in a sensible order (by atomic symbol, n, and l)
                if all([key in d.info for d in dos for key in ['symbol', 'n', 'l']]):
                    sorted_dos = sorted(dos, key=lambda x: (x.info['symbol'], x.info['n'], spdf_order[x.info['l']]))
                else:
                    sorted_dos = dos

                for d in sorted_dos:
                    if (not self.inputs.config.spin_polarized or d.info.get('spin') == 'up') \
                            and all([key in d.info for key in ['symbol', 'n', 'l']]):
                        label = f'{d.info["symbol"]} {d.info["n"]}{d.info["l"]}'
                    else:
                        label = None
                    d.plot(ax=ax_dos, xmin=xmin, xmax=xmax, vertical=True, mplargs={'label': label},
                           **dosplot_kwargs)

                # Reset color cycle so the colors of spin-up match those of spin-down
                ax_dos.set_prop_cycle(None)

            # Tweaking the DOS figure aesthetics
            maxval = 1.1 * dos_summed._weights[:, [e >= xmin and e <= xmax for e in dos_summed._energies]].max()
            if self.inputs.config.spin_polarized:
                ax_dos.set_xlim([maxval, -maxval])
                ax_dos.text(0.25, 0.10, 'up', ha='center', va='top', transform=ax_dos.transAxes)
                ax_dos.text(0.75, 0.10, 'down', ha='center', va='top', transform=ax_dos.transAxes)
            else:
                ax_dos.set_xlim([0, maxval])
            ax_dos.set_xticks([])

            # Move the legend (if there is one)
            _, labels = ax_dos.get_legend_handles_labels()
            if len(labels) > 0:
                ax_dos.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
            plt.subplots_adjust(right=0.85, wspace=0.05)

        # Saving the figure to file (as png and also in editable form)
        workflow_name = self.name
        workflow_name = workflow_name.replace('workflow', '').replace(' ', '_')
        filename = filename if filename is not None else f'{workflow_name}_bandstructure'
        filepath = self.directory / filename if self.directory is not None else filename
        legends = [ax.get_legend() for ax in axes if ax.get_legend() is not None]
        utils.savefig(fname=f'{filepath}.png', bbox_extra_artists=legends, bbox_inches='tight')

        # Close the figure
        plt.close()

    def _teardown(self):
        """Perform final tasks before the workflow completes."""
        from koopmans.io import write

        assert self.status == Status.COMPLETED
        assert self.engine is not None

        # Save workflow to file
        write(self, self.name + '.pkl')

        # Save the ML model to a separate file
        if self.inputs.ml_config.train:
            assert self.inputs.ml_model is not None
            write(self.inputs.ml_model, self.name + '_ml_model.pkl')

        # Removing tmpdirs
        if not self.engine.keep_tmpdirs:
            self.engine.remove_tmpdirs()

        # Engine-dependent teardown
        self.engine._teardown()

    def number_of_electrons(self, spin: Spin = Spin.NONE, params: PWInputParameters | KCPInputParameters | None = None) -> int:
        if params is None:
            params = self.inputs.calculator_parameters.pw
        return _number_of_electrons(self.inputs.atoms, self.inputs.pseudopotentials, self.inputs.config.spin_polarized, params, spin)


def read_atoms_dict(dct: Dict[str, Any]) -> Tuple[Atoms, List[Atoms]]:
    """Read the "atoms" block."""
    subdct: Dict[str, Any]
    if 'snapshots' in dct.get('atomic_positions', {}):
        subdct = dct.pop('atomic_positions')
        snapshots = ase_read(subdct['snapshots'], index=':')

        atoms = snapshots[0]

        # If cell parameters are manually provided, override what is provided in the xyz file
        if 'cell_parameters' in dct:
            utils.read_cell_parameters(atoms, dct['cell_parameters'])
            for s in snapshots[1:]:
                s.cell = atoms.cell
                s.pbc = atoms.pbc
    else:
        atoms = Atoms()
        snapshots = None

        readers: Dict[str, Callable] = {'cell_parameters': utils.read_cell_parameters,
                                        'atomic_positions': utils.read_atomic_positions}

        for key, reader in readers.items():
            subdct = dct.pop(key, {})
            if subdct:
                reader(atoms, subdct)
            else:
                raise ValueError(f'Please provide `{key}` in the `atoms` block')

        for block in dct:
            raise ValueError(f'Unrecognized subblock `atoms: {block}`')

    return atoms, snapshots


def spin_symmetrize(wf: Workflow, master_calc: calculators.Calc) -> Status:
    """Enforce the spin-symmetry of a spin-polarized calculation.

    Generates a series of calculators that will be equivalent to the original master_calc
    but switching back and forth between nspin=1 and nspin=2 calculations to fix spin contamination
    """
    if not isinstance(master_calc, calculators.CalculatorCanEnforceSpinSym):
        raise NotImplementedError(f'`{master_calc.__class__.__name__}` cannot enforce spin symmetry')

    # nspin=1
    calc_nspin1 = master_calc.nspin1_calculator()

    if not master_calc.from_scratch:
        # Locate the preceeding nspin=2 calculation from which the calculation is configured to read
        # (and from which we will fetch the nspin=2 wavefunctions)
        ndr_directory = str(master_calc.read_directory)
        if ndr_directory not in master_calc.linked_files.keys():
            raise ValueError('This calculator needs to be linked to a previous `nspin=2` calculation before '
                             '`run_calculator` is called')
        prev_calc_nspin2 = master_calc.linked_files[ndr_directory][0].parent_process
        if not isinstance(prev_calc_nspin2, calculators.KoopmansCPCalculator):
            raise ValueError(
                'Unexpected calculator type linked during the procedure for fixing spin contamination')

        # Wipe the linked files (except for globally-linked files)
        master_calc.linked_files = {k: v for k, v in master_calc.linked_files.items(
        ) if isinstance(v[0].parent_process, ProcessPlaceholder)}
        calc_nspin1.linked_files = {k: v for k, v in calc_nspin1.linked_files.items(
        ) if isinstance(v[0].parent_process, ProcessPlaceholder)}

        # nspin=1 dummy
        calc_nspin1_dummy = master_calc.nspin1_dummy_calculator()
        calc_nspin1_dummy.skip_qc = True
        status = wf.run_processes(calc_nspin1_dummy)
        if status != Status.COMPLETED:
            return status

        calc_nspin1.link(File(calc_nspin1_dummy, calc_nspin1_dummy.parameters.outdir),
                         calc_nspin1.parameters.outdir)

        # Copy over nspin=2 wavefunction to nspin=1 tmp directory
        process2to1 = ConvertFilesFromSpin2To1(name=None,
                                               **prev_calc_nspin2.files_to_convert_with_spin2_to_spin1)
        status = wf.run_processes(process2to1)
        if status != Status.COMPLETED:
            return status

        # Linking converted files to the nspin=1 calculation
        for f in process2to1.outputs.generated_files:
            dst_dir = calc_nspin1.parameters.outdir / \
                f'{calc_nspin1.parameters.prefix}_{calc_nspin1.parameters.ndr}.save/K00001'
            calc_nspin1.link(f, dst_dir / f.relative_path, symlink=True, overwrite=True)

    status = wf.run_processes(calc_nspin1)
    if status != Status.COMPLETED:
        return status

    # nspin=2 from scratch (dummy run for creating files of appropriate size)
    calc_nspin2_dummy = master_calc.nspin2_dummy_calculator()
    calc_nspin2_dummy.skip_qc = True
    status = wf.run_processes(calc_nspin2_dummy)
    if status != Status.COMPLETED:
        return status
    master_calc.link(File(calc_nspin2_dummy, calc_nspin2_dummy.parameters.outdir), master_calc.parameters.outdir)

    # Copy over nspin=1 wavefunction to nspin=2 tmp directory
    process1to2 = ConvertFilesFromSpin1To2(**calc_nspin1.files_to_convert_with_spin1_to_spin2)
    status = wf.run_processes(process1to2)
    if status != Status.COMPLETED:
        return status

    # nspin=2, reading in the spin-symmetric nspin=1 wavefunction
    master_calc.prepare_to_read_nspin1()
    for f in process1to2.outputs.generated_files:
        dst_dir = master_calc.parameters.outdir / \
            f'{master_calc.parameters.prefix}_{master_calc.parameters.ndr}.save/K00001'
        master_calc.link(f, dst_dir / f.relative_path, symlink=True, overwrite=True)
    status = wf.run_processes(master_calc)
    return status


def _number_of_electrons(atoms: Atoms, pseudopotentials: OrderedDict[str, UPFModel],
                         spin_polarized: bool, params: PWInputParameters | KCPInputParameters,
                         spin: Spin = Spin.NONE) -> int:
    """Return the number of electrons in a particular spin channel."""
    nelec_tot = nelec_from_pseudos(atoms, pseudopotentials)
    nelec = nelec_tot - int(params.system.tot_charge)
    if spin_polarized:
        mag = int(params.system.tot_magnetization) if params.system.tot_magnetization is not None else 0
        if spin == Spin.UP:
            nelec += mag
            nelec //= 2
        elif spin == Spin.DOWN:
            nelec -= mag
            nelec //= 2
    else:
        if spin != Spin.NONE:
            nelec //= 2
    return nelec
