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
from pydash import get
from typing_extensions import Self

# isort: off
import koopmans.mpl_config  # noqa: F401
import matplotlib.pyplot as plt
# isort: on

from ase_koopmans import Atoms
from ase_koopmans.build.supercells import make_supercell
from ase_koopmans.cell import Cell
from ase_koopmans.dft.dos import DOS
from ase_koopmans.dft.kpoints import BandPath
from ase_koopmans.io import read as ase_read
from ase_koopmans.io import write as ase_write
from ase_koopmans.io.espresso import \
    contruct_kcp_namelist as construct_namelist
from ase_koopmans.spacegroup import symmetrize
from ase_koopmans.spectrum.band_structure import BandStructure
from ase_koopmans.spectrum.doscollection import GridDOSCollection
from ase_koopmans.spectrum.dosdata import GridDOSData
from upf_tools import UPFDict

from koopmans import calculators, settings, utils
from koopmans.calculators.parameters import CalculatorsInputParameters
from koopmans.engines import Engine, LocalhostEngine
from koopmans.files import File, ParentProcessPlaceholder
from koopmans.kpoints import Kpoints
from koopmans.ml import AbstractMLModel, MLModel, OccEmpMLModels
from koopmans.processes.koopmans_cp import (ConvertFilesFromSpin1To2,
                                            ConvertFilesFromSpin2To1)
from koopmans.processes.process import Process
from koopmans.processes.protocol import ProcessProtocol
from koopmans.projections import (ExplicitProjections, ImplicitProjections,
                                  Projections)
from koopmans.pseudopotentials import (nelec_from_pseudos,
                                       nwfcs_from_projectors,
                                       nwfcs_from_pseudos,
                                       pseudopotential_library_citations)
from koopmans.pydantic import BaseModel
from koopmans.references import bib_data
from koopmans.status import Status
from koopmans.utils import Spin
from koopmans.utils.warnings import warn
from koopmans.variational_orbitals import VariationalOrbitals

from .config import Task, WorkflowConfig


class WorkflowInput(BaseModel):
    atoms: Atoms = Field('an ASE ``Atoms`` object defining the atomic positions, cell, etc')
    pseudopotentials: OrderedDict[str, UPFDict] = Field('a dictionary mapping atom labels to pseudopotential objects')
    kpoints: Kpoints
    projections: ExplicitProjections | ImplicitProjections | None = None
    calculator_parameters: CalculatorsInputParameters = Field(default_factory=CalculatorsInputParameters)
    plotting_config: settings.PlotSettingsDict = Field(default_factory=settings.PlotSettingsDict)
    ml_config: settings.MLSettingsDict = Field(default_factory=settings.MLSettingsDict)
    snapshots: List[Atoms]

    @model_validator(mode='before')
    @classmethod
    def setup_implicit_projections(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Set up the projections if they are implicit."""
        calculator_parameters = data['calculator_parameters']
        auto = [calculator_parameters[key].get('auto_projections', False)
                for key in ['w90', 'w90_up', 'w90_down'] if key in calculator_parameters]
        if any(auto):
            pseudopotentials = data['pseudopotentials']
            atoms = data['atoms']
            if data['config'].spin_polarized:
                spins = [Spin.UP, Spin.DOWN]
            else:
                spins = [Spin.NONE]
            if calculator_parameters['pw2wannier'].atom_proj_ext:
                proj_dir = calculator_parameters['pw2wannier'].atom_proj_dir
                if proj_dir is None:
                    raise ValueError('`atom_proj_dir` must be specified when using `atom_proj_ext = True`')
                projectors = OrderedDict([(k, proj_dir / p.filename.with_suffix('.dat').name)
                                         for k, p in pseudopotentials.items()])
                num_wann = [nwfcs_from_projectors(atoms, projectors) for _ in spins]
            else:
                num_wann = [nwfcs_from_pseudos(atoms, pseudopotentials) for _ in spins]
            data['projections'] = ImplicitProjections.from_block_lengths(num_wann, spins, atoms)

            # Also override various pw2wannier and w90 settings
            calculator_parameters['pw2wannier'].atom_proj = True
            for spin in spins:
                label = 'w90' if spin == Spin.NONE else f'w90_{spin.value}'
                calculator_parameters[label].guiding_centres = False
                calculator_parameters[label].dis_froz_proj = True
                if calculator_parameters[label].dis_proj_max is None:
                    calculator_parameters[label].dis_proj_max = 0.95
            data['calculator_parameters'] = calculator_parameters
        return data

    @model_validator(mode='before')
    @classmethod
    def populate_snapshots(cls, values: dict[str, Any]) -> dict[str, Any]:
        if values.get('snapshots', None) is None:
            values['snapshots'] = [values['atoms']]
        return values

    @model_validator(mode="before")
    @classmethod
    def set_default_kpoints(cls, dct: dict[str, Any]) -> dict[str, Any]:
        atoms = dct.get('atoms')
        if not isinstance(atoms, Atoms):
            raise ValueError('`atoms` must be an ASE `Atoms` object to set default k-points')
        if not isinstance(atoms.cell, Cell):
            raise ValueError('`atoms.cell` must be an ASE `Cell` object to set default k-points')
        if all(atoms.pbc):
            # By default, use ASE's default bandpath for this cell (see
            # https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#brillouin-zone-data)
            default_path = atoms.cell.bandpath(eps=1e-10).path
        else:
            default_path = 'G'
        kpoints = dct.get('kpoints', None)
        if kpoints is None:
            kpoints = Kpoints(path=default_path, gamma_only=False, cell=atoms.cell)
        elif kpoints.path is None:
            kpoints.set_path(default_path, cell=atoms.cell)
        dct['kpoints'] = kpoints
        return dct

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
        projections = data.get('projections', None)
        calculator_parameters = data['calculator_parameters']
        if projections:
            for spin in [Spin.UP, Spin.DOWN, Spin.NONE]:
                label = 'w90'
                if spin != Spin.NONE:
                    label += f'_{spin}'
                projections.exclude_bands[spin] = calculator_parameters[label].get('exclude_bands', [])
            data['projections'] = projections
        return data

    @field_validator('atoms', mode='before')
    @classmethod
    def detach_atoms_calculator(cls, value: Atoms) -> Atoms:
        if value.calc is not None:
            warn(f'You have initialized a `{cls.__name__}` object with an atoms object that '
                 'possesses a calculator. This calculator will be ignored.')
            value.calc = None
        return value

    @field_validator('calculator_parameters', mode='after')
    @classmethod
    def check_projections_provided_in_right_place(cls, value: Dict[str, settings.SettingsDict]) -> Dict[str, settings.SettingsDict]:
        # Various checks for the wannier90 blocks
        for block, params in value.items():
            if block.startswith('w90'):
                if 'projections' in params:
                    raise ValueError(f'You have provided projection information in the '
                                     f'`calculator_parameters[{block}]` argument to `{cls.__name__}`. '
                                     'Please instead specify projections via the `projections` argument')
        return value

    @field_validator('calculator_parameters', mode='after')
    @classmethod
    def warn_about_overwritten_keywords(cls, value: Dict[str, settings.SettingsDict]) -> Dict[str, settings.SettingsDict]:
        for block, params in value.items():
            if block.startswith('w90'):
                for kw in ['exclude_bands', 'num_wann', 'num_bands', 'projections']:
                    if kw in params:
                        warn(f'`{kw}` will be overwritten by the workflow; it is best to not specify this '
                             'keyword and to instead double-check the keyword in the various `.win` files '
                             'generated by the workflow.')
        return value


OutputModel = TypeVar('OutputModel', bound=BaseModel)


class Workflow(Process[WorkflowInput, OutputModel], Generic[OutputModel]):
    """Abstract base class that defines a Koopmans workflow."""

    input_model = WorkflowInput
    config: WorkflowConfig
    version: str = Field(
        'Version number used to create this workflow (becomes relevant when loading workflows from .kwf files)')
    step_counter: int = 0
    print_indent: int = 0
    status: Status = Status.NOT_STARTED
    parent_process: Process[Any, Any] | None = None
    calculations: List[calculators.CalculatorExt] = []
    processes: List[ProcessProtocol] = []
    steps: List = []
    _variational_orbitals: VariationalOrbitals | None
    ml_model: Optional[AbstractMLModel] = None

    @model_validator(mode='before')
    @classmethod
    def create_ml_model(cls, data: dict[str, Any]) -> dict[str, Any]:
        ml_config = data['inputs'].ml_config
        ml_model = data['inputs'].get('ml_model', None)
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
            data['ml_model'] = ml_model
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
        nelec = _number_of_electrons(atoms, pseudopotentials, spin_polarized, Spin.NONE, kcp_params)
        nelup = _number_of_electrons(atoms, pseudopotentials, spin_polarized, Spin.UP, kcp_params)
        neldw = _number_of_electrons(atoms, pseudopotentials, spin_polarized, Spin.DOWN, kcp_params)
        tot_charge = nelec_from_pseudos(atoms, pseudopotentials) - nelec
        tot_mag = nelup - neldw

        # Setting up the magnetic moments
        starting_mags = {k: v for k, v in calculator_parameters.get('kcp', {}).items(
        ) if k.startswith('starting_magnetization')}
        starting_mags.update(
            {k: v for k, v in calculator_parameters['pw'].items() if k.startswith('starting_magnetization')})
        if starting_mags:
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
        data['inputs']['atoms'] = atoms

        # Work out the number of bands
        nbnd = get(calculator_parameters, 'pw.nbnd', get(calculator_parameters, 'kcp.nbnd', nelec // 2 + nelec % 2))
        generated_keywords = {'nelec': nelec, 'tot_charge': tot_charge, 'tot_magnetization': tot_mag,
                              'nelup': nelup, 'neldw': neldw, 'nbnd': nbnd}

        for block, params in get(data, 'inputs.calculator_parameters', {}).items():
            # Apply auto-generated keywords
            for k, v in generated_keywords.items():
                # Skipping nbnd for kcp; if necessary it will be generated later
                if k == 'nbnd' and block == 'kcp':
                    continue
                if k in params.model_fields:
                    params[k] = v

        return data

    @model_validator(mode='before')
    @classmethod
    def distribute_generic_calculator_parameters(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Distribute generic calculator parameters (e.g. ecutwfc) to all relevant calculator blocks."""
        extras = {}
        for k, v in get(data, 'inputs.calculator_parameters', {}).items():
            if k not in CalculatorsInputParameters.model_fields:
                extras[k] = data['inputs']['calculator_parameters'].pop(k)
        data['inputs']['calculator_parameters']['extras'] = extras
        raise NotImplementedError(
            "This implementation won't work because we need these keywords to be correctly distributed before the next before validator")

    @model_validator(mode='before')
    @classmethod
    def get_pseudopotentials(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Get the pseudopotentials, ensuring that they are coerced into UPFDict objects.

        Note this is done inside `Workflow` and not `WorkflowInput` because it requires access
        to the engine and config."""
        engine = data['engine']
        atoms = get(data, 'inputs.atoms')
        pseudopotentials = get(data, 'inputs.pseudopotentials', None)
        pseudo_library = get(data, 'config.pseudo_library')

        upfs = OrderedDict()
        if pseudopotentials is None:
            # Work out the unique atom labels. This needs to be done like this rather than with set() to ensure
            # the order is consistent with how the atoms appear in the input file
            unique_labels = []
            for a in atoms:
                label_tuple: tuple[str, int] = (a.symbol, a.tag)
                if label_tuple not in unique_labels:
                    unique_labels.append(label_tuple)

            for element, tag in unique_labels:
                pseudo = engine.get_pseudopotential(library=pseudo_library, element=element)

                pseudo_type = pseudo['header']['pseudo_type']
                if pseudo_type != 'NC':
                    raise ValueError('Koopmans functionals only currently support norm-conserving pseudopotentials; '
                                     f'{pseudo.filename} is {pseudo_type}')

                symbol = element + str(tag) if tag > 0 else element
                upfs[symbol] = pseudo
        else:
            for symbol, pseudo_filename in sorted(pseudopotentials.items()):
                if isinstance(pseudo_filename, str):
                    upf = engine.get_pseudopotential(pseudo_library, filename=pseudo_filename)
                elif isinstance(pseudo_filename, UPFDict):
                    upf = pseudo_filename
                else:
                    raise ValueError(f'Invalid pseudopotential type: {pseudo_filename.__class__.__name__}')
                upfs[symbol] = upf
        data['inputs']['pseudopotentials'] = upfs
        return data

    @model_validator(mode='after')
    def ensure_w90_config_and_spin_polarization_are_consistent(self) -> Self:
        if self.config.spin_polarized != (self.inputs.calculator_parameters.w90 is None):
            raise ValueError("Inconsistent spin polarization settings: for spin-polarized calculations, "
                             "both 'w90_up' and 'w90_down' must be specified; for non-spin-polarized calculations, "
                             "'w90' must be specified.")
        return self

    @model_validator(mode='before')
    @classmethod
    def extract_projections(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Adding the projections to the workflow kwargs.

        This is unusual in that this is an attribute of the workflow but can be specified within the w90 block
        """
        calc_configs = get(data, 'inputs.calculator_parameters')
        atoms = get(data, 'inputs.atoms')
        projections = get(data, 'inputs.projections', None)
        spin_polarized = get(data, 'config.spin_polarized', False)

        if projections is None:
            w90_block_projs: list[dict[str, Any]] = []

            spins = [Spin.UP, Spin.DOWN] if spin_polarized else [Spin.NONE]
            suffixes = ['up', 'down'] if spin_polarized else ['']
            for spin, suffix in zip(spins, suffixes):
                projs = get(data, f'inputs.w90{suffix}.projections', [[]])

                # Add the spin information
                w90_block_projs += [{'spin': spin, 'projections': proj_block} for proj_block in projs]

            # Suppress type-checking because mypy does not understand that pydantic will cast blocks appropriately
            # to the desired type
            data['projections'] = ExplicitProjections(blocks=w90_block_projs, atoms=atoms)  # type: ignore[arg-type]
        return data

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

    def __repr__(self):
        entries = []

        # atoms
        entries.append(f'atoms={self.atoms.symbols}')

        # parameters
        entries.append(f'parameters={self.parameters.briefrepr()}')

        # kpoints
        entries.append(f'kpoints={self.kpoints}')

        # pseudopotentials
        entries.append(f'pseudopotentials={self.pseudopotentials}')
        return f'{self.__class__.__name__}(' + ',\n   '.join(entries) + ')'

    @utils.HasDirectory.base_directory.setter  # type: ignore
    def base_directory(self, value: Path):
        """Set the base_directory attribute with additional workflow-specific behavior."""
        old_base_directory = getattr(self, '_base_directory', None)

        super(Workflow, Workflow).base_directory.__set__(self, value)

        # If the pseudo_directory has been set, we need to update it too (because it is stored relative to
        # base_directory). Note that during workflow initialization, we set base_directory prior to defining
        # self.parameters (in which case we don't need to update pseudo_directory)
        if hasattr(self, 'parameters') and self.parameters.pseudo_directory is not None \
                and old_base_directory is not None:
            abs_pseudo_dir = (old_base_directory / self.parameters.pseudo_directory).resolve()
            assert abs_pseudo_dir.is_dir()
            self.parameters.pseudo_directory = os.path.relpath(abs_pseudo_dir, self.base_directory)

    def run(self):
        """Run the workflow."""
        assert self.engine is not None

        self.print_preamble()

        if not self.parent_process:
            bf = '**' if sys.stdout.isatty() else ''
            self.print(bf + self.name + bf)
            self.print(bf + '-' * len(self.name) + bf)
            if self.engine.from_scratch:
                self._remove_tmpdirs()
            self._run_sanity_checks()

        start_time = time.time()

        while self.status != Status.COMPLETED:
            logger = logging.getLogger(__name__)
            logger.info(f'Proceeding with workflow {self.name} ({self.__class__.__name__})')

            self.proceed()

            while self.engine.steps_are_running():
                logger.info('Steps are running; waiting for them to complete...')
                self.wait()
                if self.parameters.max_time is not None and time.time() - start_time > self.parameters.max_time:
                    self.status = Status.FAILED
                    raise ValueError('Workflow failed to complete within the specified time')

        if not self.parent_process and self.status == Status.COMPLETED:
            self._teardown()

        self.status = Status.COMPLETED

    def wait(self):
        """Wait while calculations are running.

        While waiting, the engine might perform other actions such as fetching remote files etc.
        """
        self.engine.wait(seconds=self.parameters.wait_time)

    def _pre_run(self):
        """Run checks and actions before running the workflow."""
        self.step_counter = 0

        if self.status == Status.NOT_STARTED:
            self.status = Status.RUNNING

    def proceed(self, subdirectory: Optional[str] = None, copy_outputs_to_parent: Optional[bool] = True):
        """Proceed through the workflow until a step has status `RUNNING`.

        Typically, this will involve proceeding to the next step in the workflow which has not been run,
        set it `RUNNING`, and then returning

        Note, however, for some engines this function will execute the _entire_ workflow e.g. `localhost` workflows,
        which, when requested to run a step in the workflow, will run that step and block the interpreter until the
        step is completed. For such a workflow, the status of individual steps proceed immediately from `NOT_STARTED`
        to `COMPLETE` and `proceed()` will not exit until all the steps in the workflow are `COMPLETED`.
        """
        logger = logging.getLogger(__name__)
        logger.info(f'Proceeding to the next step in the workflow {self.name} ({self.__class__.__name__})')

        self._pre_run()

        if self.parent_process:
            with self._parent_context(subdirectory, copy_outputs_to_parent):
                self._run()
        else:
            self.base_directory = Path.cwd().resolve()
            if subdirectory:
                self.directory = Path(subdirectory)
            self._run()

        self._post_run()
        logger.info(f'Exiting workflow {self.name} with status {self.status.name} ({self.__class__.__name__})')

    def _post_run(self):
        """Run checks and actions after running the workflow."""
        assert self.engine is not None
        self.engine.update_statuses()

    @abstractmethod
    def _run(self) -> None:
        """Run the workflow; to be implemented by subclasses."""
        ...

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
    def from_other(cls: Type[W], other_wf: Workflow, **kwargs: Any) -> W:
        """Create a new workflow from another workflow, copying all settings and parameters.

        e.g.
        >>> new_wf = Workflow.fromother(other_wf)
        """
        parameters = copy.deepcopy(other_wf.parameters)

        # Pass the pseudo_directory as an absolute path
        if 'pseudo_directory' in parameters:
            parameters['pseudo_directory'] = other_wf.base_directory / parameters['pseudo_directory']

        for k in list(kwargs):
            if parameters.is_valid(k):
                parameters[k] = kwargs.pop(k)
        kwargs.update(**dict(atoms=copy.deepcopy(other_wf.atoms),
                             parameters=parameters,
                             calculator_parameters=copy.deepcopy(other_wf.calculator_parameters),
                             pseudopotentials=copy.copy(other_wf.pseudopotentials),
                             kpoints=copy.deepcopy(other_wf.kpoints),
                             projections=other_wf.projections,
                             plotting=copy.deepcopy(other_wf.plotting),
                             ml=copy.deepcopy(other_wf.ml),
                             ml_model=other_wf.ml_model,
                             engine=other_wf.engine,),
                      )

        return cls(**kwargs)

    @classmethod
    def fromparent(cls: Type[W], parent_wf: Workflow, **kwargs: Any) -> W:
        """Create a subworkflow with the same configuration as the parent workflow.

        e.g.
        >>> sub_wf = Workflow.fromparent(self)
        """
        return cls.from_other(other_wf=parent_wf, parent_process=parent_wf, **kwargs)

    def _run_sanity_checks(self):
        # Check internal consistency of workflow settings
        if self.parameters.fix_spin_contamination is None:
            self.parameters.fix_spin_contamination = not self.parameters.spin_polarized
        else:
            if self.parameters.fix_spin_contamination and self.parameters.spin_polarized:
                raise ValueError('`fix_spin_contamination = True` is incompatible with `spin_polarized = True`')

        if self.parameters.method == 'dfpt':
            if self.parameters.frozen_orbitals is None:
                self.parameters.frozen_orbitals = True
            if not self.parameters.frozen_orbitals:
                raise ValueError('`frozen_orbitals` must be equal to `True` when `method` is `dfpt`')
        else:
            if self.parameters.frozen_orbitals is None:
                self.parameters.frozen_orbitals = False
            if self.parameters.frozen_orbitals:
                warn('You have requested a ΔSCF calculation with frozen orbitals. This is unusual; proceed '
                     'only if you know what you are doing')

        if self.parameters.functional == 'dft':
            self.parameters.calculate_alpha = False

        # Checking periodic image correction schemes
        if not self.parameters.calculate_alpha:
            # If we are not calculating alpha, we do not consider charged systems and therefore we don't need image
            # corrections, so we skip the following checks
            pass
        elif all(self.atoms.pbc):
            if self.parameters.method == 'dfpt':
                # For DPFT, we use gb_correction
                if self.parameters.gb_correction is None:
                    self.parameters.gb_correction = True
                if not self.parameters.gb_correction:
                    warn('Gygi-Baldereschi corrections are not being used; do this with '
                         'caution for periodic systems')

            elif self.parameters.method == 'dscf':
                # For DSCF, we use mp_correction
                if self.parameters.mp_correction is None:
                    self.parameters.mp_correction = True
                if not self.parameters.mp_correction:
                    warn('Makov-Payne corrections are not being used; do this with '
                         'caution for periodic systems')

                if self.parameters.eps_inf is None:
                    warn('`eps_inf` missing in input; it will default to 1.0. Proceed with caution for periodic '
                         'systems; consider setting `eps_inf == "auto"` to calculate it automatically.')
                    self.parameters.eps_inf = 1.0

            if self.parameters.mt_correction is None:
                self.parameters.mt_correction = False
            if self.parameters.mt_correction:
                raise ValueError('Do not use Martyna-Tuckerman corrections for periodic systems')

            # Check the value of eps_inf
            if self.parameters.eps_inf:
                if isinstance(self.parameters.eps_inf, float) and self.parameters.eps_inf < 1.0:
                    raise ValueError('`eps_inf` cannot be lower than 1.0')

            # Check symmetry of the system
            dataset = symmetrize.check_symmetry(self.atoms, 1e-6, verbose=False)
            if dataset is None or dataset['number'] not in range(195, 231):
                warn('This system is not cubic and will therefore not have a uniform dielectric tensor. However, '
                     'the image-correction schemes that are currently implemented assume a uniform dielectric. '
                     'Proceed with caution')

        else:
            if self.parameters.gb_correction is None:
                self.parameters.gb_correction = False
            if self.parameters.gb_correction:
                raise ValueError('Do not use Gygi-Baldereschi corrections for aperiodic systems')

            if self.parameters.mp_correction is None:
                self.parameters.mp_correction = False
            if self.parameters.mp_correction:
                raise ValueError('Do not use Makov-Payne corrections for aperiodic systems')

            if self.parameters.mt_correction is None:
                self.parameters.mt_correction = True
            if not self.parameters.mt_correction:
                warn('Martyna-Tuckerman corrections not applied for an aperiodic calculation; do this with '
                     'caution')

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs'] and not self.parameters.task.startswith('dft'):
            if len(self.projections) == 0:
                raise ValueError(f'In order to use `init_orbitals={self.parameters.init_orbitals}`, projections must '
                                 'be provided')
            spin_set = set([p.spin for p in self.projections])
            if self.parameters.spin_polarized:
                if spin_set != {Spin.UP, Spin.DOWN}:
                    raise ValueError('This calculation is spin-polarized; please provide spin-up and spin-down '
                                     'projections')
            else:
                if spin_set != {Spin.NONE}:
                    raise ValueError('This calculation is not spin-polarized; please do not provide spin-indexed '
                                     'projections')

        # Check the consistency between self.kpoints.gamma_only and KCP's do_wf_cmplx
        if not self.kpoints.gamma_only and self.calculator_parameters['kcp'].do_wf_cmplx is False:
            warn('`do_wf_cmplx = False` is not consistent with `gamma_only = False`. '
                 'Changing `do_wf_cmplx` to `True`')
            self.calculator_parameters['kcp'].do_wf_cmplx = True

        # Make sanity checks for the ML model
        if self.ml.predict or self.ml.train or self.ml.test:
            warn("Predicting screening parameters with machine-learning is an experimental feature; proceed with "
                 "caution")
            if self.ml_model is None:
                raise ValueError("You have requested to train or predict with a machine-learning model, but no model "
                                 "is attached to this workflow. Either set `ml:train` or `predict` to `True` when "
                                 "initializing the workflow, or directly add a model to the workflow's `ml_model` "
                                 "attribute")
            if self.ml_model.estimator_type != self.ml.estimator:
                warn(f'The estimator type of the loaded ML model (`{self.ml_model.estimator_type}`) does not '
                     f'match the estimator type specified in the Workflow settings (`{self.ml.estimator}`). '
                     'Overriding...')
                self.ml.estimator = self.ml_model.estimator_type
            if self.ml_model.descriptor_type != self.ml.descriptor:
                warn(f'The descriptor type of the loaded ML model (`{self.ml_model.descriptor_type}`) does not '
                     f'match the descriptor type specified in the Workflow settings (`{self.ml.descriptor}`). '
                     'Overriding...')
                self.ml.descriptor = self.ml_model.descriptor_type
            if [self.ml.predict, self.ml.train, self.ml.test].count(True) > 1:
                raise ValueError(
                    'Training, testing, and using the ML model are mutually exclusive; change `ml:predict` '
                    '/`ml:train`/`ml:test` so that at most one is `True`')
            if self.parameters.task not in ['singlepoint', 'trajectory']:
                raise NotImplementedError(
                    f'Using the ML-prediction for the task `{self.parameters.task}` has not yet been implemented.')
            if self.parameters.method != 'dscf':
                raise NotImplementedError(
                    f"Using the ML-prediction for the method `{self.parameters.method}` has not yet been implemented")
            if self.parameters.functional != 'ki':
                raise NotImplementedError(
                    f'Using the ML-prediction for the `{self.parameters.functional}` functional has not yet been '
                    'implemented.')
            if self.parameters.init_orbitals != 'mlwfs':
                raise NotImplementedError(
                    f'Using the ML-prediction for `{self.parameters.init_orbitals}` initial orbitals has not yet been '
                    'implemented.')
            if self.parameters.init_empty_orbitals != self.parameters.init_orbitals:
                raise NotImplementedError(
                    'Using the ML-prediction for using different initial orbitals for empty states than for occupied '
                    'states has not yet been implemented.')
            if self.parameters.spin_polarized:
                warn('Using the ML-prediction for spin-polarised systems has not yet been extensively tested.')
            if not all(self.atoms.pbc):
                warn('Using the ML-prediction for non-periodic systems has not yet been extensively tested.')
            if self.parameters.orbital_groups:
                warn('Using orbital_groups has not yet been extensively tested.')
            if not np.all(self.atoms.cell.angles() == 90.0):
                raise ValueError("The ML-workflow has only been implemented for simulation cells that have 90° angles")

            # check that each n_max, l_max, r_max and r_min are greater or equal to 0 and that r_min is smaller than
            # r_max
            assert isinstance(self.ml.n_max, int)
            assert isinstance(self.ml.l_max, int)
            assert isinstance(self.ml.r_max, float)
            assert isinstance(self.ml.r_min, float)
            if not self.ml.n_max > 0:
                raise ValueError(f"`n_max` has to be larger than zero. The provided value is `n_max = {self.ml.n_max}`")
            if not self.ml.l_max >= 0:
                raise ValueError(
                    f"`l_max` has to be equal or larger than zero. The provided value is `l_max = {self.ml.l_max}`")
            if not self.ml.r_min >= 0:
                raise ValueError(
                    f"`r_min` has to be equal or larger than zero. The provided value is `r_min = {self.ml.r_min}`")
            if self.ml.r_min < 0.5:
                warn(
                    "Small values of `r_min` (<0.5) can lead to problems in the construction of the radial basis. "
                    f"The provided value is `r_min = {self.ml.r_min}`.")
            if not self.ml.r_min < self.ml.r_max:
                raise ValueError(f"`r_min` is larger or equal to `r_max = {self.ml.r_max}`.")

        if self.parameters.initialize_with_smearing:
            if self.calculator_parameters['pw'].occupations != 'smearing':
                self.calculator_parameters['pw'].occupations = 'smearing'
            if self.calculator_parameters['pw'].smearing is None:
                self.calculator_parameters['pw'].smearing = 'cold'
            if self.calculator_parameters['pw'].degauss is None:
                self.calculator_parameters['pw'].degauss = 2e-3

    def new_calculator(self,
                       calc_type: str,
                       kpts: Optional[Union[List[int], BandPath]] = None,
                       **kwargs) -> calculators.Calc:  # type: ignore[type-var, misc]
        """Create a new calculator of a given type."""
        calc_class: calculators.CalcType

        if calc_type == 'kcp':
            calc_class = calculators.KoopmansCPCalculator
        elif calc_type == 'pw':
            calc_class = calculators.PWCalculator
        elif calc_type.startswith('w90'):
            calc_class = calculators.Wannier90Calculator
        elif calc_type == 'pw2wannier':
            calc_class = calculators.PW2WannierCalculator
        elif calc_type == 'wann2kcp':
            calc_class = calculators.Wann2KCPCalculator
        elif calc_type == 'kcw_wannier':
            calc_class = calculators.Wann2KCCalculator
        elif calc_type == 'kcw_screen':
            calc_class = calculators.KoopmansScreenCalculator
        elif calc_type == 'kcw_ham':
            calc_class = calculators.KoopmansHamCalculator
        elif calc_type == 'ph':
            calc_class = calculators.PhCalculator
        elif calc_type == 'projwfc':
            calc_class = calculators.ProjwfcCalculator
        elif calc_type == 'wjl':
            calc_class = calculators.WannierJLCalculator
        else:
            raise ValueError(f'Cound not find a calculator of type `{calc_type}`')

        # Merge calculator_parameters and kwargs, giving kwargs higher precedence
        all_kwargs: Dict[str, Any] = {'parent_process': self}
        calculator_parameters = self.calculator_parameters[calc_type]
        all_kwargs.update(**calculator_parameters)
        all_kwargs.update(**kwargs)

        # For the k-points, the Workflow has two options: self.kpoints.grid and self.kpoints.path. A calculator should
        # only ever have one of these two. By default, use the grid.
        if 'kpts' in calculator_parameters.valid:
            all_kwargs['kpts'] = kpts if kpts is not None else self.kpoints.grid

        # Add further information to the calculator as required
        for kw in ['pseudopotentials', 'gamma_only', 'kgrid', 'kpath', 'koffset', 'plotting']:
            if kw not in all_kwargs and calculator_parameters.is_valid(kw):
                val: Any
                if kw == 'kgrid':
                    val = self.kpoints.grid
                elif kw == 'kpath':
                    val = self.kpoints.path
                elif kw == 'koffset':
                    if (calc_class == calculators.PWCalculator and all_kwargs['calculation'] == 'nscf'
                            or calc_class == calculators.Wannier90Calculator) and self.kpoints.offset_nscf is not None:
                        val = self.kpoints.offset_nscf
                    else:
                        val = self.kpoints.offset
                elif kw == 'gamma_only':
                    val = self.kpoints.gamma_only
                elif kw == 'pseudopotentials':
                    val = OrderedDict(**{k: v.filename.name for k, v in self.pseudopotentials.items()})
                else:
                    val = getattr(self, kw)
                all_kwargs[kw] = val

        # Create the calculator
        calc = calc_class(atoms=copy.deepcopy(self.atoms), **all_kwargs)

        # Link the engine
        calc.engine = self.engine

        # Link the pseudopotentials if relevant
        if calculator_parameters.is_valid('pseudo_dir') or isinstance(calc, calculators.KoopmansHamCalculator):
            for pseudo in self.pseudopotentials.values():
                calc.link_pseudopotential(pseudo.filename, Path('pseudopotentials') / pseudo.filename.name)
            calc.parameters.pseudo_dir = 'pseudopotentials'

        return calc

    def primitive_to_supercell(self, matrix: Optional[npt.NDArray[np.int_]] = None, **kwargs):
        """Convert from a primitive cell to a supercell, as given by a 3x3 transformation matrix."""
        if matrix is None:
            if self.kpoints.gamma_only:
                matrix = np.identity(3, dtype=float)
            else:
                assert self.kpoints.grid is not None
                matrix = np.diag(self.kpoints.grid)
        assert np.shape(matrix) == (3, 3)
        self.atoms = make_supercell(self.atoms, matrix, **kwargs)

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

    def run_steps(self, steps: ProcessProtocol | Sequence[ProcessProtocol], additional_flags: list[str] = []) -> Status:
        """Run a sequence of steps in the workflow.

        Run either a step or sequence (i.e. list) of steps. If a list is provided, the steps must be able to be run
        in parallel.

        :param steps: The step or steps to run
        :type steps: ProcessProtocol | Sequence[ProcessProtocol]
        :return: The status of this set of steps
        :rtype: Status
        """
        assert self.engine is not None

        logger = logging.getLogger(__name__)

        # Convert steps to an immutable, covariant tuple
        if isinstance(steps, ProcessProtocol):
            steps = (steps,)
        if not isinstance(steps, tuple):
            steps = tuple(steps)

        for step in steps:
            step.parent_process = self
            step.engine = self.engine
            self.step_counter += 1

            # Check that calc.directory was not set
            if step.directory_has_been_set():
                raise ValueError(
                    f'`{step.name}.directory` should not be set manually, but it was set to `{step.directory}`')

            assert self.directory is not None
            step.directory = self.directory / f'{self.step_counter:02d}-{step.name}'

            status = self.engine.get_status(step)
            if status == Status.NOT_STARTED:
                logger.info(f'Submitting {step.directory} ({step.__class__.__name__})')
                self.engine.run(step, additional_flags=additional_flags)
            else:
                logger.info(f'Step {step.uid} is {status.name}', stacklevel=2)

        for step in steps:
            if self.engine.get_status(step) == Status.COMPLETED:
                # Load the results of the step
                logger.info(f'Loading the results for {step.directory} ({step.__class__.__name__})')
                self.engine.load_results(step)

                # Store the results of the step
                for step in steps:
                    if step not in self.steps:
                        self.steps.append(step)
                    if isinstance(step, calculators.ImplementedCalc) and step not in self.calculations:
                        self.calculations.append(step)
                        # Ensure we inherit any modifications made to the atoms object
                        if step.atoms != self.atoms:
                            self.atoms = step.atoms
                    elif isinstance(step, Process) and step not in self.processes:
                        self.processes.append(step)

        # If any of the steps are running, return RUNNING
        for step in steps:
            if self.engine.get_status(step) == Status.RUNNING:
                return Status.RUNNING
        # Otherwise, return COMPLETED
        return Status.COMPLETED

    def steps_are_running(self) -> bool:
        """Return True if any steps have status RUNNING."""
        assert self.engine is not None
        for step in self.steps:
            if self.engine.get_status(step) == Status.RUNNING:
                return True
        return False

    @contextmanager
    def _parent_context(self, subdirectory: Optional[str] = None,
                        copy_outputs_to_parent: Optional[bool] = True) -> Generator[None, None, None]:
        """Provide context for running a subworkflow within a parent workflow.

        This is the context for calling self._run(), within which self inherits relevant information from
        self.parent_process, runs, and then passes back relevant information to self.parent_process.
        """
        assert isinstance(self.parent_process, Workflow)

        # Increase the indent level
        self.print_indent = self.parent_process.print_indent + 2

        # Ensure altering self.calculator_parameters won't affect self.parent.calculator_parameters
        if self.calculator_parameters is self.parent_process.calculator_parameters:
            self.calculator_parameters = copy.deepcopy(self.parent_process.calculator_parameters)

        # Copy the lists of calculations, processes, and steps
        self.calculations = [c for c in self.parent_process.calculations]
        self.processes = [p for p in self.parent_process.processes]
        self.steps = [s for s in self.parent_process.steps]

        # Store the current length of the lists so we can tell how many calculations are new
        n_steps = len(self.steps)
        n_calculations = len(self.calculations)
        n_processes = len(self.processes)

        # Link the variational orbitals
        if self.parent_process.variational_orbitals is not None:
            if copy_outputs_to_parent:
                self.variational_orbitals = self.parent_process.variational_orbitals
            else:
                self.variational_orbitals = copy.deepcopy(self.parent_process.variational_orbitals)

        # Prepend the step counter to the subdirectory name
        subdirectory_str = self.name.replace(' ', '-').lower() if subdirectory is None else subdirectory
        self.parent_process.step_counter += 1
        subdirectory_path = Path(f'{self.parent_process.step_counter:02}-{subdirectory_str}')
        assert self.parent_process.directory is not None
        self.directory = self.parent_process.directory / subdirectory_path
        try:
            # Run the workflow
            yield
        finally:
            if copy_outputs_to_parent:
                # Updating the lists of calculations, processes, and steps
                self.parent_process.calculations += self.calculations[n_calculations:]
                self.parent_process.processes += self.processes[n_processes:]
                self.parent_process.steps += self.steps[n_steps:]

                if self.status == Status.COMPLETED:
                    if self.variational_orbitals is not None and self.parent_process.variational_orbitals is None:
                        # Copy the entire variational_orbitals object
                        self.parent_process.variational_orbitals = self.variational_orbitals

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
            engine = LocalhostEngine(**bigdct.pop('engine', {}))
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
        if self.parameters.functional in ['ki', 'kipz', 'pkipz', 'all']:
            add_ref('Dabo2010', 'One of the founding Koopmans functionals papers')
            add_ref('Borghi2014', 'One of the founding Koopmans functionals papers')

            if all(self.atoms.pbc):
                add_ref('Nguyen2018', 'Describes Koopmans functionals in periodic systems')
                if self.parameters.calculate_alpha:
                    if self.parameters.method == 'dfpt':
                        add_ref('Colonna2019', 'Introduces the DFPT method for calculating screening parameters')
                        add_ref('Colonna2022', 'Describes the algorithms underpinning the kcw.x code')
                    else:
                        add_ref('DeGennaro2022', 'Describes how to extract band structures from Koopmans functional '
                                'calculations')
                        add_ref('Borghi2015', 'Describes the algorithms underpinning the kcp.x code')
            else:
                add_ref('Borghi2015', 'Describes the algorithms underpinning the kcp.x code')

            for citation in pseudopotential_library_citations(self.parameters.pseudo_library):
                add_ref(citation,
                        f'Citation for the {self.parameters.pseudo_library.replace("_", " ")} pseudopotential library')

        utils.print_alert(
            'note', f'Please cite the papers listed in `{self.name}.bib` in work involving this calculation')
        relevant_references.to_file(self.name + '.bib')

    def print_preamble(self):
        """Print the preamble of the workflow."""
        if self.parent_process:
            return

        self.print(header())

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
        calcdct['w90'] = {}
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
            if not params_dict and not code.startswith('w90'):
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
            elif code.startswith('w90'):
                if '_' in code:
                    [spin] = [s for s in Spin if s.value == code.split('_')[1]]
                    calcdct['w90'][spin.value] = params_dict
                else:
                    spin = Spin.NONE
                    calcdct['w90'] = params_dict
                projections = self.projections.get_subset(spin) if self.projections else None
                if projections:
                    proj_kwarg = {'projections': [p.projections for p in projections]}
                else:
                    proj_kwarg = {}

                if spin != Spin.NONE:
                    calcdct['w90'][spin.value].update(**proj_kwarg)
                else:
                    calcdct['w90'].update(**proj_kwarg)
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
        defaults = {'colors': colors, 'emin': self.plotting.Emin, 'emax': self.plotting.Emax}
        for b, kwargs in zip(bs, bsplot_kwargs):
            for k, v in defaults.items():
                if k not in kwargs:
                    kwargs[k] = v
            ax_bs = b.plot(ax=ax_bs, **kwargs)

        # Move the legend (if there is one)
        if ax_bs.get_legend():
            ax_bs.legend(loc='lower left', bbox_to_anchor=(0, 1), frameon=False, ncol=min((2, len(bs))))

        if dos is None:
            axes = [ax_bs]
        else:
            # Assemble the densities of state
            dos_summed = dos.sum_by('symbol', 'n', 'l', 'spin')
            if self.parameters.spin_polarized:
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
                    if (not self.parameters.spin_polarized or d.info.get('spin') == 'up') \
                            and all([key in d.info for key in ['symbol', 'n', 'l']]):
                        label = f'{d.info["symbol"]} {d.info["n"]}{d.info["l"]}'
                    else:
                        label = None
                    d.plot_dos(ax=ax_dos, xmin=xmin, xmax=xmax, orientation='vertical', mplargs={'label': label},
                               **dosplot_kwargs)

                # Reset color cycle so the colors of spin-up match those of spin-down
                ax_dos.set_prop_cycle(None)

            # Tweaking the DOS figure aesthetics
            maxval = 1.1 * dos_summed._weights[:, [e >= xmin and e <= xmax for e in dos_summed._energies]].max()
            if self.parameters.spin_polarized:
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
        legends = [ax.get_legend() for ax in axes if ax.get_legend() is not None]
        utils.savefig(fname=f'{self.directory}/{filename}.png', bbox_extra_artists=legends, bbox_inches='tight')

        # Close the figure
        plt.close()

    def _remove_tmpdirs(self):
        """Remove temporary directories created during the workflow."""
        for calc in self.calculations:
            if calc.parameters.get('outdir', None) is None:
                continue
            outdir = File(calc, calc.parameters.outdir)
            if outdir.exists():
                outdir.unlink()

    def _teardown(self):
        """Perform final tasks before the workflow completes."""
        from koopmans.io import write

        assert self.status == Status.COMPLETED
        assert self.engine is not None

        # Save workflow to file
        write(self, self.name + '.pkl')

        # Save the ML model to a separate file
        if self.ml.train:
            assert self.ml_model is not None
            write(self.ml_model, self.name + '_ml_model.pkl')

        # Removing tmpdirs
        if not self.engine.keep_tmpdirs:
            self._remove_tmpdirs()

        # Engine-dependent teardown
        self.engine._teardown()

    def number_of_electrons(self, spin: Spin = Spin.NONE, params: Optional[settings.SettingsDict] = None) -> int:
        return _number_of_electrons(self.inputs.atoms, self.inputs.pseudopotentials, self.config.spin_polarized, spin, params)


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


def header():
    """Print the output header."""
    from koopmans import __version__

    bf = '**' if sys.stdout.isatty() else ''

    header = ["",
              bf + "koopmans" + bf,
              bf + "========" + bf,
              "",
              "*Koopmans spectral functional calculations with `Quantum ESPRESSO`*",
              "",
              f"📦 **Version:** {__version__}  ",
              "🧑 **Authors:** Edward Linscott, Nicola Colonna, Riccardo De Gennaro, Ngoc Linh Nguyen, Giovanni "
              "Borghi, Andrea Ferretti, Ismaila Dabo, and Nicola Marzari  ",
              "📚 **Documentation:** https://koopmans-functionals.org  ",
              "❓ **Support:** https://groups.google.com/g/koopmans-users  ",
              "🐛 **Report a bug:** https://github.com/epfl-theos/koopmans/issues/new"
              ]

    return '\n'.join(header)


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
        ) if isinstance(v[0].parent_process, ParentProcessPlaceholder)}
        calc_nspin1.linked_files = {k: v for k, v in calc_nspin1.linked_files.items(
        ) if isinstance(v[0].parent_process, ParentProcessPlaceholder)}

        # nspin=1 dummy
        calc_nspin1_dummy = master_calc.nspin1_dummy_calculator()
        calc_nspin1_dummy.skip_qc = True
        status = wf.run_steps(calc_nspin1_dummy)
        if status != Status.COMPLETED:
            return status

        calc_nspin1.link(File(calc_nspin1_dummy, calc_nspin1_dummy.parameters.outdir),
                         calc_nspin1.parameters.outdir)

        # Copy over nspin=2 wavefunction to nspin=1 tmp directory
        process2to1 = ConvertFilesFromSpin2To1(name=None,
                                               **prev_calc_nspin2.files_to_convert_with_spin2_to_spin1)
        status = wf.run_steps(process2to1)
        if status != Status.COMPLETED:
            return status

        # Linking converted files to the nspin=1 calculation
        for f in process2to1.outputs.generated_files:
            dst_dir = calc_nspin1.parameters.outdir / \
                f'{calc_nspin1.parameters.prefix}_{calc_nspin1.parameters.ndr}.save/K00001'
            calc_nspin1.link(f, dst_dir / f.name, symlink=True, overwrite=True)

    status = wf.run_steps(calc_nspin1)
    if status != Status.COMPLETED:
        return status

    # nspin=2 from scratch (dummy run for creating files of appropriate size)
    calc_nspin2_dummy = master_calc.nspin2_dummy_calculator()
    calc_nspin2_dummy.skip_qc = True
    status = wf.run_steps(calc_nspin2_dummy)
    if status != Status.COMPLETED:
        return status
    master_calc.link(File(calc_nspin2_dummy, calc_nspin2_dummy.parameters.outdir), master_calc.parameters.outdir)

    # Copy over nspin=1 wavefunction to nspin=2 tmp directory
    process1to2 = ConvertFilesFromSpin1To2(**calc_nspin1.files_to_convert_with_spin1_to_spin2)
    status = wf.run_steps(process1to2)
    if status != Status.COMPLETED:
        return status

    # nspin=2, reading in the spin-symmetric nspin=1 wavefunction
    master_calc.prepare_to_read_nspin1()
    for f in process1to2.outputs.generated_files:
        dst_dir = master_calc.parameters.outdir / \
            f'{master_calc.parameters.prefix}_{master_calc.parameters.ndr}.save/K00001'
        master_calc.link(f, dst_dir / f.name, symlink=True, overwrite=True)
    status = wf.run_steps(master_calc)
    return status


def _number_of_electrons(atoms: Atoms, pseudopotentials: OrderedDict[str, UPFDict], spin_polarized: bool, spin: Spin = Spin.NONE, params: Optional[settings.SettingsDict] = None) -> int:
    """Return the number of electrons in a particular spin channel."""
    nelec_tot = nelec_from_pseudos(atoms, pseudopotentials)
    if params is None:
        params = wf.calculator_parameters['pw']
    nelec = nelec_tot - params.get('tot_charge', 0)
    if spin_polarized:
        if spin == Spin.UP:
            nelec += params.tot_magnetization
            nelec //= 2
        elif spin == Spin.DOWN:
            nelec -= params.tot_magnetization
            nelec //= 2
    else:
        if spin != Spin.NONE:
            nelec //= 2
    return nelec
