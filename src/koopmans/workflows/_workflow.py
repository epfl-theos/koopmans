"""
Generic workflow object for koopmans

Written by Edward Linscott Oct 2020

Converted workflows from functions to objects Nov 2020
"""

from __future__ import annotations

import copy
import json as json_ext
import os
import re
import shutil
import sys
from abc import ABC, abstractmethod
from collections import OrderedDict
from contextlib import contextmanager
from pathlib import Path
from typing import (Any, Callable, Dict, Generator, List, Optional, Sequence,
                    Tuple, Type, TypeVar, Union)

import dill
import numpy as np
from numpy import typing as npt
from pybtex.database import BibliographyData

# isort: off
import koopmans.mpl_config
import matplotlib.pyplot as plt
# isort: on

from ase_koopmans import Atoms
from ase_koopmans.build.supercells import make_supercell
from ase_koopmans.calculators.calculator import CalculationFailed
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

from koopmans import calculators, outputs, settings, utils
from koopmans.bands import Bands
from koopmans.commands import ParallelCommandWithPostfix
from koopmans.engines import Engine, LocalhostEngine
from koopmans.files import File, LocalFile
from koopmans.kpoints import Kpoints
from koopmans.ml import AbstractMLModel, MLModel, OccEmpMLModels
from koopmans.processes import Process
from koopmans.processes.koopmans_cp import (ConvertFilesFromSpin1To2,
                                            ConvertFilesFromSpin2To1)
from koopmans.projections import ProjectionBlocks
from koopmans.pseudopotentials import (nelec_from_pseudos,
                                       pseudopotential_library_citations)
from koopmans.references import bib_data
from koopmans.status import Status
from koopmans.step import Step
from koopmans.utils import SpinType

T = TypeVar('T', bound='calculators.CalculatorExt')
W = TypeVar('W', bound='Workflow')


class Workflow(utils.HasDirectory, ABC):

    r'''
    Abstract base class that defines a Koopmans workflow

    Parameters
    ----------

    atoms : Atoms
        an ASE ``Atoms`` object defining the atomic positions, cell, etc
    pseudopotentials : OrderedDict[str, str]
        a dictionary mapping atom labels to pseudopotential filenames
    kpoints : koopmans.kpoints.Kpoints
        a dataclass defining the k-point sampling and paths
    projections : ProjectionsBlocks
        The projections to be used in the Wannierization
    name : str
        a name for the workflow
    parameters : Dict[str, Any] | koopmans.settings.WorkflowSettingsDict
        a dictionary specifying any workflow settings to use; note that a simpler alternative is to provide workflow
        settings as keyword arguments
    calculator_parameters : Dict[str, koopmans.settings.SettingsDict]
        a dictionary containing calculator-specific settings; as for the parameters, it is usually simpler to specify
        these individually as keyword arguments
    plotting : koopmans.settings.PlotSettingsDict
        a dictionary containing settings specific to plotting; again, it is usually simpler to specify these
        individually as keyword arguments
    autogenerate_settings : bool
        if True (the default), autogenerate various calculator settings; the only scenario where you do not want to do
        this is when creating a new workflow from a .kwf file
    **kwargs
        any valid workflow, calculator, or plotting settings e.g. ``{"functional": "ki", "ecutwfc": 50.0}``
    '''

    atoms: Atoms
    parameters: settings.WorkflowSettingsDict
    calculator_parameters: Dict[str, settings.SettingsDict]
    name: str
    kpoints: Kpoints
    pseudopotentials: Dict[str, UPFDict]
    pseudo_dir: Path
    projections: ProjectionBlocks
    ml_model: Optional[AbstractMLModel]
    snapshots: List[Atoms]
    version: str
    step_counter: int
    print_indent: int
    status: Status
    engine: Engine

    __slots__ = utils.HasDirectory.__slots__ + ['atoms', 'parameters', 'calculator_parameters', 'name', 'kpoints',
                                                'pseudopotentials', 'projections', 'ml_model', 'snapshots',
                                                'version', 'calculations', 'processes',
                                                'steps', 'plotting', 'ml', '_bands', 'step_counter', 'print_indent',
                                                'status']

    def __init__(self,
                 atoms: Atoms,
                 engine: Engine,
                 pseudopotentials: Dict[str, UPFDict | str] = {},
                 kpoints: Optional[Kpoints] = None,
                 projections: Optional[ProjectionBlocks] = None,
                 name: Optional[str] = None,
                 parameters: Union[Dict[str, Any], settings.WorkflowSettingsDict] = {},
                 calculator_parameters: Optional[Union[Dict[str, Dict[str, Any]],
                                                       Dict[str, settings.SettingsDict]]] = None,
                 plotting: Union[Dict[str, Any], settings.PlotSettingsDict] = {},
                 ml: Union[Dict[str, Any], settings.MLSettingsDict] = {},
                 ml_model: Optional[MLModel] = None,
                 snapshots: Optional[List[Atoms]] = None,
                 autogenerate_settings: bool = True,
                 version: Optional[str] = None,
                 parent: Optional[Workflow] = None,
                 **kwargs: Dict[str, Any]):

        self.step_counter = 0
        self.print_indent = 0
        self.status = Status.NOT_STARTED

        # Initialize the HasDirectory information (parent, base_directory, directory)
        base_directory = None if parent else Path.cwd()
        directory = None if parent else Path()
        super().__init__(parent=parent, base_directory=base_directory, directory=directory, engine=engine)
        assert self.engine == engine

        # Parsing parameters
        self.parameters = settings.WorkflowSettingsDict(**parameters)
        for key, value in kwargs.items():
            if self.parameters.is_valid(key):
                self.parameters[key] = value
        self.atoms: Atoms = atoms
        self.snapshots = snapshots if snapshots is not None else [self.atoms]
        if name:
            self.name = name
        else:
            name_with_split_acroynms = re.sub(r'([a-z])([A-Z])', r'\1 \2',
                                              self.__class__.__name__.replace('Workflow', ''))
            self.name = re.sub(r'([A-Z])([A-Z][a-z])', r'\1 \2', name_with_split_acroynms)
        self.calculations: List[calculators.Calc] = []
        self.processes: List[Process] = []
        self.steps: List = []
        self._bands: Optional[Bands] = None

        if projections is None:
            proj_list: List[List[Any]]
            spins: List[SpinType]
            if self.parameters.spin_polarized:
                proj_list = [[], []]
                spins = ['up', 'down']
            else:
                proj_list = [[]]
                spins = [None]
            self.projections = ProjectionBlocks.fromlist(proj_list, spins=spins, atoms=self.atoms)
        else:
            self.projections = projections

        self.plotting = settings.PlotSettingsDict(**plotting)
        for key, value in kwargs.items():
            if self.plotting.is_valid(key):
                self.plotting[key] = value

        self.ml = settings.MLSettingsDict(**ml)  # , task=self.parameters.task)
        for key, value in kwargs.items():
            if self.ml.is_valid(key):
                self.ml[key] = value

        # Initialize the MLModel
        if ml_model is not None:
            self.ml_model = ml_model
        elif self.ml.train:
            assert isinstance(self.ml.estimator, str)
            assert isinstance(self.ml.descriptor, str)
            if self.ml.occ_and_emp_together:
                self.ml_model = MLModel(self.ml.estimator, self.ml.descriptor, engine=self.engine)
            else:
                self.ml_model = OccEmpMLModels(self.ml.estimator, self.ml.descriptor, engine=self.engine)
        elif self.ml.predict or self.ml.test:
            if self.ml.model_file is None:
                raise ValueError('Cannot initialize a Workflow with `ml.predict = True` without providing a model via '
                                 'the `ml.model_file` setting or the `ml_model` argument')
            with open(self.ml.model_file, 'rb') as f:
                self.ml_model = dill.load(f)
            assert isinstance(self.ml_model, AbstractMLModel)

        else:
            self.ml_model = None

        if all(self.atoms.pbc):
            self.atoms.wrap(pbc=True)

        # kpoints
        if all(self.atoms.pbc):
            # By default, use ASE's default bandpath for this cell (see
            # https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#brillouin-zone-data)
            default_path = self.atoms.cell.bandpath().path
        else:
            default_path = 'G'
        if kpoints is None:
            kpoints = Kpoints(path=default_path, gamma_only=False, cell=self.atoms.cell)
        elif kpoints.path is None:
            kpoints.set_path(default_path, cell=self.atoms.cell)
        self.kpoints = kpoints

        # Pseudopotentials and pseudo_dir
        if self.parameters.pseudo_library is None:
            raise NotImplementedError()

        if self.parameters.pseudo_library is None:
            raise ValueError('No pseudopotential library was provided`')
        elif pseudopotentials:
            # Ensure pseudopotentials are converted to UPFDict objects
            self.pseudopotentials = {k: UPFDict.from_upf(
                self.parameters.pseudo_directory / v) if isinstance(v, str) else v for k, v in pseudopotentials.items()}
        else:
            self.pseudopotentials = {}
            for element, tag in set([(a.symbol, a.tag) for a in self.atoms]):
                pseudo = self.engine.get_pseudopotential(library=self.parameters.pseudo_library, element=element)

                pseudo_type = pseudo['header']['pseudo_type']
                if pseudo_type != 'NC':
                    raise ValueError('Koopmans functionals only currently support norm-conserving pseudopotentials; '
                                     f'{pseudo.filename} is {pseudo_type}')

                symbol = element + str(tag) if tag > 0 else element
                self.pseudopotentials[symbol] = pseudo

        # Make sure calculator_parameters isn't missing any entries, and every entry corresponds to
        # settings.SettingsDict objects
        calculator_parameters = sanitize_calculator_parameters(calculator_parameters) if calculator_parameters \
            is not None else generate_default_calculator_parameters()

        # For any kwargs...
        for key, value in kwargs.items():
            match = False
            # if they correspond to any valid calculator parameter, set it
            for calc_params in calculator_parameters.values():
                if calc_params.is_valid(key):
                    calc_params[key] = value
                    match = True
            # if not a calculator, workflow, or plotting keyword, raise an error
            if not match and not self.parameters.is_valid(key) and not self.plotting.is_valid(key) \
                    and not self.ml.is_valid(key):
                raise ValueError(f'`{key}` is not a valid setting')

        # Before saving the calculator_parameters, automatically generate some keywords and perform some sanity checks
        if self.parameters.task != 'ui' and autogenerate_settings:
            # Automatically calculate nelec/nelup/neldw/etc using information contained in the pseudopotential files
            # and the kcp settings
            nelec = nelec_from_pseudos(self.atoms, self.pseudopotentials)
            tot_charge = calculator_parameters['pw'].get(
                'tot_charge', calculator_parameters['kcp'].get('tot_charge', 0))
            nelec -= tot_charge
            tot_mag = calculator_parameters['pw'].get(
                'tot_magnetization', calculator_parameters['kcp'].get('tot_charge', nelec % 2))
            nelup = int(nelec / 2 + tot_mag / 2)
            neldw = int(nelec / 2 - tot_mag / 2)

            # Setting up the magnetic moments
            starting_mags = {k: v for k, v in calculator_parameters['kcp'].items(
            ) if k.startswith('starting_magnetization')}
            starting_mags.update(
                {k: v for k, v in calculator_parameters['pw'].items() if k.startswith('starting_magnetization')})
            if starting_mags:
                labels = [a.symbol + str(a.tag) if a.tag > 0 else a.symbol for a in self.atoms]
                starting_magmoms: Dict[str, float] = {}
                for i, (l, p) in enumerate(self.pseudopotentials.items()):
                    mag = starting_mags.pop(f'starting_magnetization({i + 1})', 0.0)
                    if abs(mag) < 1.0:
                        # If |mag| < 1, QE interprets this as site magnetization *per valence electron*, whereas ASE
                        # expects simply the site magnetization
                        valence = p['header']['z_valence']
                        starting_magmoms[l] = mag * valence
                    else:
                        # If |mag| >= 1, QE interprets this as site magnetization
                        starting_magmoms[l] = mag
                atoms.set_initial_magnetic_moments([starting_magmoms[l] for l in labels])
            elif tot_mag != 0:
                atoms.set_initial_magnetic_moments([tot_mag / len(atoms) for _ in atoms])

            # Work out the number of bands
            nbnd = calculator_parameters['kcp'].get('nbnd', nelec // 2 + nelec % 2)
            generated_keywords = {'nelec': nelec, 'tot_charge': tot_charge, 'tot_magnetization': tot_mag,
                                  'nelup': nelup, 'neldw': neldw, 'nbnd': nbnd}
        else:
            generated_keywords = {}
            nelec = 0

        self.calculator_parameters = generate_default_calculator_parameters()
        for block, params in calculator_parameters.items():
            # Apply auto-generated keywords
            for k, v in generated_keywords.items():
                # Skipping nbnd for kcp -- it is valid according to ASE but it is not yet properly implemented
                if k == 'nbnd' and block == 'kcp':
                    continue
                if k in params.valid and k not in params:
                    setattr(params, k, v)

            # Various checks for the wannier90 blocks
            if block.startswith('w90'):
                # If we are spin-polarized, don't store the spin-independent w90 block
                # Likewise, if we are not spin-polarized, don't store the spin-dependent w90 blocks
                if self.parameters.spin_polarized is not ('up' in block or 'down' in block):
                    continue
                if 'projections' in params or 'projections_blocks' in params:
                    raise ValueError(f'You have provided projection information in the '
                                     f'`calculator_parameters[{block}]` argument to `{self.__class__.__name__}`. '
                                     'Please instead specify projections via the `projections` argument')
                for kw in ['exclude_bands', 'num_wann', 'num_bands', 'projections']:
                    if kw in params:
                        utils.warn(f'`{kw}` will be overwritten by the workflow; it is best to not specify this '
                                   'keyword and to instead double-check the keyword in the various `.win` files '
                                   'generated by the workflow.')

            # Replace "nelec" as a unit with its numerical value, for any settings that are specified implicitly in
            # terms of nelec
            params.replace_nelec(nelec)

            # Store the sanitized parameters
            self.calculator_parameters[block] = params

        # If atoms has a calculator, overwrite the kpoints and pseudopotentials variables and then detach the
        # calculator
        if atoms.calc is not None:
            utils.warn(f'You have initialized a `{self.__class__.__name__}` object with an atoms object that '
                       'possesses a calculator. This calculator will be ignored.')
            self.atoms.calc = None

        # Adding excluded_bands info to self.projections
        if self.projections:
            for spin in ['up', 'down', None]:
                label = 'w90'
                if spin:
                    label += f'_{spin}'
                self.projections.exclude_bands[spin] = self.calculator_parameters[label].get('exclude_bands', [])

        # Version number (important when loading workflows from .kwf files)
        from koopmans import __version__
        self.version = version if version is not None else __version__
        if self.version != __version__:
            utils.warn(f'You are using version {__version__} of `koopmans`, but this workflow was generated with '
                       f'version {self.version}. Proceed with caution.')

    def __eq__(self, other: Any):
        if self.__class__ != other.__class__:
            return False
        else:
            for key in self.__slots__:
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

    def run_while(self,):
        '''
        Run the workflow
        '''
        self.print_preamble()
        attempts = 0

        if not self.parent:
            bf = '**' if sys.stdout.isatty() else ''
            self.print(bf + self.name + bf)
            self.print(bf + '-' * len(self.name) + bf)
            if self.parameters.from_scratch:
                self._remove_tmpdirs()
            self._run_sanity_checks()

        while self.status != Status.COMPLETED:
            # Reset the step counter each time we reattempt the workflow
            self.run()

            attempts += 1
            if attempts == 1000:
                self.status = Status.FAILED
                raise ValueError('Workflow failed to complete after 1000 attempts')

            self.engine.update_statuses()

        if not self.parent and self.status == Status.COMPLETED:
            self._teardown()

        self.status = Status.COMPLETED
        return

    def run(self, subdirectory: Optional[str] = None, copy_outputs_to_parent: Optional[bool] = True):
        '''
        Run the workflow
        '''

        self.step_counter = 0

        if self.status == Status.NOT_STARTED:
            self.status = Status.RUNNING

        if self.parent:
            with self._parent_context(subdirectory, copy_outputs_to_parent):
                self._run()
        else:
            if subdirectory:
                self.base_directory = Path(subdirectory).resolve()
                with utils.chdir(subdirectory):
                    self._run()
            else:
                self.base_directory = Path.cwd().resolve()
                self._run()

        self.engine.update_statuses()

    @abstractmethod
    def _run(self) -> None:
        ...

    @property
    @abstractmethod
    def output_model(self) -> Type[outputs.OutputModel]:
        ...

    @classmethod
    def from_other(cls: Type[W], other_wf: Workflow, **kwargs: Any) -> W:
        '''
        Creates a new workflow from another workflow, copying all settings and parameters
        e.g.
        >>> new_wf = Workflow.fromother(other_wf)
        '''

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
        '''
        Creates a subworkflow with the same configuration as the parent workflow
        e.g.
        >>> sub_wf = Workflow.fromparent(self)
        '''
        return cls.from_other(other_wf=parent_wf, parent=parent_wf, **kwargs)

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
                utils.warn('You have requested a ΔSCF calculation with frozen orbitals. This is unusual; proceed '
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
                    utils.warn('Gygi-Baldereschi corrections are not being used; do this with '
                               'caution for periodic systems')

            elif self.parameters.method == 'dscf':
                # For DSCF, we use mp_correction
                if self.parameters.mp_correction is None:
                    self.parameters.mp_correction = True
                if not self.parameters.mp_correction:
                    utils.warn('Makov-Payne corrections are not being used; do this with '
                               'caution for periodic systems')

                if self.parameters.eps_inf is None:
                    utils.warn('`eps_inf` missing in input; it will default to 1.0. Proceed with caution for periodic '
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
                utils.warn('This system is not cubic and will therefore not have a uniform dielectric tensor. However, '
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
                utils.warn('Martyna-Tuckerman corrections not applied for an aperiodic calculation; do this with '
                           'caution')

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs'] and not self.parameters.task.startswith('dft'):
            if len(self.projections) == 0:
                raise ValueError(f'In order to use `init_orbitals={self.parameters.init_orbitals}`, projections must '
                                 'be provided')
            spin_set = set([p.spin for p in self.projections])
            if self.parameters.spin_polarized:
                if spin_set != {'up', 'down'}:
                    raise ValueError('This calculation is spin-polarized; please provide spin-up and spin-down '
                                     'projections')
            else:
                if spin_set != {None}:
                    raise ValueError('This calculation is not spin-polarized; please do not provide spin-indexed '
                                     'projections')

        # Check the consistency between self.kpoints.gamma_only and KCP's do_wf_cmplx
        if not self.kpoints.gamma_only and self.calculator_parameters['kcp'].do_wf_cmplx is False:
            utils.warn('`do_wf_cmplx = False` is not consistent with `gamma_only = False`. '
                       'Changing `do_wf_cmplx` to `True`')
            self.calculator_parameters['kcp'].do_wf_cmplx = True

        # Make sanity checks for the ML model
        if self.ml.predict or self.ml.train or self.ml.test:
            utils.warn("Predicting screening parameters with machine-learning is an experimental feature; proceed with "
                       "caution")
            if self.ml_model is None:
                raise ValueError("You have requested to train or predict with a machine-learning model, but no model "
                                 "is attached to this workflow. Either set `ml:train` or `predict` to `True` when "
                                 "initializing the workflow, or directly add a model to the workflow's `ml_model` "
                                 "attribute")
            if self.ml_model.estimator_type != self.ml.estimator:
                utils.warn(f'The estimator type of the loaded ML model (`{self.ml_model.estimator_type}`) does not '
                           f'match the estimator type specified in the Workflow settings (`{self.ml.estimator}`). '
                           'Overriding...')
                self.ml.estimator = self.ml_model.estimator_type
            if self.ml_model.descriptor_type != self.ml.descriptor:
                utils.warn(f'The descriptor type of the loaded ML model (`{self.ml_model.descriptor_type}`) does not '
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
                utils.warn('Using the ML-prediction for spin-polarised systems has not yet been extensively tested.')
            if not all(self.atoms.pbc):
                utils.warn('Using the ML-prediction for non-periodic systems has not yet been extensively tested.')
            if self.parameters.orbital_groups:
                utils.warn('Using orbital_groups has not yet been extensively tested.')
            if not np.all(self.atoms.cell.angles() == 90.0):
                raise ValueError(f"The ML-workflow has only been implemented for simulation cells that have 90° angles")

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
                utils.warn(
                    "Small values of `r_min` (<0.5) can lead to problems in the construction of the radial basis. "
                    f"The provided value is `r_min = {self.ml.r_min}`.")
            if not self.ml.r_min < self.ml.r_max:
                raise ValueError(f"`r_min` is larger or equal to `r_max = {self.ml.r_max}`.")

    def new_calculator(self,
                       calc_type: str,
                       kpts: Optional[Union[List[int], BandPath]] = None,
                       **kwargs) -> calculators.Calc:  # type: ignore[type-var, misc]

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
        else:
            raise ValueError(f'Cound not find a calculator of type `{calc_type}`')

        # Merge calculator_parameters and kwargs, giving kwargs higher precedence
        all_kwargs: Dict[str, Any] = {'parent': self}
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
                    if (calc_class == calculators.PWCalculator and all_kwargs['calculation'] == 'nscf' or
                            calc_class == calculators.Wannier90Calculator) and self.kpoints.offset_nscf is not None:
                        val = self.kpoints.offset_nscf
                    else:
                        val = self.kpoints.offset
                elif kw == 'gamma_only':
                    val = self.kpoints.gamma_only
                elif kw == 'pseudopotentials':
                    val = {k: v.filename.name for k, v in self.pseudopotentials.items()}
                else:
                    val = getattr(self, kw)
                all_kwargs[kw] = val

        # Create the calculator
        calc = calc_class(atoms=copy.deepcopy(self.atoms), **all_kwargs)

        # Link the pseudopotentials if relevant
        if calculator_parameters.is_valid('pseudo_dir'):
            for pseudo in self.pseudopotentials.values():
                calc.link(LocalFile(pseudo.filename), Path('pseudopotentials') / pseudo.filename.name, symlink=True)
            calc.parameters.pseudo_dir = 'pseudopotentials'

        # Link the engine
        calc.engine = self.engine

        return calc

    def primitive_to_supercell(self, matrix: Optional[npt.NDArray[np.int_]] = None, **kwargs):
        # Converts to a supercell as given by a 3x3 transformation matrix
        if matrix is None:
            if self.kpoints.gamma_only:
                matrix = np.identity(3, dtype=float)
            else:
                assert self.kpoints.grid is not None
                matrix = np.diag(self.kpoints.grid)
        assert np.shape(matrix) == (3, 3)
        self.atoms = make_supercell(self.atoms, matrix, **kwargs)

    def supercell_to_primitive(self, matrix: Optional[npt.NDArray[np.int_]] = None):
        # Converts from a supercell to a primitive cell, as given by a 3x3 transformation matrix
        # The inverse of self.primitive_to_supercell()
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

    def run_steps(self, steps: Step | Sequence[Step]) -> Status:
        """Run a sequence of steps in the workflow

        Run either a step or sequence (i.e. list) of steps. If a list is provided, the steps must be able to be run
        in parallel.

        :param steps: The step or steps to run
        :type steps: Step | Sequence[Step]
        :return: The status of this set of steps
        :rtype: Status
        """
        # Convert steps to an immutable, covariant tuple
        if isinstance(steps, Step):
            steps = (steps,)
        if not isinstance(steps, tuple):
            steps = tuple(steps)

        for step in steps:
            step.parent = self
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
                # Run the step
                self.engine.run(step)

        for step in steps:
            if self.engine.get_status(step) == Status.COMPLETED:
                # Load the results of the step
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
        for step in self.steps:
            if self.engine.get_status(step) == Status.RUNNING:
                return True
        return False

    @contextmanager
    def _parent_context(self, subdirectory: Optional[str] = None,
                        copy_outputs_to_parent: Optional[bool] = True) -> Generator[None, None, None]:
        '''
        Context for calling self._run(), within which self inherits relevant information from self.parent, runs, and
        then passes back relevant information to self.parent
        '''

        assert isinstance(self.parent, Workflow)

        # Increase the indent level
        self.print_indent = self.parent.print_indent + 2

        # Ensure altering self.calculator_parameters won't affect self.parent.calculator_parameters
        if self.calculator_parameters is self.parent.calculator_parameters:
            self.calculator_parameters = copy.deepcopy(self.parent.calculator_parameters)

        # Copy the lists of calculations, processes, and steps
        self.calculations = [c for c in self.parent.calculations]
        self.processes = [p for p in self.parent.processes]
        self.steps = [s for s in self.parent.steps]

        # Store the current length of the lists so we can tell how many calculations are new
        n_steps = len(self.steps)
        n_calculations = len(self.calculations)
        n_processes = len(self.processes)

        # Link the bands
        if self.parent.bands is not None:
            if copy_outputs_to_parent:
                self.bands = self.parent.bands
            else:
                self.bands = copy.deepcopy(self.parent.bands)

        # Prepend the step counter to the subdirectory name
        subdirectory_str = self.name.replace(' ', '-').lower() if subdirectory is None else subdirectory
        self.parent.step_counter += 1
        subdirectory_path = Path(f'{self.parent.step_counter:02}-{subdirectory_str}')
        assert self.parent.directory is not None
        self.directory = self.parent.directory / subdirectory_path
        try:
            # Run the workflow
            yield
        finally:
            if copy_outputs_to_parent:
                # Updating the lists of calculations, processes, and steps
                self.parent.calculations += self.calculations[n_calculations:]
                self.parent.processes += self.processes[n_processes:]
                self.parent.steps += self.steps[n_steps:]

                if self.status == Status.COMPLETED:
                    if self.bands is not None and self.parent.bands is None:
                        # Copy the entire bands object
                        self.parent.bands = self.bands

    def print(self, *args, **kwargs):
        utils.indented_print(*args, **kwargs)

    def todict(self):
        # Shallow copy
        dct = dict(self.__dict__)

        # Adding information required by the json decoder
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @classmethod
    def fromdict(cls, dct: Dict[str, Any], **kwargs) -> Workflow:

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
    def bands(self) -> Bands | None:
        return self._bands

    @bands.setter
    def bands(self, value: Bands):
        assert isinstance(value, Bands)
        self._bands = value

    @classmethod
    def fromjson(cls, fname: str, override: Dict[str, Any] = {}, **kwargs):

        with open(fname, 'r') as fd:
            bigdct = json_ext.loads(fd.read())
        wf = cls._fromjsondct(bigdct, override, **kwargs)

        # Define the name of the workflow using the name of the json file
        wf.name = fname.replace('.json', '')

        return wf

    @classmethod
    def _fromjsondct(cls, bigdct: Dict[str, Any], override: Dict[str, Any] = {}, **kwargs):
        # Override all keywords provided explicitly
        utils.update_nested_dict(bigdct, override)

        # Loading atoms object
        atoms_dict = bigdct.pop('atoms', None)
        if atoms_dict:
            atoms, snapshots = read_atoms_dict(utils.parse_dict(atoms_dict))
        else:
            raise ValueError('Please provide an `atoms` block in the json input file')

        # Loading plot settings
        kwargs['plotting'] = settings.PlotSettingsDict(**utils.parse_dict(bigdct.pop('plotting', {})))

        # Loading workflow settings
        parameters = settings.WorkflowSettingsDict(**utils.parse_dict(bigdct.pop('workflow', {})))

        # Loading ml settings
        kwargs['ml'] = settings.MLSettingsDict(**utils.parse_dict(bigdct.pop('ml', {})))  # , task=parameters['task'])

        # Loading kpoints
        kpts = Kpoints(**utils.parse_dict(bigdct.pop('kpoints', {})), cell=atoms.cell)

        # Loading calculator-specific settings
        calcdict = bigdct.pop('calculator_parameters', {})

        # First, extract the w90 subdictionaries
        if 'w90' in calcdict:
            universal_settings = {k: v for k, v in calcdict['w90'].items() if k not in ['up', 'down']}
            for spin in ['up', 'down']:
                # Add any keywords in the spin subdictionary
                calcdict[f'w90_{spin}'] = calcdict['w90'].pop(spin, {})
                # Add any keywords in the main dictionary
                calcdict[f'w90_{spin}'].update(universal_settings)

        # Secondly, flatten the UI subdictionaries
        if 'ui' in calcdict:
            subdcts = {}
            keys = ['occ', 'emp']
            for key in keys:
                # First, we must remove the occ and emp subdicts from the UI dict
                if key in calcdict['ui']:
                    subdcts[key] = calcdict['ui'].pop(key)

            # Now, we add the ui_occ and ui_emp calculators to calculator_parameters
            for key in keys:
                if key in subdcts:
                    # Add the corresponding subdict to the rest of the UI block
                    calcdict[f'ui_{key}'] = dict(calcdict['ui'], **subdcts[key])
                else:
                    # Duplicate the UI block
                    calcdict[f'ui_{key}'] = calcdict['ui']

        # Third, flatten the kcw subdicts
        kcw_blocks = calcdict.pop('kcw', {'ham': {}, 'screen': {}, 'wannier': {}})
        calcdict.update(**{'kcw_' + k: v for k, v in kcw_blocks.items()})

        # Finally, generate a SettingsDict for every single kind of calculator, regardless of whether or not there was
        # a corresponding block in the json file
        calculator_parameters = {}
        w90_block_projs: List = []
        w90_block_spins: List[SpinType] = []
        for block, settings_class in settings_classes.items():
            # Read the block and add the resulting calculator to the calcs_dct
            dct = calcdict.pop(block, {})
            if block.startswith('ui'):
                # Dealing with redundancies in UI keywords
                if 'sc_dim' in dct and kpts.grid is not None:
                    # In this case, the sc_dim keyword is redundant
                    if kpts.grid != dct['sc_dim']:
                        raise ValueError(
                            '`sc_dim` in the `ui` block should match the kpoints provided in the `kpoints` block')
                    dct.pop('sc_dim')
                if 'kpath' in dct and kpts.path is not None:
                    if kpts.path != dct['kpath']:
                        raise ValueError('`kpath` in the `ui` block should match that provided in the `kpoints` block')
                    dct.pop('kpath')
            elif block.startswith('w90'):
                # If we are spin-polarized, don't store the spin-independent w90 block
                # Likewise, if we are not spin-polarized, don't store the spin-dependent w90 blocks
                if parameters.spin_polarized is not ('up' in block or 'down' in block):
                    continue
                projs = dct.pop('projections', [[]])
                w90_block_projs += projs
                if 'up' in block:
                    w90_block_spins += ['up' for _ in range(len(projs))]
                elif 'down' in block:
                    w90_block_spins += ['down' for _ in range(len(projs))]
                else:
                    w90_block_spins += [None for _ in range(len(projs))]

            calculator_parameters[block] = settings_class(**dct)

        # Adding the projections to the workflow kwargs (this is unusual in that this is an attribute of the workflow
        # object but it is provided in the w90 subdictionary)
        kwargs['projections'] = ProjectionBlocks.fromlist(w90_block_projs, w90_block_spins, atoms)

        kwargs['pseudopotentials'] = bigdct.pop('pseudopotentials', {})

        # Convergence
        if parameters.converge:
            conv_block = settings.ConvergenceSettingsDict(**bigdct.pop('convergence', {}))

        # Check for unexpected blocks
        for block in bigdct:
            raise ValueError(f'Unrecognized block `{block}` in the json input file')

        # Create the workflow. Note that any keywords provided in the calculator_parameters (i.e. whatever is left in
        # calcdict) are provided as kwargs
        wf = cls(atoms, snapshots=snapshots, parameters=parameters, kpoints=kpts,
                 calculator_parameters=calculator_parameters, **kwargs, **calcdict)

        if parameters.converge:
            # Wrap the workflow in a convergence outer loop
            from koopmans.workflows import ConvergenceWorkflowFactory
            return ConvergenceWorkflowFactory(wf, **conv_block)
        else:
            return wf

    def print_bib(self):
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
                add_ref(
                    citation, f'Citation for the {self.parameters.pseudo_library.replace("_", " ")} pseudopotential library')

        utils.print_alert(
            'note', f'Please cite the papers listed in `{self.name}.bib` in work involving this calculation')
        relevant_references.to_file(self.name + '.bib')

    def print_preamble(self):
        if self.parent:
            return

        self.print(header())

        self.print_bib()

    def toinputjson(self) -> Dict[str, Dict[str, Any]]:

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
                    spin = code.split('_')[1]
                    calcdct['w90'][spin] = params_dict
                else:
                    spin = None
                    calcdct['w90'] = params_dict
                projections = self.projections.get_subset(spin)
                if projections:
                    proj_kwarg = {'projections': [p.projections for p in projections]}
                else:
                    proj_kwarg = {}

                if spin:
                    calcdct['w90'][spin].update(**proj_kwarg)
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
        """
        Plots the provided band structure (and optionally also a provided DOS)
        Arguments:
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
        '''
        Removes tmpdirs
        '''

        all_outdirs = [calc.parameters.get('outdir', None) for calc in self.calculations]
        outdirs = set([o.resolve() for o in all_outdirs if o is not None and o.resolve().exists()])
        for outdir in outdirs:
            shutil.rmtree(outdir)

    def _teardown(self):
        '''
        Performs final tasks before the workflow completes
        '''

        assert self.status == Status.COMPLETED

        # Removing tmpdirs
        if not self.parameters.keep_tmpdirs:
            self._remove_tmpdirs()

    def number_of_electrons(self, spin: Optional[str] = None) -> int:
        # Return the number of electrons in a particular spin channel
        nelec_tot = nelec_from_pseudos(self.atoms, self.pseudopotentials)
        pw_params = self.calculator_parameters['pw']
        if self.parameters.spin_polarized:
            nelec = nelec_tot - pw_params.get('tot_charge', 0)
            if spin == 'up':
                nelec += pw_params.tot_magnetization
            else:
                nelec -= pw_params.tot_magnetization
            nelec = int(nelec // 2)
        else:
            nelec = nelec_tot
        return nelec


def read_atoms_dict(dct: Dict[str, Any]) -> Tuple[Atoms, List[Atoms]]:
    '''
    Reads the "atoms" block
    '''

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


def generate_default_calculator_parameters() -> Dict[str, settings.SettingsDict]:
    # Dictionary to be used as the default value for 'calculator_parameters' when initializing a workflow
    # We create this dynamically in order for the .directory attributes to make sense
    return {'kcp': settings.KoopmansCPSettingsDict(),
            'kcw_ham': settings.KoopmansHamSettingsDict(),
            'kcw_screen': settings.KoopmansScreenSettingsDict(),
            'ph': settings.PhSettingsDict(),
            'projwfc': settings.ProjwfcSettingsDict(),
            'pw': settings.PWSettingsDict(),
            'pw2wannier': settings.PW2WannierSettingsDict(),
            'wann2kcp': settings.Wann2KCPSettingsDict(),
            'ui': settings.UnfoldAndInterpolateSettingsDict(),
            'ui_occ': settings.UnfoldAndInterpolateSettingsDict(),
            'ui_emp': settings.UnfoldAndInterpolateSettingsDict(),
            'kcw_wannier': settings.Wann2KCSettingsDict(),
            'w90': settings.Wannier90SettingsDict(),
            'w90_up': settings.Wannier90SettingsDict(),
            'w90_down': settings.Wannier90SettingsDict(),
            'plot': settings.PlotSettingsDict()}


# Define which function to use to read each block
settings_classes = {'kcp': settings.KoopmansCPSettingsDict,
                    'kcw_ham': settings.KoopmansHamSettingsDict,
                    'kcw_screen': settings.KoopmansScreenSettingsDict,
                    'kcw_wannier': settings.Wann2KCSettingsDict,
                    'ph': settings.PhSettingsDict,
                    'projwfc': settings.ProjwfcSettingsDict,
                    'pw': settings.PWSettingsDict,
                    'pw2wannier': settings.PW2WannierSettingsDict,
                    'wann2kcp': settings.Wann2KCPSettingsDict,
                    'ui': settings.UnfoldAndInterpolateSettingsDict,
                    'ui_occ': settings.UnfoldAndInterpolateSettingsDict,
                    'ui_emp': settings.UnfoldAndInterpolateSettingsDict,
                    'w90': settings.Wannier90SettingsDict,
                    'w90_up': settings.Wannier90SettingsDict,
                    'w90_down': settings.Wannier90SettingsDict,
                    'plot': settings.PlotSettingsDict}


def sanitize_calculator_parameters(dct_in: Union[Dict[str, Dict], Dict[str, settings.SettingsDict]]) \
        -> Dict[str, settings.SettingsDict]:
    dct_out: Dict[str, settings.SettingsDict] = {}
    for k, cls in settings_classes.items():
        dct: Union[Dict, settings.SettingsDict] = dct_in.get(k, {})
        if isinstance(dct, settings.SettingsDict):
            dct_out[k] = dct
            if dct_out[k].directory == '':
                dct_out[k].directory = Path.cwd()
        else:
            dct_out[k] = cls(**dct, directory=Path.cwd())

    for k in dct_in.keys():
        if k not in settings_classes:
            raise ValueError(
                f'Unrecognized calculator_parameters entry `{k}`: valid options are `'
                + '`/`'.join(settings_classes.keys()) + '`')
    return dct_out


def header():
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
