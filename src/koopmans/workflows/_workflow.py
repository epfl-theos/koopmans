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
from typing import (Any, Callable, Dict, Generator, List, Optional, Tuple,
                    Type, TypeVar, Union)

import dill
import numpy as np
from numpy import typing as npt
from pybtex.database import BibliographyData

# isort: off
import koopmans.mpl_config
import matplotlib.pyplot as plt
# isort: on

from ase import Atoms
from ase.build.supercells import make_supercell
from ase.calculators.calculator import CalculationFailed
from ase.dft.dos import DOS
from ase.dft.kpoints import BandPath
from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.io.espresso import contruct_kcp_namelist as construct_namelist
from ase.spacegroup import symmetrize
from ase.spectrum.band_structure import BandStructure
from ase.spectrum.doscollection import GridDOSCollection
from ase.spectrum.dosdata import GridDOSData

from koopmans import calculators, outputs, settings, utils
from koopmans.bands import Bands
from koopmans.commands import ParallelCommandWithPostfix
from koopmans.kpoints import Kpoints
from koopmans.ml import AbstractMLModel, MLModel, OccEmpMLModels
from koopmans.processes import Process
from koopmans.processes.koopmans_cp import (ConvertFilesFromSpin1To2,
                                            ConvertFilesFromSpin2To1)
from koopmans.projections import ProjectionBlocks
from koopmans.pseudopotentials import (fetch_pseudo, nelec_from_pseudos,
                                       pseudo_database,
                                       pseudos_library_directory,
                                       valence_from_pseudo)
from koopmans.references import bib_data

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
    _pseudopotentials: Dict[str, str]
    pseudo_dir: Path
    projections: ProjectionBlocks
    ml_model: Optional[AbstractMLModel]
    snapshots: List[Atoms]
    version: str
    _step_counter: int

    __slots__ = utils.HasDirectory.__slots__ + ['atoms', 'parameters', 'calculator_parameters', 'name', 'kpoints',
                                                '_pseudopotentials', 'projections', 'ml_model', 'snapshots',
                                                'version', '_step_counter', 'calculations', 'processes',
                                                'steps', 'silent', 'print_indent', 'plotting', 'ml', '_bands']

    def __init__(self, atoms: Atoms,
                 pseudopotentials: Dict[str, str] = {},
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

        # Initialize the HasDirectory information (parent, base_directory, directory)
        super().__init__(parent)

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
        self.silent = False
        self.print_indent = 0
        self._bands = None

        if projections is None:
            proj_list: List[List[Any]]
            spins: List[Optional[str]]
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
            assert self.ml.estimator is not None
            assert self.ml.descriptor is not None
            if self.ml.occ_and_emp_together:
                self.ml_model = MLModel(self.ml.estimator, self.ml.descriptor)
            else:
                self.ml_model = OccEmpMLModels(self.ml.estimator, self.ml.descriptor)
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
        if pseudopotentials:
            self.pseudopotentials = pseudopotentials
        else:
            if self.parameters.pseudo_library is None:
                utils.warn(
                    'Neither a pseudopotential library nor a list of pseudopotentials was provided; defaulting to '
                    '`sg15_v1.2`')
                self.parameters.pseudo_library = 'sg15_v1.2'
            self.pseudopotentials = {}
            for symbol, tag in set([(a.symbol, a.tag) for a in self.atoms]):
                pseudo = fetch_pseudo(element=symbol, functional=self.parameters.base_functional,
                                      library=self.parameters.pseudo_library)
                if pseudo.kind == 'unknown':
                    utils.warn(f'You are using an unrecognized pseudopotential `{pseudo.name}`. Please note that '
                               'the current implementation of Koopmans functionals only supports norm-conserving '
                               'pseudopotentials.')
                elif pseudo.kind != 'norm-conserving':
                    raise ValueError('Koopmans functionals only currently support norm-conserving pseudopotentials; '
                                     f'{pseudo.name} is {pseudo.kind}')
                if tag > 0:
                    symbol += str(tag)
                self.pseudopotentials[symbol] = pseudo.name

        # Make sure calculator_parameters isn't missing any entries, and every entry corresponds to
        # settings.SettingsDict objects
        calculator_parameters = sanitize_calculator_parameters(calculator_parameters) if calculator_parameters \
            is not None else generate_default_calculator_parameters()

        # Work out the pseudopotential directory. If using a pseudo_library this is straightforward, if not...
        #  1. try to locating the directory as currently specified by the calculator
        #  2. if that fails, check if $ESPRESSO_PSEUDO is set
        #  3. if that fails, raise an error
        if self.parameters.pseudo_directory is None:
            if self.parameters.pseudo_library:
                pseudo_dir = pseudos_library_directory(self.parameters.pseudo_library, self.parameters.base_functional)
                for params in calculator_parameters.values():
                    if params.get('pseudo_dir', pseudo_dir).resolve() != pseudo_dir:
                        raise ValueError(
                            '`pseudo_dir` and `pseudo_library` are conflicting; please do not provide `pseudo_dir`')
            elif 'pseudo_dir' in calculator_parameters['kcp'] or 'pseudo_dir' in calculator_parameters['pw']:
                pseudo_dir = calculator_parameters['kcp'].get(
                    'pseudo_dir', calculator_parameters['pw'].get('pseudo_dir'))
                assert isinstance(pseudo_dir, Path)
            elif 'ESPRESSO_PSEUDO' in os.environ:
                pseudo_dir = Path(os.environ['ESPRESSO_PSEUDO'])
            else:
                pseudo_dir = Path.cwd()

            self.parameters.pseudo_directory = os.path.relpath(pseudo_dir.resolve(), self.base_directory)
        elif self.parameters.pseudo_directory.is_absolute():
            self.parameters.pseudo_directory = os.path.relpath(self.parameters.pseudo_directory, self.base_directory)

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
            nelec = nelec_from_pseudos(self.atoms, self.pseudopotentials,
                                       self.base_directory / self.parameters.pseudo_directory)
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
                        valence = valence_from_pseudo(p, self.base_directory / self.parameters.pseudo_directory)
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

        # Initialize the step counter
        self._step_counter = 0

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

    def run(self, subdirectory: Optional[str] = None, from_scratch: Optional[bool] = None) -> None:
        '''
        Run the workflow
        '''

        self.print_preamble()

        if not self.parent:
            bf = '**' if sys.stdout.isatty() else ''
            self.print(bf + self.name + bf)
            self.print(bf + '-' * len(self.name) + bf)
            if self.parameters.from_scratch:
                self._remove_tmpdirs()
            self._run_sanity_checks()

        if self.parent:
            with self._parent_context(subdirectory, from_scratch):
                self.print(f'- **{self.name}**', style='heading')
                self._run()
        else:
            if subdirectory:
                self.base_directory = Path(subdirectory).resolve()
                with utils.chdir(subdirectory):
                    self._run()
            else:
                self.base_directory = Path.cwd().resolve()
                self._run()

        self.print_conclusion()

        if not self.parent:
            self._teardown()

    @abstractmethod
    def _run(self) -> None:
        ...

    @property
    @abstractmethod
    def output_model(self) -> Type[outputs.OutputModel]:
        ...

    @property
    def pseudopotentials(self) -> OrderedDict[str, str]:
        # Always return the pseudopotentials in alphabetical order
        out = OrderedDict()
        for k in sorted(self._pseudopotentials):
            out[k] = self._pseudopotentials[k]
        return out

    @pseudopotentials.setter
    def pseudopotentials(self, value: Dict[str, str]):
        self._pseudopotentials = value

    def get_step_by_uuid(self, uuid: str):
        for step in self.steps:
            if step.uuid == uuid:
                return step
        raise ValueError(f'No step with UUID {uuid} found in workflow')

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
                             pseudopotentials=copy.deepcopy(other_wf.pseudopotentials),
                             kpoints=copy.deepcopy(other_wf.kpoints),
                             projections=other_wf.projections,
                             plotting=copy.deepcopy(other_wf.plotting),
                             ml=copy.deepcopy(other_wf.ml),
                             ml_model=other_wf.ml_model))

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

        # Check pseudopotentials exist
        if not os.path.isdir(self.base_directory / self.parameters.pseudo_directory):
            raise NotADirectoryError('The pseudopotential directory you provided '
                                     f'(`{self.base_directory / self.parameters.pseudo_directory}`) does not exist')
        if self.parameters.task != 'ui':
            for pseudo in self.pseudopotentials.values():
                if not (self.base_directory / self.parameters.pseudo_directory / pseudo).exists():
                    raise FileNotFoundError(f'`{self.base_directory / self.parameters.pseudo_directory / pseudo}` '
                                            'does not exist. Please double-check your pseudopotential settings')

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
        elif calc_type == 'wann2kc':
            calc_class = calculators.Wann2KCCalculator
        elif calc_type == 'kc_screen':
            calc_class = calculators.KoopmansScreenCalculator
        elif calc_type == 'kc_ham':
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
                else:
                    val = getattr(self, kw)
                all_kwargs[kw] = val

        # Create the calculator
        calc = calc_class(atoms=copy.deepcopy(self.atoms), **all_kwargs)

        # Link the pseudopotentials if relevant
        if calculator_parameters.is_valid('pseudo_dir'):
            for pseudo in self.pseudopotentials.values():
                self.link(None, self.base_directory / self.parameters.pseudo_directory /
                          pseudo, calc, Path('pseudopotentials') / pseudo)
            calc.parameters.pseudo_dir = 'pseudopotentials'

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

    def run_calculator(self, master_calc: calculators.Calc, enforce_spin_symmetry: bool = False):
        ''' Run a calculator.

        If enforce_spin_symmetry is True, the calculation will be run with spin symmetry enforced.
        Ultimately this wraps self.run_calculators

        :param master_calc: the calculator to run
        :param enforce_spin_symmetry: whether to enforce spin symmetry
        '''

        if enforce_spin_symmetry:
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
                prev_calc_nspin2 = master_calc.linked_files[ndr_directory][0]
                if not isinstance(prev_calc_nspin2, calculators.KoopmansCPCalculator):
                    raise ValueError(
                        'Unexpected calculator type linked during the procedure for fixing spin contamination')

                # Wipe the linked files (except for globally-linked files)
                master_calc.linked_files = {k: v for k, v in master_calc.linked_files.items() if v[0] is None}
                calc_nspin1.linked_files = {k: v for k, v in calc_nspin1.linked_files.items() if v[0] is None}

                # nspin=1 dummy
                calc_nspin1_dummy = master_calc.nspin1_dummy_calculator()
                calc_nspin1_dummy.skip_qc = True
                self.run_calculator(calc_nspin1_dummy)
                self.link(calc_nspin1_dummy, calc_nspin1_dummy.parameters.outdir,
                          calc_nspin1, calc_nspin1.parameters.outdir)

                # Copy over nspin=2 wavefunction to nspin=1 tmp directory
                process2to1 = ConvertFilesFromSpin2To1(name=None,
                                                       **prev_calc_nspin2.files_to_convert_with_spin2_to_spin1)
                self.run_process(process2to1)
                for f in process2to1.outputs.generated_files:
                    dst_dir = calc_nspin1.parameters.outdir / \
                        f'{calc_nspin1.parameters.prefix}_{calc_nspin1.parameters.ndr}.save/K00001'
                    self.link(process2to1, f, calc_nspin1, dst_dir / f.name, symlink=True, overwrite=True)

            self.run_calculator(calc_nspin1)

            # nspin=2 from scratch (dummy run for creating files of appropriate size)
            calc_nspin2_dummy = master_calc.nspin2_dummy_calculator()
            calc_nspin2_dummy.skip_qc = True
            self.run_calculator(calc_nspin2_dummy)
            self.link(calc_nspin2_dummy, calc_nspin2_dummy.parameters.outdir,
                      master_calc, master_calc.parameters.outdir)

            # Copy over nspin=1 wavefunction to nspin=2 tmp directory
            process1to2 = ConvertFilesFromSpin1To2(**calc_nspin1.files_to_convert_with_spin1_to_spin2)
            self.run_process(process1to2)

            # nspin=2, reading in the spin-symmetric nspin=1 wavefunction
            master_calc.prepare_to_read_nspin1()
            for f in process1to2.outputs.generated_files:
                dst_dir = master_calc.parameters.outdir / \
                    f'{master_calc.parameters.prefix}_{master_calc.parameters.ndr}.save/K00001'
                self.link(process1to2, f, master_calc, dst_dir / f.name, symlink=True, overwrite=True)
            self.run_calculator(master_calc)

        else:

            self.run_calculators([master_calc])

        return

    def _pre_run_calculator(self, qe_calc: calculators.Calc) -> bool:
        """Perform operations that need to occur before a calculation is run

        :param qe_calc: The calculator to run
        :return: Whether the calculation should be run
        """

        # Check that calc.directory was not set
        if qe_calc.directory_has_been_set():
            raise ValueError(
                f'`calc.directory` should not be set manually, but it was set to `{qe_calc.directory}`')

        # Prepend calc.directory with a counter
        self._step_counter += 1
        qe_calc.directory = Path(f'{self._step_counter:02}-{qe_calc.prefix}')

        # Check that another calculation hasn't already been run in this directory
        if any([calc.absolute_directory == qe_calc.absolute_directory for calc in self.calculations]):
            raise ValueError(
                f'A calculation has already been run in `{qe_calc.directory}`; this should not happen')

        # If an output file already exists, check if the run completed successfully
        if not self.parameters.from_scratch:

            calc_file = qe_calc.directory / qe_calc.prefix

            if calc_file.with_suffix(qe_calc.ext_out).is_file():

                is_complete = self.load_old_calculator(qe_calc)
                if is_complete:
                    if not self.silent:
                        self.print(
                            f'- ⏭️  Not running `{os.path.relpath(qe_calc.directory)}` as it is already complete  ')

                    # Check the convergence of the calculation
                    qe_calc.check_convergence()

                    if isinstance(qe_calc, calculators.ReturnsBandStructure):
                        qe_calc.generate_band_structure()

                    if isinstance(qe_calc, calculators.ProjwfcCalculator):
                        qe_calc.generate_dos()

                    if isinstance(qe_calc, calculators.PhCalculator):
                        qe_calc.read_dynG()

                    return False

            # If we reach this point, then this and subsequent calculations should be performed from scratch
            self.parameters.from_scratch = True

        # Remove the directory if it already exists
        if qe_calc.directory.exists():
            utils.remove(qe_calc.directory)

        # Update postfix if relevant
        if self.parameters.npool:
            if isinstance(qe_calc.command, ParallelCommandWithPostfix):
                qe_calc.command.postfix = f'-npool {self.parameters.npool}'

        return True

    def _run_calculators(self, calcs: List[calculators.Calc]) -> None:
        """Run a list of calculators, without doing anything else (other than printing messages.)

        Any pre- or post-processing of each calculator should have been done in _pre_run_calculator and
        _post_run_calculator.

        :param calcs: The calculators to run
        """

        for calc in calcs:
            assert calc.directory is not None
            dir_str = str(os.path.relpath(calc.directory))
            if sys.stdout.isatty():
                self.print(f'- 🖥️  Running `{dir_str}`...', end='\r', flush=True)

            try:
                calc.calculate()
            except CalculationFailed:
                self.print(f'- ❌ `{dir_str}` failed     ')
                raise

            if not self.silent:
                self.print(f'- ✅ `{dir_str}` completed  ')

            # If we reached here, all future calculations should be performed from scratch
            self.parameters.from_scratch = True

        return

    def _post_run_calculator(self, calc: calculators.Calc) -> None:
        """
        Perform any operations that need to be performed after a calculation has run.
        """
        # Store the calculator
        self.calculations.append(calc)
        self.steps.append(calc)

        # Ensure we inherit any modifications made to the atoms object
        if calc.atoms != self.atoms:
            self.atoms = calc.atoms

        return

    def run_calculators(self, calcs: List[calculators.Calc]):
        '''
        Run a list of *independent* calculators (default implementation is to run them in sequence)
        '''

        calcs_to_run = []
        for calc in calcs:
            calc.parent = self
            proceed = self._pre_run_calculator(calc)
            if proceed:
                calcs_to_run.append(calc)

        self._run_calculators(calcs_to_run)

        for calc in calcs_to_run:
            self._post_run_calculator(calc)

    def run_process(self, process: Process):
        '''
        Run a Process
        '''

        process.parent = self
        self._step_counter += 1
        process.directory = Path(f'{self._step_counter:02}-{process.name}')

        if not self.parameters.from_scratch and process.is_complete():
            self.print(f'- ⏭️  Not running `{os.path.relpath(process.directory)}` as it is already complete  ')
            process.load_outputs()
        else:
            if sys.stdout.isatty():
                self.print('- 🖥️  Running `' + os.path.relpath(process.directory) + '`...', end='\r')
            process.run()
            self.print('- ✅ `' + os.path.relpath(process.directory) + '` completed  ')

        self.processes.append(process)
        self.steps.append(process)

    def load_old_calculator(self, qe_calc: calculators.Calc) -> bool:
        # This is a separate function so that it can be monkeypatched by the test suite
        assert qe_calc.directory is not None
        old_calc = qe_calc.__class__.fromfile(qe_calc.directory / qe_calc.prefix)

        if old_calc.is_complete():
            # If it is complete, load the results
            qe_calc.results = old_calc.results

            # Load kpts if relevant
            if hasattr(old_calc, 'kpts'):
                qe_calc.kpts = old_calc.kpts

            self.calculations.append(qe_calc)
            self.steps.append(qe_calc)

        return old_calc.is_complete()

    def link(self, src_calc: utils.HasDirectory | None, src_path: Path | str, dest_calc: calculators.Calc,
             dest_path: Path | str, symlink=False, recursive_symlink=False, overwrite=False) -> None:
        """Link a file from one calculator to another

        Paths must be provided relative to the the calculator's directory i.e. calc.directory, unless src_calc is None

        :param src_calc: The prior calculation or process which generated the file (None if it was not generated by a
            process/calculation)
        :type src_calc: calculators.Calc | Process | None
        :param src_path: The path of the file, relative to the directory of the calculation/process
        :type src_path: Path
        :param dest_calc: The calculator which requires this file
        :type dest_calc: calculators.Calc
        :param dest_path: The path (relative to dest_calc's directory) where this file should appear
        :type dest_path: Path
        :param symlink: Whether to create a symlink instead of copying the file. Only do this if the file is not going
            to be modified by the dest_calc
        :type symlink: bool
        :param recursive_symlink: If true, when linking a folder, link its contents individually and not the folder
            itself
        :type recursive_symlink: bool
        :param overwrite: If true, allow linking to a destination that has already been linked to, in which case the
            link will be overwritten
        :type overwrite: bool
        """

        if isinstance(src_path, str):
            src_path = Path(src_path)
        if isinstance(dest_path, str):
            dest_path = Path(dest_path)

        if src_path.is_absolute() and src_calc is not None:
            raise ValueError(f'`src_path` in `{self.__class__.__name__}.link()` must be a relative path if a '
                             f'`src_calc` is provided')
        if src_calc is None and not src_path.is_absolute():
            raise ValueError(f'`src_path` in `{self.__class__.__name__}.link()` must be an absolute path if '
                             '`src_calc` is not provided')
        if dest_path.is_absolute():
            raise ValueError(f'`dest_path` in `{self.__class__.__name__}.link()` must be a relative path')

        dest_calc.link_file(src_calc, src_path, dest_path, symlink=symlink,  # type: ignore
                            recursive_symlink=recursive_symlink, overwrite=overwrite)

    def print(self, text: str = '', style='body', parse_asterisks=True, flush=True, wrap=True, **kwargs: Any):
        if sys.stdout.isatty():
            if style == 'heading' and '**' not in text:
                text = f'**{text}**'
            if parse_asterisks:
                while '**' in text:
                    text = text.replace('**', '\033[1m', 1).replace('**', '\033[0m', 1)
                while '*' in text:
                    text = text.replace('*', '\033[3m', 1).replace('*', '\033[0m', 1)
        if style == 'body':
            utils.indented_print(text, self.print_indent + 2, flush=flush, wrap=wrap, **kwargs)
        elif style == 'heading':
            assert kwargs.get('end', '\n') == '\n'
            utils.indented_print(text, self.print_indent, flush=flush, wrap=wrap, **kwargs)
        else:
            raise ValueError(f'Invalid choice `{style}` for style; must be `heading`/`body`')

    @contextmanager
    def _parent_context(self, subdirectory: Optional[str] = None,
                        from_scratch: Optional[bool] = None) -> Generator[None, None, None]:
        '''
        Context for calling self._run(), within which self inherits relevant information from self.parent, runs, and
        then passes back relevant information to self.parent
        '''

        assert self.parent is not None

        # Increase the indent level
        self.print_indent = self.parent.print_indent + 2

        # Ensure altering self.calculator_parameters won't affect self.parent.calculator_parameters
        if self.calculator_parameters is self.parent.calculator_parameters:
            self.calculator_parameters = copy.deepcopy(self.parent.calculator_parameters)

        # Setting from_scratch to a non-None value will override the value of subworkflow.from_scratch...
        if from_scratch is None:
            self.parameters.from_scratch = self.parent.parameters.from_scratch
        else:
            self.parameters.from_scratch = from_scratch

        # Link the lists of calculations and processes
        self.calculations = self.parent.calculations
        self.processes = self.parent.processes

        # Link the bands
        if self.parent.bands is not None:
            self.bands = self.parent.bands

        subdirectory = self.name.replace(' ', '-').lower() if subdirectory is None else subdirectory
        try:
            # Prepend the step counter to the subdirectory name
            self.parent._step_counter += 1
            subdirectory_path = Path(f'{self.parent._step_counter:02}-{subdirectory}')
            if self.parameters.from_scratch and subdirectory_path.is_dir():
                shutil.rmtree(subdirectory_path)

            # Store the directory
            self.directory = subdirectory_path

            # Run the workflow
            with utils.chdir(subdirectory_path):
                yield
        finally:
            # ... and will prevent inheritance of from_scratch
            if from_scratch is None:
                self.parent.parameters.from_scratch = self.parameters.from_scratch

            if self.bands is not None and self.parent.bands is None:
                # Copy the entire bands object
                self.parent.bands = self.bands

        self.parent.steps.append(self.steps)

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
                 pseudopotentials=dct.pop('_pseudopotentials'),
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

        # Third, flatten the kc_wann subdicts
        kc_wann_blocks = calcdict.pop('kc_wann', {'kc_ham': {}, 'kc_screen': {}, 'wann2kc': {}})
        calcdict.update(**kc_wann_blocks)

        # Finally, generate a SettingsDict for every single kind of calculator, regardless of whether or not there was
        # a corresponding block in the json file
        calculator_parameters = {}
        w90_block_projs: List = []
        w90_block_spins: List[Union[str, None]] = []
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

            psp_lib = self.parameters.pseudo_library
            if psp_lib is not None:
                psp_subset = [p for p in pseudo_database if p.functional == self.parameters.base_functional
                              and p.library == psp_lib]
                citations = set([c for psp in psp_subset for c in psp.citations if psp.element in self.atoms.symbols])

                for citation in citations:
                    add_ref(citation, f'Citation for the {psp_lib.replace("_", " ")} pseudopotential library')

        utils.print_alert(
            'note', f'Please cite the papers listed in `{self.name}.bib` in work involving this calculation')
        relevant_references.to_file(self.name + '.bib')

    def print_preamble(self):
        if self.parent:
            return

        self.print(header())

        self.print_bib()

    def print_conclusion(self):
        from koopmans.io import write

        if self.parent:
            return

        # Save workflow to file
        write(self, self.name + '.pkl')

        # Save the ML model to a separate file
        if self.ml.train:
            assert self.ml_model is not None
            write(self.ml_model, self.name + '_ml_model.pkl')

        # Print farewell message
        self.print('\n**Workflow complete** 🎉')

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

            elif code in ['pw2wannier', 'wann2kc', 'kc_screen', 'kc_ham', 'projwfc', 'wann2kcp', 'ph']:
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
            bsplot_kwargs = [bsplot_kwargs]
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
        utils.savefig(fname=filename + '.png', bbox_extra_artists=legends, bbox_inches='tight')

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

        # Removing tmpdirs
        if not self.parameters.keep_tmpdirs:
            self._remove_tmpdirs()

    def number_of_electrons(self, spin: Optional[str] = None) -> int:
        # Return the number of electrons in a particular spin channel
        nelec_tot = nelec_from_pseudos(self.atoms, self.pseudopotentials,
                                       self.base_directory / self.parameters.pseudo_directory)
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
            'kc_ham': settings.KoopmansHamSettingsDict(),
            'kc_screen': settings.KoopmansScreenSettingsDict(),
            'ph': settings.PhSettingsDict(),
            'projwfc': settings.ProjwfcSettingsDict(),
            'pw': settings.PWSettingsDict(),
            'pw2wannier': settings.PW2WannierSettingsDict(),
            'wann2kcp': settings.Wann2KCPSettingsDict(),
            'ui': settings.UnfoldAndInterpolateSettingsDict(),
            'ui_occ': settings.UnfoldAndInterpolateSettingsDict(),
            'ui_emp': settings.UnfoldAndInterpolateSettingsDict(),
            'wann2kc': settings.Wann2KCSettingsDict(),
            'w90': settings.Wannier90SettingsDict(),
            'w90_up': settings.Wannier90SettingsDict(),
            'w90_down': settings.Wannier90SettingsDict(),
            'plot': settings.PlotSettingsDict()}


# Define which function to use to read each block
settings_classes = {'kcp': settings.KoopmansCPSettingsDict,
                    'kc_ham': settings.KoopmansHamSettingsDict,
                    'kc_screen': settings.KoopmansScreenSettingsDict,
                    'wann2kc': settings.Wann2KCSettingsDict,
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
