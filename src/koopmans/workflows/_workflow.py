"""
Generic workflow object for koopmans
Written by Edward Linscott Oct 2020
Converted workflows from functions to objects Nov 2020
"""

from __future__ import annotations

import copy
import json as json_ext
import operator
import os
import shutil
import subprocess
from abc import ABC, abstractmethod
from contextlib import contextmanager
from functools import reduce
from pathlib import Path
from types import ModuleType
from typing import (Any, Callable, Dict, Generator, List, Optional, Type,
                    TypeVar, Union)

import numpy as np
from numpy import typing as npt
from pybtex.database import BibliographyData

# isort: off
import koopmans.mpl_config
import matplotlib.pyplot as plt
# isort: on

import ase
from ase import Atoms
from ase.build.supercells import make_supercell
from ase.calculators.calculator import CalculationFailed
from ase.dft.dos import DOS
from ase.dft.kpoints import BandPath
from ase.io.espresso import contruct_kcp_namelist as construct_namelist
from ase.spacegroup import symmetrize
from ase.spectrum.band_structure import BandStructure
from ase.spectrum.doscollection import GridDOSCollection
from ase.spectrum.dosdata import GridDOSData

from koopmans import calculators, settings, utils
from koopmans.bands import Bands
from koopmans.commands import ParallelCommandWithPostfix
from koopmans.kpoints import Kpoints
from koopmans.ml import MLModel
from koopmans.projections import ProjectionBlocks
from koopmans.pseudopotentials import (fetch_pseudo, nelec_from_pseudos,
                                       pseudo_database,
                                       pseudos_library_directory,
                                       valence_from_pseudo)
from koopmans.references import bib_data

T = TypeVar('T', bound='calculators.CalculatorExt')
W = TypeVar('W', bound='Workflow')


class Workflow(ABC):

    r'''
    Abstract base class that defines a Koopmans workflow
    Parameters
    ----------
    atoms : Atoms
        an ASE ``Atoms`` object defining the atomic positions, cell, etc
    pseudopotentals : Dict[str, str]
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
    parent: Optional[Workflow]

    def __init__(self, atoms: Atoms,
                 pseudopotentials: Dict[str, str] = {},
                 kpoints: Optional[Kpoints] = None,
                 projections: Optional[ProjectionBlocks] = None,
                 name: str = 'koopmans_workflow',
                 parameters: Union[Dict[str, Any], settings.WorkflowSettingsDict] = {},
                 calculator_parameters: Optional[Union[Dict[str, Dict[str, Any]],
                                                       Dict[str, settings.SettingsDict]]] = None,
                 plotting: Union[Dict[str, Any], settings.PlotSettingsDict] = {},
                 ml: Union[Dict[str, Any], settings.MLSettingsDict] = {},
                 autogenerate_settings: bool = True,
                 **kwargs: Dict[str, Any]):

        # Parsing parameters
        self.parameters = settings.WorkflowSettingsDict(**parameters)
        for key, value in kwargs.items():
            if self.parameters.is_valid(key):
                self.parameters[key] = value
        self.atoms: Atoms = atoms
        self.name = name
        self.calculations: List[calculators.Calc] = []
        self.silent = False
        self.print_indent = 1

        if projections is None:
            proj_list: List[List[Any]]
            spins: List[Optional[str]]
            if self.parameters.spin_polarized:
                proj_list = [[], [], [], []]
                fillings = [True, True, False, False]
                spins = ['up', 'down', 'up', 'down']
            else:
                proj_list = [[], []]
                fillings = [True, False]
                spins = [None, None]
            self.projections = ProjectionBlocks.fromprojections(
                proj_list, fillings=fillings, spins=spins, atoms=self.atoms)
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
        if self.ml.use_ml:
            if self.ml.occ_and_emp_together:
                self.ml.ml_model = MLModel(self.ml.type_of_ml_model)
            else:
                self.ml.ml_model_occ = MLModel(self.ml.type_of_ml_model)
                self.ml.ml_model_emp = MLModel(self.ml.type_of_ml_model)

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
                    'sg15_v1.2')
                self.parameters.pseudo_library = 'sg15_v1.2'
            self.pseudopotentials = {}
            for symbol, tag in set([(a.symbol, a.tag) for a in self.atoms]):
                pseudo = fetch_pseudo(element=symbol, functional=self.parameters.base_functional,
                                      library=self.parameters.pseudo_library)
                if pseudo.kind == 'unknown':
                    utils.warn(f'You are using an unrecognized pseudopotential {pseudo.name}. Please note that '
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
                            '"pseudo_dir" and "pseudo_library" are conflicting; please do not provide "pseudo_dir"')
            elif 'pseudo_dir' in calculator_parameters['kcp'] or 'pseudo_dir' in calculator_parameters['pw']:
                pseudo_dir = calculator_parameters['kcp'].get(
                    'pseudo_dir', calculator_parameters['pw'].get('pseudo_dir'))
                assert isinstance(pseudo_dir, Path)
            elif 'ESPRESSO_PSEUDO' in os.environ:
                pseudo_dir = Path(os.environ['ESPRESSO_PSEUDO'])
            else:
                pseudo_dir = Path.cwd()

            self.parameters.pseudo_directory = pseudo_dir.resolve()

        # Before saving the calculator_parameters, automatically generate some keywords and perform some sanity checks
        if self.parameters.task != 'ui' and autogenerate_settings:
            # Automatically calculate nelec/nelup/neldw/etc using information contained in the pseudopotential files
            # and the kcp settings
            nelec = nelec_from_pseudos(self.atoms, self.pseudopotentials, self.parameters.pseudo_directory)

            tot_charge = calculator_parameters['kcp'].get('tot_charge', 0)
            nelec -= tot_charge
            tot_mag = calculator_parameters['kcp'].get('tot_magnetization', nelec % 2)
            nelup = int(nelec / 2 + tot_mag / 2)
            neldw = int(nelec / 2 - tot_mag / 2)

            # Setting up the magnetic moments
            if 'starting_magnetization(1)' in calculator_parameters['kcp']:
                labels = [s + str(t) if t > 0 else s for s, t in zip(atoms.symbols, atoms.get_tags())]
                starting_magmoms = {}
                for i, (l, p) in enumerate(self.pseudopotentials.items()):
                    # ASE uses absolute values; QE uses the fraction of the valence
                    frac_mag = calculator_parameters['kcp'].pop(f'starting_magnetization({i + 1})', 0.0)
                    valence = valence_from_pseudo(p, self.parameters.pseudo_directory)
                    starting_magmoms[l] = frac_mag * valence
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
                    raise ValueError(f'You have provided projection information in the calculator_parameters[{block}] '
                                     f'argument to {self.__class__.__name__}. Please instead specify projections '
                                     'via the "projections" argument')
                for kw in ['exclude_bands', 'num_wann', 'num_bands', 'projections']:
                    if kw in params:
                        utils.warn(f'{kw} will be overwritten by the workflow; it is best to not specify this keyword '
                                   'and to instead double-check the keyword in the various .win files '
                                   'generated by the workflow.')

            # Automatically parsing algebraic settings. We need to provide nelec explicitly for calculators such as
            # PWCalculators, which don't have nelec as a setting but it is still useful for specifying settings
            # algebraically
            params.parse_algebraic_settings(nelec=nelec)

            # Store the sanitized parameters
            self.calculator_parameters[block] = params

        # If atoms has a calculator, overwrite the kpoints and pseudopotentials variables and then detach the calculator
        if atoms.calc is not None:
            utils.warn(f'You have initialized a {self.__class__.__name__} object with an atoms object that possesses '
                       'a calculator. This calculator will be ignored.')
            self.atoms.calc = None

        # Initialize self.parent
        self.parent: Optional[Workflow] = None

        # For any kwargs...
        for key, value in kwargs.items():
            match = False
            # if they correspond to any valid calculator parameter, set it
            for calc_params in self.calculator_parameters.values():
                if calc_params.is_valid(key):
                    calc_params[key] = value
                    match = True
            # if not a calculator, workflow, or plotting keyword, raise an error
            if not match and not self.parameters.is_valid(key) and not self.plotting.is_valid(key) and not self.ml.is_valid(key):
                raise ValueError(f'{key} is not a valid setting')

    def __eq__(self, other: Any):
        if isinstance(other, Workflow):
            return self.__dict__ == other.__dict__
        return False

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

    def run(self, subdirectory: Optional[Union[str, Path]] = None, from_scratch: Optional[bool] = None) -> None:
        '''
        Run the workflow
        '''

        self.print_preamble()

        if not self.parent:
            if self.parameters.from_scratch:
                self._remove_tmpdirs()
            self._run_sanity_checks()

        if self.parent:
            with self._parent_context(subdirectory, from_scratch):
                self._run()
        else:
            self._run()

        self.print_conclusion()

        if not self.parent:
            self._teardown()

    @abstractmethod
    def _run(self) -> None:
        ...

    @property
    def pseudopotentials(self) -> Dict[str, str]:
        return self._pseudopotentials

    @pseudopotentials.setter
    def pseudopotentials(self, value: Dict[str, str]):
        self._pseudopotentials = value

    @classmethod
    def fromparent(cls: Type[W], parent_wf: Workflow, **kwargs: Any) -> W:
        '''
        Creates a subworkflow with the same configuration as the parent workflow
        e.g.
        >>> sub_wf = Workflow.fromparent(self)
        '''

        parameters = copy.deepcopy(parent_wf.parameters)
        parameter_kwargs = {k: v for k, v in kwargs.items() if parameters.is_valid(k)}
        other_kwargs = {k: v for k, v in kwargs.items() if not parameters.is_valid(k)}
        parameters.update(**parameter_kwargs)

        child_wf = cls(atoms=copy.deepcopy(parent_wf.atoms),
                       parameters=parameters,
                       calculator_parameters=copy.deepcopy(parent_wf.calculator_parameters),
                       name=copy.deepcopy(parent_wf.name),
                       pseudopotentials=copy.deepcopy(parent_wf.pseudopotentials),
                       kpoints=copy.deepcopy(parent_wf.kpoints),
                       projections=copy.deepcopy(parent_wf.projections),
                       plotting=copy.deepcopy(parent_wf.plotting),
                       ml=copy.deepcopy(parent_wf.ml),
                       **other_kwargs)
        child_wf.parent = parent_wf
        return child_wf

    def _run_sanity_checks(self):
        # Check internal consistency of workflow settings
        if self.parameters.fix_spin_contamination is None:
            self.parameters.fix_spin_contamination = not self.parameters.spin_polarized
        else:
            if self.parameters.fix_spin_contamination and self.parameters.spin_polarized:
                raise ValueError('fix_spin_contamination = True is incompatible with spin_polarized = True')

        if self.parameters.method == 'dfpt':
            if self.parameters.frozen_orbitals is None:
                self.parameters.frozen_orbitals = True
            if not self.parameters.frozen_orbitals:
                raise ValueError('"frozen_orbitals" must be equal to True when "method" is "dfpt"')
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
                    utils.warn('eps_inf missing in input; it will default to 1.0. Proceed with caution for periodic '
                               'systems; consider setting eps_inf == "auto" to calculate it automatically.')
                    self.parameters.eps_inf = 1.0

            if self.parameters.mt_correction is None:
                self.parameters.mt_correction = False
            if self.parameters.mt_correction:
                raise ValueError('Do not use Martyna-Tuckerman corrections for periodic systems')

            # Check the value of eps_inf
            if self.parameters.eps_inf:
                if isinstance(self.parameters.eps_inf, float) and self.parameters.eps_inf < 1.0:
                    raise ValueError('eps_inf cannot be lower than 1.0')

            # Check symmetry of the system
            dataset = symmetrize.check_symmetry(self.atoms, 1e-6, verbose=False)
            if dataset['number'] not in range(195, 231):
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
                raise ValueError(f'In order to use init_orbitals={self.parameters.init_orbitals}, projections must be '
                                 'provided')
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
            utils.warn('In KCP do_wf_cmplx = False is not consistent with gamma_only = False. '
                       'Changing do_wf_cmplx to True')
            self.calculator_parameters['kcp'].do_wf_cmplx = True

        # Check pseudopotentials exist
        if not os.path.isdir(self.parameters.pseudo_directory):
            raise NotADirectoryError(
                f'The pseudopotential directory you provided ({self.parameters.pseudo_directory}) does not exist')
        if self.parameters.task != 'ui':
            for pseudo in self.pseudopotentials.values():
                if not (self.parameters.pseudo_directory / pseudo).exists():
                    raise FileNotFoundError(f'{self.parameters.pseudo_directory / pseudo} does not exist. Please '
                                            'double-check your pseudopotential settings')

        # Make sanity checks for the ML model
        if self.ml.use_ml:
            utils.warn("Predicting screening parameters with machine-learning is an experimental feature; proceed with caution")
            if self.parameters.task not in ['trajectory', 'convergence_ml']:
                raise NotImplementedError(
                    f'Using the ML-prediction for the {self.parameter.task}-task has not yet been implemented.')
            if self.parameters.method != 'dscf':
                raise NotImplementedError(
                    f"Using the ML-prediction for the {self.parameters.method}-method has not yet been implemented")
            if self.parameters.functional != 'ki':
                raise NotImplementedError(
                    f'Using the ML-prediction for the {self.parameters.functional}-functional has not yet been implemented.')
            if self.parameters.init_orbitals != 'mlwfs':
                raise NotImplementedError(
                    f'Using the ML-prediction for {self.parameters.init_orbitals}-init orbitals has not yet been implemented.')
            if self.parameters.init_empty_orbitals != self.parameters.init_orbitals:
                raise NotImplementedError(
                    f'Using the ML-prediction for using different init orbitals for empty states than for occupied states has not yet been implemented.')
            if self.parameters.spin_polarized:
                utils.warn(f'Using the ML-prediction for spin-polarised systems has not yet been extensively tested.')
            if not all(self.atoms.pbc):
                utils.warn(f'Using the ML-prediction for non-periodic systems has not yet been extensively tested.')
            if self.parameters.orbital_groups:
                utils.warn('Using orbital_groups has not yet been extensively tested.')
            if not np.all(self.atoms.cell.angles() == 90.0):
                raise ValueError(f"The ML-workflow has only been implemented for simulation cells that have 90° angles")

            if self.parameters.task is not 'convergence_ml':
                assert isinstance(self.ml.n_max, int)
                assert isinstance(self.ml.l_max, int)
                assert isinstance(self.ml.r_min, float)
                assert isinstance(self.ml.r_max, float)

            # convert now each parameter to a list to be able to run the same checks irregardless of the task

            def convert_to_list(param, type):
                if isinstance(param, type):  # if param is an int or a float convert it for the checks to a list
                    return [param]
                else:  # if param is not an int or a float check that it is a list of ints / floats
                    assert(isinstance(param, list))
                    for value in param:
                        assert(isinstance(value, type))
                    return param

            n_maxs = convert_to_list(self.ml.n_max, int)
            l_maxs = convert_to_list(self.ml.l_max, int)
            r_mins = convert_to_list(self.ml.r_min, float)
            r_maxs = convert_to_list(self.ml.r_max, float)

            # check that each n_max, l_max, r_max and r_min are greater or equal to 0 and that r_min is smaller than r_max
            for n_max in n_maxs:
                if not n_max > 0:
                    raise ValueError(f"n_max has to be larger than zero. The provided value is n_max={n_max}")
            for l_max in l_maxs:
                if not l_max >= 0:
                    raise ValueError(f"l_max has to be equal or larger than zero. The provided value is l_max={l_max}")
            for r_min in r_mins:
                if not r_min >= 0:
                    raise ValueError(f"r_min has to be equal or larger than zero. The provided value is r_min={r_min}")
                if r_min < 0.5:
                    utils.warn(
                        f"Small values of r_min (<0.5) can lead to problems in the construction of the radial basis. The provided value is r_min={r_min}.")
            for r_max in r_maxs:
                if not any(r_min < r_max for r_min in r_mins):
                    raise ValueError(f"All provided values of r_min are larger or equal to r_max={r_max}.")

            # for the convergence_ml task we want to have each parameter in a list form
            if self.parameters.task == 'convergence_ml':

                implemented_quantities_of_interest = ['alphas', 'evs']
                self.ml.quantities_of_interest = convert_to_list(self.ml.quantities_of_interest, str)

                for qoi in self.ml.quantities_of_interest:
                    if qoi not in implemented_quantities_of_interest:
                        raise NotImplementedError(
                            "Performing the convergence_analysis w.r.t. {qoi} has not yet been implement.")

    def new_calculator(self,
                       calc_type: str,
                       directory: Optional[Path] = None,
                       kpts: Optional[Union[List[int], BandPath]] = None,
                       **kwargs) -> T:  # type: ignore[type-var, misc]

        calc_class: Type[T]

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
        elif calc_type.startswith('ui'):
            calc_class = calculators.UnfoldAndInterpolateCalculator
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
            raise ValueError(f'Cound not find a calculator of type {calc_type}')

        # Merge calculator_parameters and kwargs, giving kwargs higher precedence
        all_kwargs: Dict[str, Any] = {}
        calculator_parameters = self.calculator_parameters[calc_type]
        all_kwargs.update(**calculator_parameters)
        all_kwargs.update(**kwargs)

        # For the k-points, the Workflow has two options: self.kpoints.grid and self.kpoints.path. A calculator should
        # only ever have one of these two. By default, use the grid.
        if 'kpts' in calculator_parameters.valid:
            all_kwargs['kpts'] = kpts if kpts is not None else self.kpoints.grid

        # Add further information to the calculator as required
        for kw in ['pseudopotentials', 'pseudo_dir', 'gamma_only', 'kgrid', 'kpath', 'koffset', 'plotting']:
            if kw not in all_kwargs and calculator_parameters.is_valid(kw):
                val: Any
                if kw == 'kgrid':
                    val = self.kpoints.grid
                elif kw == 'kpath':
                    val = self.kpoints.path
                elif kw == 'koffset':
                    val = self.kpoints.offset
                elif kw == 'gamma_only':
                    val = self.kpoints.gamma_only
                elif kw == 'pseudo_dir':
                    val = self.parameters.pseudo_directory
                else:
                    val = getattr(self, kw)
                all_kwargs[kw] = val

        # Create the calculator
        calc = calc_class(atoms=copy.deepcopy(self.atoms), **all_kwargs)

        # Add the directory if provided
        if directory is not None:
            calc.directory = directory

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

    def run_calculator(self, master_qe_calc: calculators.Calc, enforce_ss=False):
        '''
        Wrapper for run_calculator_single that manages the optional enforcing of spin symmetry
        '''

        if enforce_ss:
            if not isinstance(master_qe_calc, calculators.CalculatorCanEnforceSpinSym):
                raise NotImplementedError(f'{master_qe_calc.__class__.__name__} cannot enforce spin symmetry')

            if not master_qe_calc.from_scratch:
                # PBE with nspin=1 dummy
                qe_calc = master_qe_calc.nspin1_dummy_calculator()
                qe_calc.skip_qc = True
                self.run_calculator_single(qe_calc)
                # Copy over nspin=2 wavefunction to nspin=1 tmp directory (if it has not been done already)
                if self.parameters.from_scratch:
                    master_qe_calc.convert_wavefunction_2to1()

            # PBE with nspin=1
            qe_calc = master_qe_calc.nspin1_calculator()
            self.run_calculator_single(qe_calc)

            # PBE from scratch with nspin=2 (dummy run for creating files of appropriate size)
            qe_calc = master_qe_calc.nspin2_dummy_calculator()
            qe_calc.skip_qc = True
            self.run_calculator_single(qe_calc)

            # Copy over nspin=1 wavefunction to nspin=2 tmp directory (if it has not been done already)
            if self.parameters.from_scratch:
                master_qe_calc.convert_wavefunction_1to2()

            # PBE with nspin=2, reading in the spin-symmetric nspin=1 wavefunction
            master_qe_calc.prepare_to_read_nspin1()
            self.run_calculator_single(master_qe_calc)

        else:

            self.run_calculator_single(master_qe_calc)

        return

    def run_calculator_single(self, qe_calc: calculators.Calc):
        # Runs qe_calc.calculate with additional checks

        # If an output file already exists, check if the run completed successfully
        verb = 'Running'
        if not self.parameters.from_scratch:

            calc_file = qe_calc.directory / qe_calc.prefix

            if calc_file.with_suffix(qe_calc.ext_out).is_file():
                verb = 'Rerunning'

                is_complete = self.load_old_calculator(qe_calc)
                if is_complete:
                    if not self.silent:
                        self.print(f'Not running {os.path.relpath(calc_file)} as it is already complete')

                    # Check the convergence of the calculation
                    qe_calc.check_convergence()

                    if isinstance(qe_calc, calculators.ReturnsBandStructure):
                        qe_calc.generate_band_structure()

                    if isinstance(qe_calc, calculators.ProjwfcCalculator):
                        qe_calc.generate_dos()

                    if isinstance(qe_calc, calculators.PhCalculator):
                        qe_calc.read_dynG()
                    return

        if not self.silent:
            dir_str = os.path.relpath(qe_calc.directory) + '/'
            self.print(f'{verb} {dir_str}{qe_calc.prefix}...', end='', flush=True)

        # Update postfix if relevant
        if self.parameters.npool:
            if isinstance(qe_calc.command, ParallelCommandWithPostfix):
                qe_calc.command.postfix = f'-npool {self.parameters.npool}'

        try:
            qe_calc.calculate()
        except CalculationFailed:
            self.print(' failed')
            raise

        if not self.silent:
            self.print(' done')

        # Store the calculator
        self.calculations.append(qe_calc)

        # If we reached here, all future calculations should be performed from scratch
        self.parameters.from_scratch = True

        return

    def load_old_calculator(self, qe_calc: calculators.Calc) -> bool:
        # This is a separate function so that it can be monkeypatched by the test suite
        old_calc = qe_calc.__class__.fromfile(qe_calc.directory / qe_calc.prefix)

        if old_calc.is_complete():
            # If it is complete, load the results
            qe_calc.results = old_calc.results

            # Load kpts if relevant
            if hasattr(old_calc, 'kpts'):
                qe_calc.kpts = old_calc.kpts

            # Load bandstructure if present, too
            if isinstance(qe_calc, calculators.UnfoldAndInterpolateCalculator):
                qe_calc.read_bands()
                # If the band structure file does not exist, we must re-run
                if 'band structure' not in qe_calc.results:
                    return False

            self.calculations.append(qe_calc)

        return old_calc.is_complete()

    def print(self, text: str = '', style: str = 'body', **kwargs: Any):
        if style == 'body':
            utils.indented_print(str(text), self.print_indent + 1, **kwargs)
        else:
            if style == 'heading':
                underline = '='
            elif style == 'subheading':
                underline = '-'
            else:
                raise ValueError(f'Invalid choice "{style}" for style; must be heading/subheading/body')
            assert kwargs.get('end', '\n') == '\n'
            utils.indented_print()
            utils.indented_print(str(text), self.print_indent, **kwargs)
            utils.indented_print(underline * len(text), self.print_indent, **kwargs)

    @contextmanager
    def _parent_context(self, subdirectory: Optional[Union[str, Path]] = None,
                        from_scratch: Optional[bool] = None) -> Generator[None, None, None]:
        '''
        Context for calling self._run(), within which self inherits relevant information from self.parent, runs, and
        then passes back relevant information to self.parent
        '''

        assert self.parent is not None

        # Automatically pass along the name of the overall workflow
        if self.name is None:
            self.name = self.parent.name

        # Increase the indent level
        self.print_indent = self.parent.print_indent + 1

        # Ensure altering self.calculator_parameters won't affect self.parent.calculator_parameters
        if self.calculator_parameters is self.parent.calculator_parameters:
            self.calculator_parameters = copy.deepcopy(self.parent.calculator_parameters)

        # Setting from_scratch to a non-None value will override the value of subworkflow.from_scratch...
        if from_scratch is None:
            self.parameters.from_scratch = self.parent.parameters.from_scratch
        else:
            self.parameters.from_scratch = from_scratch

        # Link the list of calculations
        self.calculations = self.parent.calculations

        # Link the ML_Model
        if self.ml.use_ml:
            if self.parent.ml.occ_and_emp_together:
                self.ml.ml_model = self.parent.ml.ml_model
            else:
                self.ml.ml_model_occ = self.parent.ml.ml_model_occ
                self.ml.ml_model_emp = self.parent.ml.ml_model_emp

        # Link the bands
        if hasattr(self.parent, 'bands'):
            self.bands = copy.deepcopy(self.parent.bands)
            # Only include the most recent screening parameter and wipe the error history
            for b in self.bands:
                if len(b.alpha_history) > 0:
                    b.alpha_history = [b.alpha]
                b.error_history = []

        try:
            if subdirectory is not None:
                # Ensure subdirectory is a Path
                subdirectory = Path(subdirectory)
                # Update directories
                for key in self.calculator_parameters.keys():
                    params = self.calculator_parameters[key]
                    for setting in params.are_paths:
                        if setting == 'pseudo_dir':
                            continue
                        path = getattr(params, setting, None)
                        if path is not None and Path.cwd() in path.parents:
                            new_path = subdirectory.resolve() / os.path.relpath(path)
                            setattr(params, setting, new_path)

                # Run the workflow
                with utils.chdir(subdirectory):
                    yield
            else:
                yield
        finally:
            # ... and will prevent inheritance of from_scratch
            if from_scratch is None:
                self.parent.parameters.from_scratch = self.parameters.from_scratch

            # Copy back over the bands
            if hasattr(self, 'bands'):
                if hasattr(self.parent, 'bands'):
                    # Add the alpha and error history
                    for b, b_sub in zip(self.parent.bands, self.bands):
                        b.alpha_history += b_sub.alpha_history[1:]
                        b.error_history += b_sub.error_history
                else:
                    # Copy the entire bands object
                    self.parent.bands = self.bands

            # Make sure any updates to the projections are passed along
            self.parent.projections = self.projections

    def todict(self):
        # Shallow copy
        dct = dict(self.__dict__)

        # Adding information required by the json decoder
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @classmethod
    def fromdict(cls, dct: Dict[str, Any], **kwargs) -> Workflow:
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
    def bands(self):
        if not hasattr(self, '_bands'):
            raise AttributeError('Bands have not been initialized')
        return self._bands

    @bands.setter
    def bands(self, value: Bands):
        assert isinstance(value, Bands)
        self._bands = value

    @classmethod
    def fromjson(cls, fname: str, override: Dict[str, Any] = {}):

        with open(fname, 'r') as fd:
            bigdct = json_ext.loads(fd.read())
        wf = cls._fromjsondct(bigdct, override)

        # Define the name of the workflow using the name of the json file
        wf.name = fname.replace('.json', '')
        return wf

    @classmethod
    def _fromjsondct(cls, bigdct: Dict[str, Any], override: Dict[str, Any] = {}):

        # Override all keywords provided explicitly
        utils.update_nested_dict(bigdct, override)

        kwargs: Dict[str, Any] = {}

        # Loading atoms object
        atoms_dict = bigdct.pop('atoms', None)
        if atoms_dict:
            atoms = read_atoms_dict(utils.parse_dict(atoms_dict))
        else:
            raise ValueError('Please provide an "atoms" block in the json input file')

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

        # First, extract the nested w90 subdictionaries
        if 'w90' in calcdict:
            for filling in ['occ', 'emp']:
                for spin in ['up', 'down']:
                    # Add any keywords in the filling:spin subsubdictionary
                    subsubdct = calcdict['w90'].get(filling, {}).get(spin, {})
                    calcdict[f'w90_{filling}_{spin}'] = subsubdct
                    # Add any keywords in the filling subdictionary
                    subdct = {k: v for k, v in calcdict['w90'].get(filling, {}).items() if k not in ['up', 'down']}
                    calcdict[f'w90_{filling}_{spin}'].update(subdct)
                    # Add any keywords in the main dictionary
                    dct = {k: v for k, v in calcdict['w90'].items() if k not in ['occ', 'emp']}
                    calcdict[f'w90_{filling}_{spin}'].update(dct)
                # Also create a spin-independent set of parameters
                calcdict[f'w90_{filling}'] = {}
                calcdict[f'w90_{filling}'].update(subdct)
                calcdict[f'w90_{filling}'].update(dct)
            # Finally, remove the nested w90 entry
            del calcdict['w90']

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
        w90_block_filling: List[bool] = []
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
                            'sc_dim in the UI block should match the kpoints provided in the kpoints block')
                    dct.pop('sc_dim')
                if 'kpath' in dct and kpts.path is not None:
                    if kpts.path != dct['kpath']:
                        raise ValueError('kpath in the UI block should match that provided in the kpoints block')
                    dct.pop('kpath')
            elif block.startswith('w90'):
                # If we are spin-polarized, don't store the spin-independent w90 block
                # Likewise, if we are not spin-polarized, don't store the spin-dependent w90 blocks
                if parameters.spin_polarized is not ('up' in block or 'down' in block):
                    continue
                if 'projections' in dct and 'projections_blocks' in dct:
                    raise ValueError(f'You have provided both "projections" and "projections_block" for {block} but '
                                     'these keywords are mutually exclusive')
                elif 'projections_blocks' in dct:
                    projs = dct.pop('projections_blocks')
                else:
                    projs = [dct.pop('projections', [])]
                w90_block_projs += projs
                w90_block_filling += ['occ' in block for _ in range(len(projs))]
                if 'up' in block:
                    w90_block_spins += ['up' for _ in range(len(projs))]
                elif 'down' in block:
                    w90_block_spins += ['down' for _ in range(len(projs))]
                else:
                    w90_block_spins += [None for _ in range(len(projs))]

            calculator_parameters[block] = settings_class(**dct)

        # Adding the projections to the workflow kwargs (this is unusual in that this is an attribute of the workflow
        # object but it is provided in the w90 subdictionary)
        kwargs['projections'] = ProjectionBlocks.fromprojections(
            w90_block_projs, w90_block_filling, w90_block_spins, atoms)

        # Check for unexpected blocks
        for block in bigdct:
            raise ValueError(f'Unrecognized block "{block}" in the json input file')

        # Create the workflow. Note that any keywords provided in the calculator_parameters (i.e. whatever is left in
        # calcdict) are provided as kwargs
        return cls(atoms, parameters=parameters, kpoints=kpts, calculator_parameters=calculator_parameters, **kwargs,
                   **calcdict)

    def print_header(self):
        print(header())

    def print_bib(self):
        relevant_references = BibliographyData()

        def add_ref(bibkey: str, note: str):
            if bibkey not in bib_data.entries:
                raise ValueError(f'Could not find bibliography entry for {bibkey}')
            else:
                entry = bib_data.entries[bibkey]
                entry.fields['note'] = note
                relevant_references.add_entry(bibkey, entry)

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

        print(f'\n Please cite the papers listed in {self.name}.bib in work involving this calculation')
        relevant_references.to_file(self.name + '.bib')

    def print_preamble(self):
        if self.parent:
            return

        self.print_header()

        self.print_bib()

    def print_conclusion(self):
        from koopmans.io import write

        if self.parent:
            return

        # Save workflow to file
        write(self, self.name + '.kwf')

        # Print farewell message
        print('\n Workflow complete')

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
                nested_keys = code.split('_')[1:]
                # The following very opaque code fills out the nested dictionary with the list of nested keys
                for i, k in enumerate(nested_keys):
                    parent_level = reduce(operator.getitem, nested_keys[:i], calcdct['w90'])
                    if k not in parent_level:
                        reduce(operator.getitem, nested_keys[:i], calcdct['w90'])[k] = {}
                reduce(operator.getitem, nested_keys[:-1], calcdct['w90'])[k] = params_dict
                # Projections
                filling = nested_keys[0] == 'occ'
                if len(nested_keys) == 2:
                    spin = nested_keys[1]
                else:
                    spin = None
                projections = self.projections.get_subset(filling, spin)
                if len(projections) > 1:
                    proj_kwarg = {'projections_blocks': [p.projections for p in projections]}
                elif len(projections) == 1:
                    proj_kwarg = {'projections': projections[0].projections}
                else:
                    proj_kwarg = {}
                reduce(operator.getitem, nested_keys[:-1], calcdct['w90'])[k].update(**proj_kwarg)
            else:
                raise NotImplementedError(
                    f'Writing of {params.__class__.__name__} with write_json is not yet implemented')

        other_blocks: Dict[str, Any] = {'plotting': self.plotting}
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
            raise ValueError('The "bs" and "bsplot_kwargs" arguments to plot_bandstructure() should be the same length')
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
        workflow_name = self.__class__.__name__.lower()
        for s in ['workflow', 'mock', 'benchgen', 'stumbling', 'check']:
            workflow_name = workflow_name.replace(s, '')
        filename = filename if filename is not None else f'{self.name}_{workflow_name}_bandstructure'
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


def header():
    from koopmans import __version__

    header = [r"  _",
              r" | | _____   ___  _ __  _ __ ___   __ _ _ __  ___",
              r" | |/ / _ \ / _ \| '_ \| '_ ` _ \ / _` | '_ \/ __|",
              r" |   < (_) | (_) | |_) | | | | | | (_| | | | \__ \ ",
              r" |_|\_\___/ \___/| .__/|_| |_| |_|\__,_|_| |_|___/",
              f"                 |_|",
              "",
              " Koopmans spectral functional calculations with Quantum ESPRESSO",
              "",
              f" version {__version__}",
              "",
              " Written by Edward Linscott, Riccardo De Gennaro, and Nicola Colonna"]
    return '\n'.join(header)


def read_atoms_dict(dct: Dict[str, Any]):
    '''
    Reads the "atoms" block
    '''

    atoms = Atoms()

    readers: Dict[str, Callable] = {'cell_parameters': utils.read_cell_parameters,
                                    'atomic_positions': utils.read_atomic_positions}

    for key, reader in readers.items():
        subdct: Dict[str, Any] = dct.pop(key, {})
        if subdct:
            reader(atoms, subdct)
        else:
            raise ValueError(f'Please provide "{key}" in the atoms block')

    for block in dct:
        raise ValueError(f'Unrecognized subblock atoms: "{block}"')

    return atoms


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
            'w90_occ': settings.Wannier90SettingsDict(),
            'w90_emp': settings.Wannier90SettingsDict(),
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
                    'w90_occ': settings.Wannier90SettingsDict,
                    'w90_emp': settings.Wannier90SettingsDict,
                    'w90_occ_up': settings.Wannier90SettingsDict,
                    'w90_emp_up': settings.Wannier90SettingsDict,
                    'w90_occ_down': settings.Wannier90SettingsDict,
                    'w90_emp_down': settings.Wannier90SettingsDict,
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
                f'Unrecognized calculator_parameters entry "{k}": valid options are '
                '/'.join(settings_classes.keys()))
    return dct_out
