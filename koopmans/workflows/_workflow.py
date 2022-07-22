"""

Generic workflow object for koopmans

Written by Edward Linscott Oct 2020
Converted workflows from functions to objects Nov 2020

"""

import copy
import json as json_ext
import operator
import os
import shutil
import subprocess
from abc import ABC, abstractmethod
from functools import reduce
from pathlib import Path
from socket import IP_DROP_MEMBERSHIP
from types import ModuleType
from typing import Any, Dict, List, Optional, Type, TypeVar, Union

import matplotlib.pyplot as plt
import numpy as np
from numpy import typing as npt
from pybtex.database import BibliographyData

import ase
import koopmans.mpl_config
from ase import Atoms
from ase.build.supercells import make_supercell
from ase.calculators.calculator import CalculationFailed
from ase.calculators.espresso import Espresso_kcp
from ase.dft.dos import DOS
from ase.dft.kpoints import BandPath
from ase.io.espresso import cell_to_ibrav
from ase.io.espresso import contruct_kcp_namelist as construct_namelist
from ase.io.espresso import ibrav_to_cell, kcp_keys
from ase.spectrum.band_structure import BandStructure
from ase.spectrum.doscollection import GridDOSCollection
from ase.spectrum.dosdata import GridDOSData
from koopmans import calculators, settings, utils
from koopmans.bands import Bands
from koopmans.commands import ParallelCommandWithPostfix
from koopmans.ml_utils._ml_models import RidgeRegression
from koopmans.projections import ProjectionBlocks
from koopmans.pseudopotentials import (fetch_pseudo, nelec_from_pseudos,
                                       pseudo_database,
                                       pseudos_library_directory,
                                       valence_from_pseudo)
from koopmans.references import bib_data

T = TypeVar('T', bound='calculators.CalculatorExt')


class Workflow(ABC):

    def __init__(self, atoms: Atoms,
                 parameters: settings.SettingsDict = settings.WorkflowSettingsDict(),
                 master_calc_params: Optional[Union[Dict[str, Dict[str, Any]],
                                                    Dict[str, settings.SettingsDict]]] = None,
                 name: str = 'koopmans_workflow',
                 pseudopotentials: Dict[str, str] = {},
                 pseudo_dir: Optional[Path] = None,
                 gamma_only: Optional[bool] = False,
                 kgrid: Optional[List[int]] = [1, 1, 1],
                 koffset: Optional[List[int]] = [0, 0, 0],
                 kpath: Optional[Union[BandPath, str]] = None,
                 kpath_density: int = 10,
                 projections: Optional[ProjectionBlocks] = None,
                 plot_params: Union[Dict[str, Any], settings.PlotSettingsDict] = {},
                 autogenerate_settings: bool = True):

        # Parsing parameters
        self.parameters: Union[settings.MLSettingsDict, settings.WorkflowSettingsDict]
        if hasattr(parameters, 'use_ml') and parameters.use_ml:
            self.parameters = settings.MLSettingsDict(**parameters)
        else:
            self.parameters = settings.WorkflowSettingsDict(**parameters)

        self.atoms = atoms
        self.name = name
        self.calculations: List[calculators.Calc] = []
        self.silent = False
        self.print_indent = 1
        self.gamma_only = gamma_only
        if self.gamma_only:
            self.kgrid = None
            self.koffset = [0, 0, 0]
        else:
            self.kgrid = kgrid
            self.koffset = koffset

        if projections is None:
            proj_list: List[List[Any]]
            spins: List[Optional[str]]
            if self.parameters.spin_polarised:
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

        self.plot_params = settings.PlotSettingsDict(**plot_params)

        if 'periodic' in parameters:
            # If "periodic" was explicitly provided, override self.atoms.pbc
            self.atoms.pbc = self.parameters.periodic
        else:
            # If "periodic" was not explicitly provided, use the value from self.atoms.pbc
            self.parameters.periodic = all(self.atoms.pbc)

        if all(self.atoms.pbc):
            self.atoms.wrap(pbc=True)

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
                    utils.warn(f'You are using an unrecognised pseudopotential {pseudo.name}. Please note that '
                               'the current implementation of Koopmans functionals only supports norm-conserving '
                               'pseudopotentials.')
                elif pseudo.kind != 'norm-conserving':
                    raise ValueError('Koopmans functionals only currently support norm-conserving pseudopotentials; '
                                     f'{pseudo.name} is {pseudo.kind}')
                if tag > 0:
                    symbol += str(tag)
                self.pseudopotentials[symbol] = pseudo.name

        # Make sure master_calc_params isn't missing any entries, and every entry corresponds to settings.SettingsDict
        # objects
        master_calc_params = sanitise_master_calc_params(
            master_calc_params) if master_calc_params is not None else generate_default_master_calc_params()

        # Work out the pseudopotential directory. If using a pseudo_library this is straightforward, if not...
        #  1. try to locating the directory as currently specified by the calculator
        #  2. if that fails, check if $ESPRESSO_PSEUDO is set
        #  3. if that fails, raise an error
        if pseudo_dir is not None:
            pass
        elif self.parameters.pseudo_library:
            pseudo_dir = pseudos_library_directory(self.parameters.pseudo_library, self.parameters.base_functional)
            for params in master_calc_params.values():
                if params.get('pseudo_dir', pseudo_dir).resolve() != pseudo_dir:
                    raise ValueError(
                        '"pseudo_dir" and "pseudo_library" are conflicting; please do not provide "pseudo_dir"')
        elif 'pseudo_dir' in master_calc_params['kcp'] or 'pseudo_dir' in master_calc_params['pw']:
            pseudo_dir = master_calc_params['kcp'].get('pseudo_dir', master_calc_params['pw'].get('pseudo_dir'))
            assert isinstance(pseudo_dir, Path)
        elif 'ESPRESSO_PSEUDO' in os.environ:
            pseudo_dir = Path(os.environ['ESPRESSO_PSEUDO'])
        else:
            pseudo_dir = Path.cwd()

        self.pseudo_dir = pseudo_dir

        # Before saving the master_calc_params, automatically generate some keywords and perform some sanity checks
        if self.parameters.task != 'ui' and autogenerate_settings:
            # Automatically calculate nelec/nelup/neldw/etc using information contained in the pseudopotential files
            # and the kcp settings
            nelec = nelec_from_pseudos(self.atoms, self.pseudopotentials, self.pseudo_dir)

            tot_charge = master_calc_params['kcp'].get('tot_charge', 0)
            nelec -= tot_charge
            tot_mag = master_calc_params['kcp'].get('tot_magnetization', nelec % 2)
            nelup = int(nelec / 2 + tot_mag / 2)
            neldw = int(nelec / 2 - tot_mag / 2)

            # Setting up the magnetic moments
            if 'starting_magnetization(1)' in master_calc_params['kcp']:
                labels = [s + str(t) if t > 0 else s for s, t in zip(atoms.symbols, atoms.get_tags())]
                starting_magmoms = {}
                for i, (l, p) in enumerate(self.pseudopotentials.items()):
                    # ASE uses absoulte values; QE uses the fraction of the valence
                    frac_mag = master_calc_params['kcp'].pop(f'starting_magnetization({i + 1})', 0.0)
                    valence = valence_from_pseudo(p, self.pseudo_dir)
                    starting_magmoms[l] = frac_mag * valence
                atoms.set_initial_magnetic_moments([starting_magmoms[l] for l in labels])
            elif tot_mag != 0:
                atoms.set_initial_magnetic_moments([tot_mag / len(atoms) for _ in atoms])

            # Work out the number of bands
            nbnd = master_calc_params['kcp'].get('nbnd', nelec // 2 + nelec % 2)
            generated_keywords = {'nelec': nelec, 'tot_charge': tot_charge, 'tot_magnetization': tot_mag,
                                  'nelup': nelup, 'neldw': neldw, 'nbnd': nbnd, 'pseudo_dir': self.pseudo_dir}
        else:
            generated_keywords = {}
            nelec = 0

        self.master_calc_params = generate_default_master_calc_params()
        for block, params in master_calc_params.items():
            # Apply auto-generated keywords
            for k, v in generated_keywords.items():
                # Skipping nbnd for kcp -- it is valid according to ASE but it is not yet properly implemented
                if k == 'nbnd' and block == 'kcp':
                    continue
                if k in params.valid and k not in params:
                    setattr(params, k, v)

            # Various checks for the wannier90 blocks
            if block.startswith('w90'):
                # If we are spin-polarised, don't store the spin-independent w90 block
                # Likewise, if we are not spin-polarised, don't store the spin-dependent w90 blocks
                if self.parameters.spin_polarised is not ('up' in block or 'down' in block):
                    continue
                if 'projections' in params or 'projections_blocks' in params:
                    raise ValueError(f'You have provided projection information in the master_calc_params[{block}] '
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

            # Store the sanitised parameters
            self.master_calc_params[block] = params

        # Generate a default kpath
        if kpath is None:
            if self.parameters.periodic:
                # By default, use ASE's default bandpath for this cell (see
                # https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#brillouin-zone-data)
                kpath = self.atoms.cell.bandpath().path
            else:
                kpath = 'G'

        # Convert the kpath to a BandPath object
        if isinstance(kpath, str):
            self.kpath = utils.convert_kpath_str_to_bandpath(kpath, self.atoms.cell, kpath_density)
        else:
            self.kpath = kpath

        # If atoms has a calculator, overwrite the kpoints and pseudopotentials variables and then detach the calculator
        if atoms.calc is not None:
            utils.warn(f'You have initialised a {self.__class__.__name__} object with an atoms object that possesses '
                       'a calculator. This calculator will be ignored.')
            self.atoms.calc = None

        # Records whether or not this workflow is a subworkflow of another
        self._is_a_subworkflow = False

        # Initialize the RidgeRegression() model
        if self.parameters.use_ml:
            self.ml_model = RidgeRegression()

    def __eq__(self, other: Any):
        if isinstance(other, Workflow):
            return self.__dict__ == other.__dict__
        return False

    def run(self) -> None:
        self.print_preamble()
        if not self._is_a_subworkflow:
            self._run_sanity_checks()
        self._run()
        self.print_conclusion()
        if not self._is_a_subworkflow:
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

    @property
    def gamma_only(self) -> Optional[bool]:
        return self._gamma_only

    @gamma_only.setter
    def gamma_only(self, value: Optional[bool]):
        self._gamma_only = value

    @property
    def kgrid(self):
        return self._kgrid

    @kgrid.setter
    def kgrid(self, value: List[int]):
        self._kgrid = value

    @property
    def koffset(self):
        return self._koffset

    @koffset.setter
    def koffset(self, value: List[int]):
        self._koffset = value

    @property
    def kpath(self):
        return self._kpath

    @kpath.setter
    def kpath(self, value: Union[str, BandPath]):
        if isinstance(value, str):
            raise NotImplementedError()
        self._kpath = value

    @property
    def wf_kwargs(self) -> Dict[str, Any]:
        # Returns a kwargs designed to be used to initialise another workflow with the same configuration as this one
        # i.e.
        # > sub_wf = Workflow(**self.wf_kwargs)
        return {'atoms': copy.deepcopy(self.atoms),
                'parameters': copy.deepcopy(self.parameters),
                'master_calc_params': copy.deepcopy(self.master_calc_params),
                'name': copy.deepcopy(self.name),
                'pseudopotentials': copy.deepcopy(self.pseudopotentials),
                'pseudo_dir': copy.deepcopy(self.pseudo_dir),
                'gamma_only': copy.deepcopy(self.gamma_only),
                'kgrid': copy.deepcopy(self.kgrid),
                'kpath': copy.deepcopy(self.kpath),
                'projections': copy.deepcopy(self.projections),
                'plot_params': copy.deepcopy(self.plot_params)}

    def _run_sanity_checks(self):
        # Check internal consistency of workflow settings
        if self.parameters.fix_spin_contamination is None:
            self.parameters.fix_spin_contamination = not self.parameters.spin_polarised
        else:
            if self.parameters.fix_spin_contamination and self.parameters.spin_polarised:
                raise ValueError('fix_spin_contamination = True is incompatible with spin_polarised = True')

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

        if self.parameters.periodic:
            if self.parameters.gb_correction is None:
                self.parameters.gb_correction = True

            if self.parameters.mp_correction:
                if self.parameters.eps_inf is None:
                    raise ValueError('eps_inf missing in input; needed when mp_correction is true')
                elif self.parameters.eps_inf < 1.0:
                    raise ValueError('eps_inf cannot be lower than 1')
            else:
                utils.warn('Makov-Payne corrections not applied for a periodic calculation; do this with '
                           'caution')

            if self.parameters.mt_correction is None:
                self.parameters.mt_correction = False
            if self.parameters.mt_correction:
                raise ValueError('Do not use Martyna-Tuckerman corrections for periodic systems')

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

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs']:
            if len(self.projections) == 0:
                raise ValueError(f'In order to use init_orbitals={self.parameters.init_orbitals}, projections must be '
                                 'provided')
            spin_set = set([p.spin for p in self.projections])
            if self.parameters.spin_polarised:
                if spin_set != {'up', 'down'}:
                    raise ValueError('This calculation is spin-polarised; please provide spin-up and spin-down '
                                     'projections')
            else:
                if spin_set != {None}:
                    raise ValueError('This calculation is not spin-polarised; please do not provide spin-indexed '
                                     'projections')

        # Check the consistency between self.gamma_only and KCP's do_wf_cmplx
        if not self.gamma_only and self.master_calc_params['kcp'].do_wf_cmplx is False:
            utils.warn('In KCP do_wf_cmplx = False is not consistent with gamma_only = False. '
                       'Changing do_wf_cmplx to True')
            self.master_calc_params['kcp'].do_wf_cmplx = True

        # Check pseudopotentials exist
        if not os.path.isdir(self.pseudo_dir):
            raise NotADirectoryError(f'The pseudo_dir you provided ({self.pseudo_dir}) does not exist')
        if self.parameters.task != 'ui':
            for pseudo in self.pseudopotentials.values():
                if not (self.pseudo_dir / pseudo).exists():
                    raise FileNotFoundError(
                        f'{self.pseudo_dir / pseudo} does not exist. Please double-check your pseudopotential settings')

        # Make sanity checks for the ML model
        if self.parameters.use_ml:
            if self.parameters.task != 'trajectory':
                raise NotImplementedError(
                    f'Using the ML-prediction for the {self.parameter.task}-task has not yet been implemented.')
            if self.parameters.method != 'dscf':  # TODO Yannick: implement also with DFPT
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
            if self.parameters.spin_polarised:
                utils.warn(f'Using the ML-prediction for spin-polarised systems has not yet been implemented.')
            if not self.parameters.periodic:
                utils.warn(f'Using the ML-prediction for non-periodic systems has not yet been extensively tested.')
            if self.parameters.orbital_groups:
                utils.warn('Using orbital_groups has not yet been extensively tested.')
            if not np.all(self.atoms.cell.angles() == 90.0):
                raise ValueError(f"The ML-workflow has only been implemented for simulation cells that have 90° angles")

    def new_calculator(self,
                       calc_type: str,
                       directory: Optional[Path] = None,
                       kpts: Optional[Union[List[int], BandPath]] = None,
                       **kwargs) -> T:

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
        elif calc_type == 'projwfc':
            calc_class = calculators.ProjwfcCalculator
        else:
            raise ValueError(f'Cound not find a calculator of type {calc_type}')

        # Merge master_calc_params and kwargs, giving kwargs higher precedence
        all_kwargs: Dict[str, Any] = {}
        master_calc_params = self.master_calc_params[calc_type]
        all_kwargs.update(**master_calc_params)
        all_kwargs.update(**kwargs)

        # For the k-points, the Workflow has two options: self.kgrid and self.kpath. A calculator should only ever
        # have one of these two. By default, use the kgrid.
        if 'kpts' in master_calc_params.valid:
            all_kwargs['kpts'] = kpts if kpts is not None else self.kgrid

        # Add pseudopotential and kpt information to the calculator as required
        for kw in ['pseudopotentials', 'pseudo_dir', 'gamma_only', 'kgrid', 'kpath', 'koffset', 'plot_params']:
            if kw not in all_kwargs and kw in master_calc_params.valid:
                all_kwargs[kw] = getattr(self, kw)

        # Create the calculator
        calc = calc_class(atoms=copy.deepcopy(self.atoms), **all_kwargs)

        # Add the directory if provided
        if directory is not None:
            calc.directory = directory

        return calc

    def update_celldms(self):
        # Update celldm(*) to match the current self.atoms.cell
        for k, params in self.master_calc_params.items():
            if params.get('ibrav', 0) != 0:
                celldms = cell_to_ibrav(self.atoms.cell, params.ibrav)
                self.master_calc_params[k].update(**celldms)

    def primitive_to_supercell(self, matrix: Optional[npt.NDArray[np.int_]] = None, **kwargs):
        # Converts to a supercell as given by a 3x3 transformation matrix
        if matrix is None:
            matrix = np.diag(self.kgrid) if not self.gamma_only else np.identity(3, dtype=float)
        assert np.shape(matrix) == (3, 3)
        self.atoms = make_supercell(self.atoms, matrix, **kwargs)

        self.update_celldms()

    def supercell_to_primitive(self, matrix: Optional[npt.NDArray[np.int_]] = None):
        # Converts from a supercell to a primitive cell, as given by a 3x3 transformation matrix
        # The inverse of self.primitive_to_supercell()
        if matrix is None:
            matrix = np.diag(self.kgrid)
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

        self.update_celldms()

    def run_calculator(self, master_qe_calc, enforce_ss=False):
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

    def run_calculator_single(self, qe_calc):
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
        except CalculationFailed as e:
            self.print(' failed')
            raise CalculationFailed(e)

        if not self.silent:
            self.print(' done')

        # Store the calculator
        self.calculations.append(qe_calc)

        # If we reached here, all future calculations should be performed from scratch
        self.parameters.from_scratch = True

        return

    def load_old_calculator(self, qe_calc):
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

    def print(self, text='', style='body', **kwargs):
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

    def run_subworkflow(self, workflow, subdirectory=None, from_scratch=None, **kwargs):
        '''
        Runs a workflow object, taking care of inheritance of several important properties

        '''

        # When testing, make sure the sub-workflow has access to the benchmark
        if hasattr(self, 'benchmark'):
            workflow.benchmark = self.benchmark

        # Automatically pass along the name of the overall workflow
        if workflow.name is None:
            workflow.name = self.name

        # Increase the indent level
        workflow.print_indent = self.print_indent + 1

        # Don't print out the header, generate .bib and .kwf files etc. for subworkflows
        workflow._is_a_subworkflow = True

        # Ensure altering workflow.master_calc_params won't affect self.master_calc_params
        if workflow.master_calc_params is self.master_calc_params:
            workflow.master_calc_params = copy.deepcopy(self.master_calc_params)

        # Setting from_scratch to a non-None value will override the value of subworkflow.from_scratch...
        if from_scratch is None:
            workflow.parameters.from_scratch = self.parameters.from_scratch
        else:
            workflow.parameters.from_scratch = from_scratch

        # Link the list of calculations
        workflow.calculations = self.calculations

        # Link the ML_Model
        if self.parameters.use_ml:
            workflow.ml_model = self.ml_model

        # Link the bands
        if hasattr(self, 'bands'):
            workflow.bands = copy.deepcopy(self.bands)
            # Only include the most recent screening parameter and wipe the error history
            for b in workflow.bands:
                if len(b.alpha_history) > 0:
                    b.alpha_history = [b.alpha]
                b.error_history = []

        if subdirectory is not None:
            # Update directories
            for key in workflow.master_calc_params.keys():
                params = workflow.master_calc_params[key]
                for setting in params.are_paths:
                    if setting == 'pseudo_dir':
                        continue
                    path = getattr(params, setting, None)
                    if path is not None and Path.cwd() in path.parents:
                        new_path = Path(subdirectory).resolve() / os.path.relpath(path)
                        setattr(params, setting, new_path)

            # Run the workflow
            with utils.chdir(subdirectory):
                workflow.run(**kwargs)
        else:
            workflow.run(**kwargs)

        # ... and will prevent inheritance of from_scratch
        if from_scratch is None:
            self.parameters.from_scratch = workflow.parameters.from_scratch

        # Copy back over the bands
        if hasattr(workflow, 'bands'):
            if hasattr(self, 'bands'):
                # Add the alpha and error history
                for b, b_sub in zip(self.bands, workflow.bands):
                    b.alpha_history += b_sub.alpha_history[1:]
                    b.error_history += b_sub.error_history
            else:
                # Copy the entire bands object
                self.bands = workflow.bands

        # Make sure any updates to the projections are passed along
        self.projections = workflow.projections

    def todict(self):
        # Shallow copy
        dct = dict(self.__dict__)

        # Adding information required by the json decoder
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @classmethod
    def fromdict(cls, dct: Dict):
        wf = cls(atoms=dct.pop('atoms'),
                 parameters=dct.pop('parameters'),
                 master_calc_params=dct.pop('master_calc_params'),
                 pseudopotentials=dct.pop('_pseudopotentials'),
                 gamma_only=dct.pop('_gamma_only'),
                 kgrid=dct.pop('_kgrid'),
                 kpath=dct.pop('_kpath'),
                 projections=dct.pop('projections'),
                 autogenerate_settings=False)

        for k, v in dct.items():
            setattr(wf, k, v)

        return wf

    @property
    def bands(self):
        if not hasattr(self, '_bands'):
            raise AttributeError('Bands have not been initialised')
        return self._bands

    @bands.setter
    def bands(self, value):
        assert isinstance(value, Bands)
        self._bands = value

    @classmethod
    def fromjson(cls, fname: str):

        with open(fname, 'r') as fd:
            bigdct = json_ext.loads(fd.read())
        wf = cls._fromjsondct(bigdct)
        wf.name = fname.replace('.json', '')
        return wf

    @classmethod
    def _fromjsondct(cls, bigdct: Dict[str, Any]):

        # Deal with the nested w90 subdictionaries
        if 'w90' in bigdct:
            for filling in ['occ', 'emp']:
                for spin in ['up', 'down']:
                    # Add any keywords in the filling:spin subsubdictionary
                    subsubdct = bigdct['w90'].get(filling, {}).get(spin, {})
                    bigdct[f'w90_{filling}_{spin}'] = subsubdct
                    # Add any keywords in the filling subdictionary
                    subdct = {k: v for k, v in bigdct['w90'].get(filling, {}).items() if k not in ['up', 'down']}
                    bigdct[f'w90_{filling}_{spin}'].update(subdct)
                    # Add any keywords in the main dictionary
                    dct = {k: v for k, v in bigdct['w90'].items() if k not in ['occ', 'emp']}
                    bigdct[f'w90_{filling}_{spin}'].update(dct)
                # Also create a spin-independent set of parameters
                bigdct[f'w90_{filling}'] = {}
                bigdct[f'w90_{filling}'].update(subdct)
                bigdct[f'w90_{filling}'].update(dct)
            # Finally, remove the nested w90 entry
            del bigdct['w90']

        # Deal with UI subdicts
        if 'ui' in bigdct:
            subdcts = {}
            keys = ['occ', 'emp']
            for key in keys:
                # First, we must remove the occ and emp subdicts from the UI dict
                if key in bigdct['ui']:
                    subdcts[key] = bigdct['ui'].pop(key)

            # Now, we add the ui_occ and ui_emp calculators to master_calc_params
            for key in keys:
                if key in subdcts:
                    # Add the corresponding subdict to the rest of the UI block
                    bigdct[f'ui_{key}'] = dict(bigdct['ui'], **subdcts[key])
                else:
                    # Duplicate the UI block
                    bigdct[f'ui_{key}'] = bigdct['ui']

        # Deal with kc_wann subdicts
        kc_wann_blocks = bigdct.pop('kc_wann', {'kc_ham': {}, 'kc_screen': {}, 'wann2kc': {}})
        bigdct.update(**kc_wann_blocks)

        # Loading plot settings
        plot_params = settings.PlotSettingsDict(**utils.parse_dict(bigdct.get('plot', {})))

        # Loading workflow settings
        parameters: Union[settings.MLSettingsDict, settings.WorkflowSettingsDict]
        if 'ML' in bigdct:
            if bigdct['ML']['use_ml']:  # If the user wants to use the ML model, a bigger dictionary needs to be loaded to the parameters
                parameters = settings.MLSettingsDict(**utils.parse_dict(bigdct.get('workflow', {})),
                                                     **utils.parse_dict(bigdct['ML']))
            else:
                parameters = settings.WorkflowSettingsDict(**utils.parse_dict(bigdct.get('workflow', {})))
            bigdct.pop('ML')  # remove ML from the master-calc params
        else:
            parameters = settings.WorkflowSettingsDict(**utils.parse_dict(bigdct.get('workflow', {})))

        # Check for unexpected blocks
        for block in bigdct:
            if block not in list(settings_classes.keys()) + ['workflow', 'setup']:
                raise ValueError(f'Unrecognised block "{block}" in json input file; '
                                 'valid options are workflow/' + '/'.join(settings_classes.keys()))

        # Load default values
        if 'setup' in bigdct:
            atoms, setup_parameters, workflow_kwargs = read_setup_dict(bigdct['setup'], parameters.task)
            del bigdct['setup']
        elif parameters.task != 'ui':
            raise ValueError('You must provide a "setup" block in the input file, specifying atomic positions, atomic '
                             'species, etc.')
        else:
            # Create dummy objects
            atoms = Atoms()
            setup_parameters = {}
            workflow_kwargs = {}

        # Loading calculator-specific settings. We generate a SettingsDict for every single kind of calculator,
        # regardless of whether or not there was a corresponding block in the json file
        master_calc_params = {}
        w90_block_projs: List = []
        w90_block_filling: List[bool] = []
        w90_block_spins: List[Union[str, None]] = []
        for block, settings_class in settings_classes.items():
            # Read the block and add the resulting calculator to the calcs_dct
            dct = bigdct.get(block, {})
            if block.startswith('ui'):
                # Dealing with redundancies in UI keywords
                if 'sc_dim' in dct and 'kpts' in workflow_kwargs:
                    # In this case, the sc_dim keyword is redundant
                    if workflow_kwargs['kpts'] != dct['sc_dim']:
                        raise ValueError('sc_dim in the UI block should match the kpoints provided in the setup block')
                    dct.pop('sc_dim')
                if 'kpath' in dct and 'kpath' in workflow_kwargs:
                    if workflow_kwargs['kpath'] != dct['kpath']:
                        raise ValueError('kpath in the UI block should match that provided in the setup block')
                    dct.pop('kpath')
            elif block.startswith('w90'):
                # If we are spin-polarised, don't store the spin-independent w90 block
                # Likewise, if we are not spin-polarised, don't store the spin-dependent w90 blocks
                if parameters.spin_polarised is not ('up' in block or 'down' in block):
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

            master_calc_params[block] = settings_class(**dct)
            master_calc_params[block].update(
                **{k: v for k, v in setup_parameters.items() if master_calc_params[block].is_valid(k)})

        # Adding the projections to the workflow kwargs (this is unusual in that this is an attribute of the workflow
        # object but it is provided in the w90 subdictionary)
        workflow_kwargs['projections'] = ProjectionBlocks.fromprojections(
            w90_block_projs, w90_block_filling, w90_block_spins, atoms)

        workflow_kwargs['plot_params'] = plot_params

        return cls(atoms, parameters, master_calc_params, **workflow_kwargs)

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

            if self.parameters.periodic:
                add_ref('Nguyen2018', 'Describes Koopmans functionals in periodic systems')
                if self.parameters.calculate_alpha:
                    if self.parameters.method == 'dfpt':
                        add_ref('Colonna2019', 'Introduces the DFPT method for calculating screening parameters')
                        add_ref('Colonna2022', 'Describes the algorithms underpinning the kcw.x code')
                    else:
                        add_ref('DeGennaro2021', 'Describes how to extract band structures from Koopmans functional '
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
        if self._is_a_subworkflow:
            return

        self.print_header()

        self.print_bib()

    def print_conclusion(self):
        from koopmans.io import write

        if self._is_a_subworkflow:
            return

        # Save workflow to file
        write(self, self.name + '.kwf')

        # Print farewell message
        print('\n Workflow complete')

    def toinputjson(self) -> Dict[str, Dict]:

        bigdct: Dict[str, Dict] = {}

        bigdct['workflow'] = {}

        # "workflow" block (not printing any values that match the defaults except for core keywords)
        for k, v in self.parameters.items():
            if v is None:
                continue
            if isinstance(v, Path):
                v = str(v)
            default = self.parameters.defaults.get(k, None)
            if v != default or k in ['task', 'functional']:
                bigdct['workflow'][k] = v

        # "setup" block
        # Working out ibrav
        ibrav = self.master_calc_params['kcp'].get('ibrav', self.master_calc_params['pw'].get('ibrav', 0))

        bigdct['setup'] = {}

        # cell parameters
        if ibrav == 0:
            bigdct['setup']['cell_parameters'] = utils.construct_cell_parameters_block(self.atoms)

        # atomic positions
        bigdct['setup']['atomic_positions'] = utils.construct_atomic_positions_block(self.atoms, ibrav != 0)

        # k-points
        bigdct['setup']['k_points'] = {'kgrid': self.kgrid, 'kpath': self.kpath.path}

        # Populating calculator-specific blocks
        bigdct['w90'] = {}
        bigdct['ui'] = {}
        for code, params in self.master_calc_params.items():
            # Remove default settings (ensuring we switch to using relative paths to check this)
            tmp, params.use_relative_paths = params.use_relative_paths, True
            params_dict = {k: v for k, v in params.items() if params.defaults.get(k, None) != v}
            params.use_relative_paths = tmp

            # convert Paths to strings
            for k in params_dict:
                if isinstance(params_dict[k], Path):
                    params_dict[k] = str(params_dict[k])

            # pseudo directory belongs in setup, not elsewhere
            pseudo_dir = params_dict.pop('pseudo_dir', None)
            if pseudo_dir is not None and self.parameters.pseudo_library is None:
                bigdct['setup']['control'] = {'pseudo_dir': str(pseudo_dir)}

            # If the params_dict is empty, don't add a block for this calculator
            if not params_dict and not code.startswith('w90'):
                continue

            if code in ['pw', 'kcp']:
                bigdct[code] = {}

                # Populate bigdct with the settings
                input_data = construct_namelist(params_dict)
                for key, block in input_data.items():

                    if len(block) > 0:
                        bigdct[code][key] = {k: v for k, v in dict(
                            block).items() if v is not None}

            elif code in ['pw2wannier', 'wann2kc', 'kc_screen', 'kc_ham', 'projwfc', 'wann2kcp', 'plot']:
                bigdct[code] = params_dict
            elif code.startswith('ui_'):
                bigdct['ui'][code.split('_')[-1]] = params_dict
            elif code == 'ui':
                bigdct['ui'].update(**params_dict)
            elif code.startswith('w90'):
                nested_keys = code.split('_')[1:]
                # The following very opaque code fills out the nested dictionary with the list of nested keys
                for i, k in enumerate(nested_keys):
                    parent_level = reduce(operator.getitem, nested_keys[:i], bigdct['w90'])
                    if k not in parent_level:
                        reduce(operator.getitem, nested_keys[:i], bigdct['w90'])[k] = {}
                reduce(operator.getitem, nested_keys[:-1], bigdct['w90'])[k] = params_dict
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
                reduce(operator.getitem, nested_keys[:-1], bigdct['w90'])[k].update(**proj_kwarg)
            else:
                raise NotImplementedError(
                    f'Writing of {params.__class__.__name__} with write_json is not yet implemented')

        return bigdct

    def plot_bandstructure(self,
                           bs: Union[BandStructure, List[BandStructure]],
                           dos: Optional[Union[GridDOSCollection, DOS]] = None,
                           filename: Optional[str] = None,
                           bsplot_kwargs: Union[Dict[str, Any], List[Dict[str, Any]]] = {},
                           dosplot_kwargs: Dict[str, Any] = {}):
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

        # Sanitise input
        if isinstance(bs, BandStructure):
            bs = [bs]
        if isinstance(bsplot_kwargs, dict):
            bsplot_kwargs = [bsplot_kwargs]
        if len(bs) != len(bsplot_kwargs):
            raise ValueError('The "bs" and "bsplot_kwargs" arguments to plot_bandstructure() should be the same length')
        spins: List[Optional[str]]
        if isinstance(dos, DOS):
            if self.parameters.spin_polarised:
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
        defaults = {'colors': colors, 'emin': self.plot_params.Emin, 'emax': self.plot_params.Emax}
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
            if self.parameters.spin_polarised:
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

                label: Union[str, None]
                for d in sorted_dos:
                    if (not self.parameters.spin_polarised or d.info.get('spin') == 'up') \
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
            if self.parameters.spin_polarised:
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
        filename = filename if filename is not None else f'{self.name}_bandstructure'
        legends = [ax.get_legend() for ax in axes if ax.get_legend() is not None]
        utils.savefig(fname=filename + '.png', bbox_extra_artists=legends, bbox_inches='tight')

    def _teardown(self):
        '''
        Performs final tasks before the workflow completes
        '''

        # Removing tmpdirs
        if not self.parameters.keep_tmpdirs:
            all_outdirs = [calc.parameters.get('outdir', None) for calc in self.calculations]
            outdirs = set([o.resolve() for o in all_outdirs if o is not None and o.resolve().exists()])
            for outdir in outdirs:
                shutil.rmtree(outdir)


def get_version(module):
    if isinstance(module, ModuleType):
        module = module.__path__[0]
    with utils.chdir(module):
        version_label = subprocess.check_output(["git", "describe", "--always", "--tags"]).strip()
    return version_label.decode("utf-8")


def header():

    koopmans_version = get_version(os.path.dirname(__file__))
    ase_version = get_version(ase)
    qe_version = get_version((calculators.bin_directory / 'pw.x').resolve().parents[2])

    header = [r"  _                                                ",
              r" | | _____   ___  _ __  _ __ ___   __ _ _ __  ___  ",
              r" | |/ / _ \ / _ \| '_ \| '_ ` _ \ / _` | '_ \/ __| ",
              r" |   < (_) | (_) | |_) | | | | | | (_| | | | \__ \ ",
              r" |_|\_\___/ \___/| .__/|_| |_| |_|\__,_|_| |_|___/ ",
              r"                 |_|                               ",
              "",
              " Koopmans spectral functional calculations with Quantum ESPRESSO",
              "",
              " Written by Edward Linscott, Riccardo De Gennaro, and Nicola Colonna",
              "",
              f" using QE version {qe_version}, workflow manager version {koopmans_version}, and ASE version "
              f"{ase_version}"
              ""]
    return '\n'.join(header)


def read_setup_dict(dct: Dict[str, Any], task: str):
    '''

    Reads the "setup" block. This block uses the same syntax as kcp

    '''

    calc = Espresso_kcp(atoms=Atoms())

    compulsory_block_readers = {'atomic_positions': utils.read_atomic_positions}

    for block, subdct in dct.items():
        if block in compulsory_block_readers or block in ['cell_parameters', 'k_points', 'atomic_species']:
            # We will read these afterwards
            continue
        elif block in kcp_keys:
            for key, value in subdct.items():
                if value == "":
                    continue

                # Force pseudo_dir to be an absolute path
                if key == 'pseudo_dir' and value[0] != '/':
                    value = os.path.abspath(value) + '/'

                try:
                    value = json_ext.loads(value)
                except (TypeError, json_ext.decoder.JSONDecodeError) as e:
                    pass
                calc.parameters[key] = value
        else:
            raise ValueError(f'Unrecognised block "setup:{block}" in the input file')

    # Calculating the simulation cell
    cell = None
    if 'cell_parameters' in dct:
        subdct = dct['cell_parameters']
        cell = utils.read_cell_parameters(calc, subdct)

    # Generating cell if it is missing
    if cell is None:
        _, cell = ibrav_to_cell(calc.parameters)

    # Attaching the cell to the calculator
    calc.atoms = Atoms(cell=cell)

    # Handling atomic species
    if 'atomic_species' in dct:
        utils.read_atomic_species(calc, dct['atomic_species'])

    # Calculating kpoints
    psps_and_kpts: Dict[str, Any] = {}
    if 'k_points' in dct:
        psps_and_kpts.update(**dct['k_points'])

    if task != 'ui':
        def read_compulsory_block(block_name, extract_function):
            if block_name in dct:
                subdct = dct[block_name]
                extract_function(calc, subdct)
                del dct[block_name]
            else:
                raise ValueError(f'{block_name} not found in "setup" block')

        for block_name, extract_function in compulsory_block_readers.items():
            read_compulsory_block(block_name, extract_function)

    # Separamting the output into atoms, parameters, and psp+kpoint information
    atoms = calc.atoms
    atoms.calc = None
    parameters = calc.parameters
    if 'pseudopotentials' in parameters:
        psps_and_kpts['pseudopotentials'] = parameters.pop('pseudopotentials')

    return atoms, parameters, psps_and_kpts


def generate_default_master_calc_params():
    # Dictionary to be used as the default value for 'master_calc_params' when initialising a workflow
    # We create this dynamically in order for the .directory attributes to make sense
    return {'kcp': settings.KoopmansCPSettingsDict(),
            'kc_ham': settings.KoopmansHamSettingsDict(),
            'kc_screen': settings.KoopmansScreenSettingsDict(),
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


def sanitise_master_calc_params(dct_in: Union[Dict[str, Dict], Dict[str, settings.SettingsDict]]) \
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
                f'Unrecognised master_calc_params entry "{k}": valid options are ' + '/'.join(settings_classes.keys()))
    return dct_out
