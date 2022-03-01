"""

Generic workflow object for koopmans

Written by Edward Linscott Oct 2020
Converted workflows from functions to objects Nov 2020

"""

import os
import copy
import operator
from functools import reduce
import subprocess
from pathlib import Path
import json as json_ext
import numpy as np
from numpy import typing as npt
from types import ModuleType
from typing import Optional, Dict, List, Type, Union, Any, TypeVar
from pybtex.database import BibliographyData
import ase
from ase import Atoms
from ase.build.supercells import make_supercell
from ase.dft.kpoints import BandPath
from ase.calculators.espresso import Espresso_kcp
from ase.calculators.calculator import CalculationFailed
from ase.io.espresso.utils import cell_to_ibrav, ibrav_to_cell
from ase.io.espresso.koopmans_cp import KEYS as kcp_keys, construct_namelist
from koopmans.pseudopotentials import nelec_from_pseudos, pseudos_library_directory, pseudo_database, fetch_pseudo
from koopmans import utils, settings
import koopmans.calculators as calculators
from koopmans.commands import ParallelCommandWithPostfix
from koopmans.bands import Bands
from koopmans.projections import ProjectionBlocks
from koopmans.references import bib_data
from abc import ABC, abstractmethod


T = TypeVar('T', bound='calculators.CalculatorExt')


class Workflow(ABC):

    def __init__(self, atoms: Atoms,
                 parameters: settings.SettingsDict = settings.WorkflowSettingsDict(),
                 master_calc_params: Union[Dict[str, Dict], Dict[str, settings.SettingsDict]
                                           ] = settings.default_master_calc_params,
                 name: str = 'koopmans_workflow',
                 pseudopotentials: Dict[str, str] = {},
                 gamma_only: Optional[bool] = False,
                 kgrid: Optional[List[int]] = [1, 1, 1],
                 koffset: Optional[List[int]] = [0, 0, 0],
                 kpath: Optional[Union[BandPath, str]] = None,
                 kpath_density: int = 10,
                 projections: Optional[ProjectionBlocks] = None):

        # Parsing parameters
        self.parameters = settings.WorkflowSettingsDict(**parameters)
        self.atoms = atoms
        self.name = name
        self.calculations: List[calculators.CalculatorExt] = []
        self.silent = False
        self.print_indent = 1
        self.gamma_only = gamma_only
        if self.gamma_only:
            if kgrid != [1, 1, 1]:
                utils.warn(f'You have initialised kgrid to {kgrid}, not compatible with gamma_only=True; '
                           'kgrid is set equal to [1, 1, 1]')
            if koffset != [0, 0, 0]:
                utils.warn(f'You have initialised koffset to {koffset}, not compatible with gamma_only=True; '
                           'koffset is set equal to [0, 0, 0]')
            self.kgrid = [1, 1, 1]
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
        elif self.parameters.pseudo_library:
            self.pseudopotentials = {}
            for symbol in set(self.atoms.symbols):
                pseudo = fetch_pseudo(element=symbol, functional=self.parameters.base_functional,
                                      library=self.parameters.pseudo_library)
                if pseudo.kind == 'unknown':
                    utils.warn(f'You are using an unrecognised pseudopotential {pseudo.name}. Please note that '
                               'the current implementation of Koopmans functionals only supports norm-conserving '
                               'pseudopotentials.')
                elif pseudo.kind != 'norm-conserving':
                    raise ValueError('Koopmans functionals only currently support norm-conserving pseudopotentials; '
                                     f'{pseudo.name} is {pseudo.kind}')
                self.pseudopotentials[symbol] = pseudo.name
        else:
            self.pseudopotentials = pseudopotentials

        # Make sure master_calc_params isn't missing any entries, and every entry corresponds to settings.SettingsDict
        # objects
        master_calc_params = sanitise_master_calc_params(master_calc_params)

        # Work out the pseudopotential directory. If using a pseudo_library this is straightforward, if not...
        #  1. try to locating the directory as currently specified by the calculator
        #  2. if that fails, check if $ESPRESSO_PSEUDO is set
        #  3. if that fails, raise an error
        if self.parameters.pseudo_library:
            pseudo_dir = pseudos_library_directory(self.parameters.pseudo_library, self.parameters.base_functional)
            for params in master_calc_params.values():
                if params.get('pseudo_dir', pseudo_dir).resolve() != pseudo_dir:
                    raise ValueError(
                        '"pseudo_dir" and "pseudo_library" are conflicting; please do not provide "pseudo_dir"')
        elif 'pseudo_dir' in master_calc_params['kcp'] or 'pseudo_dir' in master_calc_params['pw']:
            pseudo_dir = master_calc_params['kcp'].get('pseudo_dir', master_calc_params['pw'].get('pseudo_dir', None))
            if not os.path.isdir(pseudo_dir):
                raise NotADirectoryError(f'The pseudo_dir you provided ({pseudo_dir}) does not exist')
        elif 'ESPRESSO_PSEUDO' in os.environ:
            pseudo_dir = Path(os.environ['ESPRESSO_PSEUDO'])
        else:
            pseudo_dir = Path.cwd()

        if self.parameters.task != 'ui':
            for pseudo in self.pseudopotentials.values():
                if not (pseudo_dir / pseudo).exists():
                    raise FileNotFoundError(
                        f'{pseudo_dir / pseudo} does not exist. Please double-check your pseudopotential settings')

        # Before saving the master_calc_params, automatically generate some keywords and perform some sanity checks
        if self.parameters.task != 'ui':
            # Automatically calculate nelec/nelup/neldw/etc using information contained in the pseudopotential files
            # and the kcp settings
            nelec = nelec_from_pseudos(self.atoms, self.pseudopotentials, pseudo_dir)

            tot_charge = master_calc_params['kcp'].get('tot_charge', 0)
            nelec -= tot_charge
            tot_mag = master_calc_params['kcp'].get('tot_magnetization', nelec % 2)
            nelup = int(nelec / 2 + tot_mag / 2)
            neldw = int(nelec / 2 - tot_mag / 2)
            if tot_mag != 0:
                atoms.set_initial_magnetic_moments([tot_mag / len(atoms) for _ in atoms])

            # Work out the number of filled and empty bands
            n_filled = nelec // 2 + nelec % 2
            n_empty = master_calc_params['kcp'].get('empty_states_nbnd', 0)
            nbnd = n_filled + n_empty
            generated_keywords = {'nelec': nelec, 'tot_charge': tot_charge, 'tot_magnetization': tot_mag,
                                  'nelup': nelup, 'neldw': neldw, 'nbnd': nbnd, 'pseudo_dir': pseudo_dir}
        else:
            generated_keywords = {}
            nelec = 0

        self.master_calc_params = settings.default_master_calc_params.copy()
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

        # Records whether or not this workflow is a subworkflow of another
        self._is_a_subworkflow = False

    def __eq__(self, other):
        if isinstance(other, Workflow):
            return self.__dict__ == other.__dict__
        return False

    def run(self) -> None:
        self.print_preamble()
        self._run()
        self.print_conclusion()

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
    def gamma_only(self):
        return self._gamma_only

    @gamma_only.setter
    def gamma_only(self, value: bool):
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
    def wf_kwargs(self):
        # Returns a kwargs designed to be used to initialise another workflow with the same configuration as this one
        # i.e.
        # > sub_wf = Workflow(**self.wf_kwargs)
        return {'atoms': copy.deepcopy(self.atoms),
                'parameters': copy.deepcopy(self.parameters),
                'master_calc_params': copy.deepcopy(self.master_calc_params),
                'name': copy.deepcopy(self.name),
                'pseudopotentials': copy.deepcopy(self.pseudopotentials),
                'gamma_only': copy.deepcopy(self.gamma_only),
                'kgrid': copy.deepcopy(self.kgrid),
                'kpath': copy.deepcopy(self.kpath),
                'projections': copy.deepcopy(self.projections)}

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
        elif calc_type.startswith('ui'):
            calc_class = calculators.UnfoldAndInterpolateCalculator
        elif calc_type == 'wann2kc':
            calc_class = calculators.Wann2KCCalculator
        elif calc_type == 'kc_screen':
            calc_class = calculators.KoopmansScreenCalculator
        elif calc_type == 'kc_ham':
            calc_class = calculators.KoopmansHamCalculator
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
        for kw in ['pseudopotentials', 'gamma_only', 'kgrid', 'kpath', 'koffset']:
            if kw not in all_kwargs and kw in master_calc_params.valid:
                all_kwargs[kw] = getattr(self, kw)

        # Create the calculator
        calc = calc_class(atoms=copy.deepcopy(self.atoms), **all_kwargs)

        # Add the directory if provided
        if directory is not None:
            calc.directory = directory

        return calc

    def convert_wavefunction_2to1(self, nspin2_tmpdir: Path, nspin1_tmpdir: Path):

        if not self.parameters.from_scratch:
            return

        for directory in [nspin2_tmpdir, nspin1_tmpdir]:
            if not directory.is_dir():
                raise OSError(f'{directory} not found')

        for wfile in ['evc0.dat', 'evc0_empty1.dat', 'evcm.dat', 'evc.dat', 'evcm.dat', 'hamiltonian.xml',
                      'eigenval.xml', 'evc_empty1.dat', 'lambda01.dat', 'lambdam1.dat']:
            if '1.' in wfile:
                prefix, suffix = wfile.split('1.')
            else:
                prefix, suffix = wfile.split('.')

            file_out = nspin1_tmpdir / wfile
            file_in = nspin2_tmpdir / f'{prefix}1.{suffix}'

            if file_in.is_file():

                with open(file_in, 'rb') as fd:
                    contents = fd.read()

                contents = contents.replace(b'nk="2"', b'nk="1"')
                contents = contents.replace(b'nspin="2"', b'nspin="1"')

                with open(file_out, 'wb') as fd:
                    fd.write(contents)

    def convert_wavefunction_1to2(self, nspin1_tmpdir: Path, nspin2_tmpdir: Path):

        if not self.parameters.from_scratch:
            return

        for directory in [nspin2_tmpdir, nspin1_tmpdir]:
            if not directory.is_dir():
                raise OSError(f'{directory} not found')

        for wfile in ['evc0.dat', 'evc0_empty1.dat', 'evcm.dat', 'evc.dat', 'evcm.dat', 'hamiltonian.xml',
                      'eigenval.xml', 'evc_empty1.dat', 'lambda01.dat']:
            if '1.' in wfile:
                prefix, suffix = wfile.split('1.')
            else:
                prefix, suffix = wfile.split('.')

            file_in = nspin1_tmpdir / wfile

            if file_in.is_file():
                with open(file_in, 'rb') as fd:
                    contents = fd.read()

                contents = contents.replace(b'nk="1"', b'nk="2"')
                contents = contents.replace(b'nspin="1"', b'nspin="2"')

                file_out = nspin2_tmpdir / f'{prefix}1.{suffix}'
                with open(file_out, 'wb') as fd:
                    fd.write(contents)

                contents = contents.replace(b'ik="1"', b'ik="2"')
                contents = contents.replace(b'ispin="1"', b'ispin="2"')

                file_out = nspin2_tmpdir / f'{prefix}2.{suffix}'
                with open(file_out, 'wb') as fd:
                    fd.write(contents)

    def update_celldms(self):
        # Update celldm(*) to match the current self.atoms.cell
        for k, params in self.master_calc_params.items():
            if 'ibrav' in params:
                celldms = cell_to_ibrav(self.atoms.cell, params.ibrav)
                self.master_calc_params[k].update(**celldms)

    def primitive_to_supercell(self, matrix: Optional[npt.NDArray[np.int_]] = None, **kwargs):
        # Converts to a supercell as given by a 3x3 transformation matrix
        if matrix is None:
            matrix = np.diag(self.kgrid)
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
            if not isinstance(master_qe_calc, calculators.KoopmansCPCalculator):
                raise NotImplementedError('Workflow.run_calculator(..., enforce_ss = True) needs to be generalised to '
                                          'other calculator types')
            # Create a copy of the calculator object (to avoid modifying the input)
            qe_calc = copy.deepcopy(master_qe_calc)
            nspin2_tmpdir = master_qe_calc.parameters.outdir / \
                f'{master_qe_calc.parameters.prefix}_{master_qe_calc.parameters.ndw}.save/K00001'

            if master_qe_calc.parameters.restart_mode == 'restart':
                # PBE with nspin=1 dummy
                qe_calc.prefix += '_nspin1_dummy'
                qe_calc.parameters.do_outerloop = False
                qe_calc.parameters.do_outerloop_empty = False
                qe_calc.parameters.nspin = 1
                if hasattr(qe_calc, 'alphas'):
                    qe_calc.alphas = [qe_calc.alphas[0]]
                if hasattr(qe_calc, 'filling'):
                    qe_calc.filling = [qe_calc.filling[0]]
                qe_calc.parameters.nelup = None
                qe_calc.parameters.neldw = None
                qe_calc.parameters.tot_magnetization = None
                qe_calc.parameters.ndw, qe_calc.parameters.ndr = 98, 98
                qe_calc.parameters.restart_mode = 'from_scratch'
                qe_calc.skip_qc = True
                self.run_calculator_single(qe_calc)
                # Copy over nspin=2 wavefunction to nspin=1 tmp directory (if it has not been done already)
                nspin1_tmpdir = qe_calc.parameters.outdir / \
                    f'{qe_calc.parameters.prefix}_{qe_calc.parameters.ndw}.save/K00001'
                self.convert_wavefunction_2to1(nspin2_tmpdir, nspin1_tmpdir)

            # PBE with nspin=1
            qe_calc = copy.deepcopy(master_qe_calc)
            qe_calc.prefix += '_nspin1'
            qe_calc.parameters.nspin = 1
            qe_calc.parameters.nelup = None
            qe_calc.parameters.neldw = None
            if hasattr(qe_calc, 'alphas'):
                qe_calc.alphas = [qe_calc.alphas[0]]
            if hasattr(qe_calc, 'filling'):
                qe_calc.filling = [qe_calc.filling[0]]
            qe_calc.parameters.tot_magnetization = None
            qe_calc.parameters.ndw, qe_calc.parameters.ndr = 98, 98
            nspin1_tmpdir = qe_calc.parameters.outdir / \
                f'{qe_calc.parameters.prefix}_{qe_calc.parameters.ndw}.save/K00001'
            self.run_calculator_single(qe_calc)

            # PBE from scratch with nspin=2 (dummy run for creating files of appropriate size)
            qe_calc = copy.deepcopy(master_qe_calc)
            qe_calc.prefix += '_nspin2_dummy'
            qe_calc.parameters.restart_mode = 'from_scratch'
            qe_calc.parameters.do_outerloop = False
            qe_calc.parameters.do_outerloop_empty = False
            qe_calc.parameters.ndw = 99
            qe_calc.skip_qc = True
            self.run_calculator_single(qe_calc)

            # Copy over nspin=1 wavefunction to nspin=2 tmp directory (if it has not been done already)
            nspin2_tmpdir = qe_calc.parameters.outdir / \
                f'{qe_calc.parameters.prefix}_{qe_calc.parameters.ndw}.save/K00001'
            self.convert_wavefunction_1to2(nspin1_tmpdir, nspin2_tmpdir)

            # PBE with nspin=2, reading in the spin-symmetric nspin=1 wavefunction
            master_qe_calc.prefix += '_nspin2'
            master_qe_calc.parameters.restart_mode = 'restart'
            master_qe_calc.parameters.ndr = 99
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
                    return

        if not self.silent:
            dir_str = os.path.relpath(qe_calc.directory) + '/'
            self.print(f'{verb} {dir_str}{qe_calc.prefix}...', end='', flush=True)

        # Update postfix if relevant
        if self.parameters.npool:
            if isinstance(qe_calc.command, ParallelCommandWithPostfix):
                qe_calc.command.postfix = f'-npool {self.parameters.npool}'

        qe_calc.calculate()

        if not qe_calc.is_complete():
            self.print(' failed')
            raise CalculationFailed(
                f'{qe_calc.directory}/{qe_calc.prefix} failed; check the Quantum ESPRESSO output file for more details')

        if not self.silent:
            self.print(' done')

        # Check spin-up and spin-down eigenvalues match
        if 'eigenvalues' in qe_calc.results and isinstance(qe_calc, calculators.KoopmansCPCalculator):
            if qe_calc.is_converged() and qe_calc.parameters.do_outerloop and qe_calc.parameters.nspin == 2 \
                    and qe_calc.parameters.tot_magnetization == 0 and not qe_calc.parameters.fixed_state \
                    and len(qe_calc.results['eigenvalues']) > 0:
                rms_eigenval_difference = np.sqrt(
                    np.mean(np.diff(qe_calc.results['eigenvalues'], axis=0)**2))
                if rms_eigenval_difference > 0.05:
                    utils.warn('Spin-up and spin-down eigenvalues differ substantially')

        # Store the calculator
        self.calculations.append(qe_calc)

        # Print quality control
        if self.parameters.print_qc and not qe_calc.skip_qc:
            for result in qe_calc.results_for_qc:
                val = qe_calc.results.get(result, None)
                if val:
                    self.print_qc_keyval(result, val, qe_calc)

        # If we reached here, all future calculations should be performed from scratch
        self.parameters.from_scratch = True

        return

    def load_old_calculator(self, qe_calc):
        # This is a separate function so that it can be monkeypatched by the test suite
        old_calc = qe_calc.__class__.fromfile(qe_calc.directory / qe_calc.prefix)

        if old_calc.is_complete():
            # If it is complete, load the results
            qe_calc.results = old_calc.results

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

    def print_qc_keyval(self, key, value, calc=None):
        '''
        Prints out a quality control message for testcode to evaluate, and if a calculator is provided, stores the
        result in calc.qc_results
        '''

        if calc is None:
            calc_str = ''
        else:
            calc_str = f'{calc.prefix}_'

        if not isinstance(value, list):
            self.print(f'<QC> {calc_str}{key} {value}')

        if calc is not None:
            calc.qc_results[key] = value

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
                 projections=dct.pop('projections'))

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

        # Check for unexpected blocks
        for block in bigdct:
            if block not in list(settings_classes.keys()) + ['workflow', 'setup']:
                raise ValueError(f'Unrecognised block "{block}" in json input file; '
                                 'valid options are workflow/' + '/'.join(settings_classes.keys()))

        # Loading workflow settings
        parameters = settings.WorkflowSettingsDict(**utils.parse_dict(bigdct.get('workflow', {})))

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
                **{k: v for k, v in setup_parameters.items() if k in master_calc_params[block].valid})

        # Adding the projections to the workflow kwargs (this is unusual in that this is an attribute of the workflow
        # object but it is provided in the w90 subdictionary)
        workflow_kwargs['projections'] = ProjectionBlocks.fromprojections(
            w90_block_projs, w90_block_filling, w90_block_spins, atoms)

        name = fname.replace('.json', '')

        return cls(atoms, parameters, master_calc_params, name, **workflow_kwargs)

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
        if self.atoms.has('labels'):
            labels = self.atoms.get_array('labels')
        else:
            labels = self.atoms.get_chemical_symbols()
        if ibrav == 0:
            bigdct['setup']['atomic_positions'] = {'positions': [
                [label] + [str(x) for x in pos] for label, pos in zip(labels, self.atoms.get_positions())],
                'units': 'angstrom'}
        else:
            bigdct['setup']['atomic_positions'] = {'positions': [
                [label] + [str(x) for x in pos] for label, pos in zip(labels, self.atoms.get_scaled_positions())],
                'units': 'crystal'}

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
            if not params_dict:
                continue

            if code in ['pw', 'kcp']:
                bigdct[code] = {}

                # Populate bigdct with the settings
                input_data = construct_namelist(params_dict)
                for key, block in input_data.items():

                    if len(block) > 0:
                        bigdct[code][key] = {k: v for k, v in dict(
                            block).items() if v is not None}

            elif code in ['pw2wannier', 'wann2kc', 'kc_screen', 'kc_ham']:
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
                else:
                    proj_kwarg = {'projections': projections[0].projections}
                reduce(operator.getitem, nested_keys[:-1], bigdct['w90'])[k].update(**proj_kwarg)

            else:
                raise NotImplementedError(
                    f'Writing of {params.__class__.__name__} with write_json is not yet implemented')

        return bigdct


def get_version(module):
    if isinstance(module, ModuleType):
        module = module.__path__[0]
    with utils.chdir(module):
        version_label = subprocess.check_output(["git", "describe", "--always", "--tags"]).strip()
    return version_label.decode("utf-8")


def header():

    koopmans_version = get_version(os.path.dirname(__file__))
    ase_version = get_version(ase)
    qe_version = get_version(calculators.qe_bin_directory)

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


# Define which function to use to read each block
settings_classes = {'kcp': settings.KoopmansCPSettingsDict,
                    'kc_ham': settings.KoopmansHamSettingsDict,
                    'kc_screen': settings.KoopmansScreenSettingsDict,
                    'wann2kc': settings.Wann2KCSettingsDict,
                    'pw': settings.PWSettingsDict,
                    'pw2wannier': settings.PW2WannierSettingsDict,
                    'ui': settings.UnfoldAndInterpolateSettingsDict,
                    'ui_occ': settings.UnfoldAndInterpolateSettingsDict,
                    'ui_emp': settings.UnfoldAndInterpolateSettingsDict,
                    'w90_occ': settings.Wannier90SettingsDict,
                    'w90_emp': settings.Wannier90SettingsDict,
                    'w90_occ_up': settings.Wannier90SettingsDict,
                    'w90_emp_up': settings.Wannier90SettingsDict,
                    'w90_occ_down': settings.Wannier90SettingsDict,
                    'w90_emp_down': settings.Wannier90SettingsDict}


def sanitise_master_calc_params(dct_in: Union[Dict[str, Dict], Dict[str, settings.SettingsDict]]) \
        -> Dict[str, settings.SettingsDict]:
    dct_out: Dict[str, settings.SettingsDict]
    for k, cls in settings_classes.items():
        dct: Union[Dict, settings.SettingsDict] = dct_in.get(k, {})
        if isinstance(dct, settings.SettingsDict):
            dct_out[k] = dct
            if dct_out[k].directory == '':
                dct_out[k].directory = Path().resolve()
        else:
            dct_out[k] = cls(**dct, directory=Path().resolve())

    for k in dct_in.keys():
        if k not in settings_classes:
            raise ValueError(
                f'Unrecognised master_calc_params entry "{k}": valid options are ' + '/'.join(settings_classes.keys()))
    return dct_out
