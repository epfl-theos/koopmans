"""

Generic workflow object for koopmans

Written by Edward Linscott Oct 2020
Converted workflows from functions to objects Nov 2020

"""

import os
import copy
from pathlib import Path
from ase.io.espresso.utils import cell_to_ibrav
import numpy as np
import numpy.typing as npt
from typing import Optional, Dict, List, Type, Union, Any, TypeVar
from ase import Atoms
from ase.build.supercells import make_supercell
from ase.dft.kpoints import BandPath
from ase.calculators.calculator import CalculationFailed
from koopmans.pseudopotentials import nelec_from_pseudos
from koopmans import utils, settings
import koopmans.calculators as calculators
from koopmans.commands import ParallelCommandWithPostfix
from koopmans.bands import Bands
from koopmans.projections import ProjectionBlocks


T = TypeVar('T', bound='calculators.CalculatorExt')


class Workflow(object):

    def __init__(self, atoms: Atoms,
                 parameters: settings.SettingsDict = settings.WorkflowSettingsDict(),
                 master_calc_params: Dict[str, settings.SettingsDict] = settings.default_master_calc_params,
                 name: str = 'koopmans_workflow',
                 pseudopotentials: Dict[str, str] = {},
                 gamma_only: Optional[bool] = False,
                 kgrid: Optional[List[int]] = [1, 1, 1],
                 koffset: Optional[List[int]] = [0, 0, 0],
                 kpath: Optional[BandPath] = None,
                 projections: Optional[ProjectionBlocks] = None):

        # Parsing parameters
        self.parameters = settings.WorkflowSettingsDict(**parameters)

        self.master_calc_params = master_calc_params
        self.atoms = atoms
        self.name = name
        self.calculations: List[calculators.CalculatorExt] = []
        self.silent = False
        self.print_indent = 1
        self.pseudopotentials = pseudopotentials
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
        self.projections = ProjectionBlocks([], self.atoms) if projections is None else projections

        if 'periodic' in parameters:
            # If "periodic" was explicitly provided, override self.atoms.pbc
            self.atoms.pbc = self.parameters.periodic
        else:
            # If "periodic" was not explicitly provided, use the value from self.atoms.pbc
            self.parameters.periodic = all(self.atoms.pbc)

        if all(self.atoms.pbc):
            self.atoms.wrap(pbc=True)

        # We rely on kcp_params.nelec so make sure this has been initialised
        if 'nelec' not in self.master_calc_params['kcp']:
            pseudo_dir = self.master_calc_params['kcp'].get('pseudo_dir', None)
            self.master_calc_params['kcp'].nelec = nelec_from_pseudos(self.atoms, self.pseudopotentials, pseudo_dir)

        # Generate the kpath
        if kpath is None:
            if self.parameters.periodic:
                # By default, use ASE's default bandpath for this cell (see
                # https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#brillouin-zone-data)
                bl = self.atoms.cell.get_bravais_lattice()
                path = bl.bandpath().path
                npoints = 10 * len(path) - 9 - 29 * path.count(',')
                self.kpath = bl.bandpath(npoints=npoints)
            else:
                self.kpath = BandPath(path='G', cell=self.atoms.cell)
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
            if self.gamma_only:
                all_kwargs['kpts'] = None
            else:
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

    def run(self):
        raise NotImplementedError('This workflow class has not implemented the run() function')

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
                if qe_calc.parameters.do_orbdep:
                    qe_calc.alphas = [qe_calc.alphas[0]]
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
            if qe_calc.parameters.do_orbdep:
                qe_calc.alphas = [qe_calc.alphas[0]]
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
                 kpath=dct.pop('_kpath'))

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
