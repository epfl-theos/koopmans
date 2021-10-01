"""

Generic workflow object for koopmans

Written by Edward Linscott Oct 2020
Converted workflows from functions to objects Nov 2020

"""

import os
import copy
from pathlib import Path
import numpy as np
from typing import Optional, Dict, List, Type, Union, Any
from ase import Atoms
from ase.build import make_supercell
from ase.dft.kpoints import BandPath
from ase.calculators.calculator import CalculationFailed
from ase.calculators.espresso import EspressoWithBandstructure
from koopmans import pseudopotentials, utils, settings
import koopmans.calculators as calculators
from koopmans.commands import ParallelCommandWithPostfix
from koopmans.bands import Bands

valid_settings = [
    settings.Setting('task',
                     'Task to perform',
                     str, 'singlepoint', ('singlepoint', 'convergence', 'wannierise', 'environ_dscf', 'ui')),
    settings.Setting('functional',
                     'orbital-density-dependent-functional/density-functional to use',
                     str, 'ki', ('ki', 'kipz', 'pkipz', 'dft', 'all')),
    settings.Setting('calculate_alpha',
                     'whether or not to calculate the screening parameters ab-initio',
                     bool, True, (True, False)),
    settings.Setting('method',
                     'the method to calculate the screening parameters: either with ΔSCF or DFPT',
                     str, 'dscf', ('dscf', 'dfpt')),
    settings.Setting('init_orbitals',
                     'which orbitals to use as an initial guess for the variational orbitals',
                     str, 'pz', ('pz', 'kohn-sham', 'mlwfs', 'projwfs', 'from old ki')),
    settings.Setting('init_empty_orbitals',
                     'which orbitals to use as an initial guess for the empty variational orbitals '
                     '(defaults to the same value as "init_orbitals")',
                     str, 'same', ('same', 'pz', 'kohn-sham', 'mlwfs', 'projwfs', 'from old ki')),
    settings.Setting('frozen_orbitals',
                     "if True, freeze the variational orbitals for the duration of the calculation once they've been "
                     "initialised",
                     bool, None, (True, False)),
    settings.Setting('periodic',
                     'whether or not the system is periodic',
                     bool, False, (True, False)),
    settings.Setting('npool',
                     'Number of pools for parallelising over kpoints (should be commensurate with the k-point grid)',
                     int, None, None),
    settings.Setting('gb_correction',
                     'if True, apply the Gygi-Baldereschi scheme to deal with the q->0 divergence of the Coulomb '
                     'interation for periodic systems',
                     bool, None, (True, False)),
    settings.Setting('mp_correction',
                     'if True, apply the Makov-Payne correction for charged periodic systems',
                     bool, False, (True, False)),
    settings.Setting('mt_correction',
                     'if True, apply the Martyna-Tuckerman correction for charged aperiodic systems',
                     bool, None, (True, False)),
    settings.Setting('eps_inf',
                     'dielectric constant of the system used by the Gygi-Baldereschi and Makov-Payne corrections',
                     float, None, None),
    settings.Setting('n_max_sc_steps',
                     'maximum number of self-consistency steps for calculating alpha',
                     int, 1, None),
    settings.Setting('alpha_conv_thr',
                     'convergence threshold for |Delta E_i - epsilon_i|; if below this '
                     'threshold, the corresponding alpha value is not updated',
                     (float, str), 1e-3, None),
    settings.Setting('alpha_guess',
                     'starting guess for alpha (overridden if alpha_from_file is true)',
                     (float, list), 0.6, None),
    settings.Setting('alpha_from_file',
                     'if True, uses the file_alpharef.txt from the base directory as a '
                     'starting guess',
                     bool, False, (True, False)),
    settings.Setting('print_qc',
                     'if True, prints out strings for the purposes of quality control',
                     bool, False, (True, False)),
    settings.Setting('from_scratch',
                     'if True, will delete any preexisting workflow and start again; '
                     'if False, will resume a workflow from where it was last up to',
                     bool, False, (True, False)),
    settings.Setting('orbital_groups',
                     'a list of integers the same length as the total number of bands, '
                     'denoting which bands to assign the same screening parameter to',
                     list, None, None),
    settings.Setting('orbital_groups_self_hartree_tol',
                     'when calculating alpha parameters, the code will group orbitals '
                     'together only if their self-Hartree energy is within this '
                     'threshold',
                     float, None, None),
    settings.Setting('enforce_spin_symmetry',
                     'if True, the spin-up and spin-down wavefunctions will be forced '
                     'to be the same',
                     bool, True, (True, False)),
    settings.Setting('check_wannierisation',
                     'if True, checks the Im/Re ratio and generates a plot of the interpolated band structure',
                     bool, False, (True, False)),
    settings.Setting('convergence_observable',
                     'System observable of interest which we converge',
                     str, 'total energy', ('homo energy', 'lumo energy', 'total energy')),
    settings.Setting('convergence_threshold',
                     'Convergence threshold for the system observable of interest',
                     (str, float), None, None),
    settings.Setting('convergence_parameters',
                     'The observable of interest will be converged with respect to this/these '
                     'simulation parameter(s)',
                     (list, str), ['ecutwfc'], None),
    settings.Setting('eps_cavity',
                     'a list of epsilon_infinity values for the cavity in dscf calculations',
                     list, [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20], None)]


class Workflow(object):

    def __init__(self, atoms: Atoms, workflow_settings: settings.SettingsDict,
                 master_calc_params: Dict[str, settings.SettingsDict], name: str,
                 pseudopotentials: Optional[Dict[str, str]] = None,
                 kpts: Optional[List[int]] = [1, 1, 1],
                 koffset: Optional[List[int]] = [0, 0, 0],
                 kpath: Optional[BandPath] = None):

        self.master_calc_params = master_calc_params
        self.atoms = atoms
        self.name = name
        self.calculations: List[calculators.ExtendedCalculator] = []
        self.silent = False
        self.print_indent = 1
        if pseudopotentials is not None:
            self.pseudopotentials = pseudopotentials
        self.kpts = kpts
        self.koffset = koffset
        self.kpath = BandPath(path='G', cell=self.atoms.cell) if kpath is None else kpath

        # If atoms has a calculator, overwrite the kpoints and pseudopotentials variables and then detach the calculator
        if atoms.calc is not None:
            utils.warn(f'You have initialised a {self.__class__.__name__} object with an atoms object that possesses '
                       'a calculator. This calculator will be ignored.')
            self.atoms.calc = None

        # Parsing workflow_settings
        self.parameters = settings.SettingsDictWithChecks(
            settings=valid_settings, physicals=['alpha_conv_thr', 'convergence_threshold'], **workflow_settings)

        # Check internal consistency of workflow settings
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
    def kpts(self):
        return self._kpts

    @kpts.setter
    def kpts(self, value: List[int]):
        self._kpts = value

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
                'workflow_settings': copy.deepcopy(self.parameters),
                'master_calc_params': copy.deepcopy(self.master_calc_params),
                'name': copy.deepcopy(self.name),
                'pseudopotentials': copy.deepcopy(self.pseudopotentials),
                'kpts': copy.deepcopy(self.kpts),
                'koffset': copy.deepcopy(self.koffset),
                'kpath': copy.deepcopy(self.kpath)}

    def new_calculator(self, calc_type: str, directory: Optional[Path] = None, **kwargs) -> calculators.ExtendedCalculator:
        calc_class: Type[calculators.ExtendedCalculator]

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

        # Add pseudopotential and kpt information to the calculator as required
        for kw in ['pseudopotentials', 'kpts', 'koffset', 'kpath']:
            if kw not in all_kwargs and kw in master_calc_params.valid:
                all_kwargs[kw] = getattr(self, kw)

        # Create the calculator
        calc = calc_class(atoms=copy.deepcopy(self.atoms), **all_kwargs)

        # Add the directory if provided
        if directory is not None:
            calc.directory = directory

        if calc_type == 'kcp' and 'pseudopotentials' not in calc.parameters:
            raise ValueError()

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

            for i in range(1, 3):
                to_delete = nspin2_tmpdir / f'{prefix}{i}.{suffix}'
                if to_delete != file_out and to_delete.is_file():
                    to_delete.unlink()

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

    def primitive_to_supercell(self, matrix, **kwargs):
        # Converts to a supercell as given by a 3x3 transformation matrix
        assert np.shape(matrix) == (3, 3)
        self.atoms = make_supercell(self.atoms, matrix, **kwargs)

    def supercell_to_primitive(self, matrix, **kwargs):
        # Converts from a supercell to a primitive cell, as given by a 3x3 transformation matrix
        # The inverse of self.primitive_to_supercell()
        assert np.shape(matrix) == (3, 3)
        raise NotImplementedError('Yet to write this function')

    def run_calculator(self, master_qe_calc, enforce_ss=False):
        '''
        Wrapper for run_calculator_single that manages the optional enforcing of spin symmetry
        '''

        if enforce_ss:
            # Create a copy of the calculator object (to avoid modifying the input)
            qe_calc = copy.deepcopy(master_qe_calc)
            nspin2_tmpdir = master_qe_calc.parameters.outdir / \
                f'{master_qe_calc.parameters.prefix}_{master_qe_calc.parameters.ndw}.save/K00001'

            if master_qe_calc.parameters.restart_mode == 'restart':
                # PBE with nspin=1 dummy
                qe_calc.prefix += '_nspin1_dummy'
                qe_calc.parameters.do_outerloop = False
                qe_calc.parameters.do_outerloop_empty = False
                qe_calc.parameters.nspin, qe_calc.parameters.nelup, qe_calc.parameters.neldw, qe_calc.parameters.tot_magnetization = 1, None, None, None
                qe_calc.parameters.ndw, qe_calc.parameters.ndr = 98, 98
                qe_calc.parameters.restart_mode = 'from_scratch'
                qe_calc.parameters.skip_qc = True
                self.run_calculator_single(qe_calc)
                # Copy over nspin=2 wavefunction to nspin=1 tmp directory (if it has not been done already)
                nspin1_tmpdir = qe_calc.parameters.outdir / \
                    f'{qe_calc.parameters.prefix}_{qe_calc.parameters.ndw}.save/K00001'
                self.convert_wavefunction_2to1(nspin2_tmpdir, nspin1_tmpdir)

            # PBE with nspin=1
            qe_calc = copy.deepcopy(master_qe_calc)
            qe_calc.prefix += '_nspin1'
            qe_calc.parameters.nspin, qe_calc.parameters.nelup, qe_calc.parameters.neldw, qe_calc.parameters.tot_magnetization = 1, None, None, None
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
            qe_calc.parameters.skip_qc = True
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

            calc_file = f'{qe_calc.directory}/{qe_calc.parameters.name}'

            if os.path.isfile(calc_file + qe_calc.ext_out):
                verb = 'Rerunning'

                is_complete = self.load_old_calculator(qe_calc)

                if is_complete:
                    if not self.silent:
                        self.print(f'Not running {os.path.relpath(calc_file)} as it is already complete')
                    return

        # Write out screening parameters to file
        if qe_calc.parameters.get('do_orbdep', False) or isinstance(qe_calc, calculators.KoopmansHamCalculator):
            qe_calc.write_alphas()

        if not self.silent:
            dir_str = os.path.relpath(qe_calc.directory) + '/'
            self.print(f'{verb} {dir_str}{qe_calc.prefix}...', end='', flush=True)

        # Update postfix if relevant
        if self.parameters.npool:
            if isinstance(qe_calc.command, ParallelCommandWithPostfix):
                qe_calc.command.postfix = f'-npool {self.npool}'

        qe_calc.calculate()

        if not qe_calc.is_complete():
            self.print(' failed')
            raise CalculationFailed()

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
        old_calc = qe_calc.__class__(qe_files=f'{qe_calc.directory}/{qe_calc.name}')

        if old_calc.is_complete():
            # If it is complete, load the results
            qe_calc.results = old_calc.results

            # Load bandstructure if present, too
            if isinstance(qe_calc, calculators.UnfoldAndInterpolateCalculator):
                qe_calc.read_bands()
                # If the band structure file does not exist, we must re-run
                if 'band structure' not in qe_calc.results:
                    return False
            raise ValueError('Need to work out how to restructure the code below')
            # elif isinstance(qe_calc.calc, EspressoWithBandstructure):
            #     if not isinstance(qe_calc, calculators.EspressoCalculator) or qe_calc.calculation == 'bands':
            #         qe_calc.band_structure()

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

        # Make sure the subworkflow has inherited the pseudos and kpoints data
        workflow.pseudopotentials = self.pseudopotentials
        workflow.kpts = self.kpts
        workflow.koffset = self.koffset
        workflow.kpath = self.kpath

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
        try:
            workflow.bands = self.bands
            # Only include the most recent values
            workflow.bands.alpha_history = workflow.bands.alphas
            workflow.bands.error_history = []
        except AttributeError:
            pass

        if subdirectory is not None:
            # Update directories
            for key in workflow.master_calc_params.keys():
                params = workflow.master_calc_params[key]
                for setting in params.are_paths:
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
        try:
            self.bands = workflow.bands
        except AttributeError:
            pass

    def todict(self):
        # Shallow copy
        dct = dict(self.__dict__)

        # Adding information required by the json decoder
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @classmethod
    def fromdict(cls, dct):
        raise NotImplementedError('Need to rewrite fromdict using classmethod syntax')

    @property
    def bands(self):
        if not hasattr(self, '_bands'):
            raise AttributeError('Bands have not been initialised')
        return self._bands

    @bands.setter
    def bands(self, value):
        assert isinstance(value, Bands)
        self._bands = value
