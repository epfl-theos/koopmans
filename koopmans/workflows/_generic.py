"""

Generic workflow object for koopmans

Written by Edward Linscott Oct 2020
Converted workflows from functions to objects Nov 2020

"""

import os
import copy
from pathlib import Path
import numpy as np
from ase.calculators.calculator import CalculationFailed
from ase.calculators.espresso import EspressoWithBandstructure
from koopmans import utils, settings
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

    def __init__(self, workflow_settings=None, calcs_dct=None, name=None, dct={}):
        if dct:
            assert workflow_settings is None, f'If using the "dct" argument to initialise {self.__class__.__name__}, '
            'do not use any other arguments'
            assert calcs_dct is None, f'If using the "dct" argument to initialise {self.__class__.__name__}, do not '
            'use any other arguments'
            self.fromdict(dct)
        else:
            assert not dct, f'If using the "dct" argument to initialise {self.__class__.__name__}, do not use any '
            'other arguments'
            self.master_calcs = calcs_dct
            self.name = name
            self.all_calcs = []
            self.silent = False
            self.print_indent = 1

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

        # Update postfix for all relevant calculators
        if self.parameters.npool:
            for calc in self.master_calcs.values():
                if isinstance(calc.command, ParallelCommandWithPostfix):
                    calc.command.postfix = f'-npool {self.npool}'

    def new_calculator(self, calc_type, **kwargs):
        if calc_type in self.master_calcs:
            calc = copy.deepcopy(self.master_calcs[calc_type])
            calc.parameters.update(**kwargs)
            return calc
        else:
            raise ValueError(f'Could not find a calculator of type {calc_type}')

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

    def run_calculator(self, master_qe_calc, enforce_ss=False):
        '''
        Wrapper for run_calculator_single that manages the optional enforcing of spin symmetry
        '''

        if enforce_ss:
            # Create a copy of the calculator object (to avoid modifying the input)
            qe_calc = copy.deepcopy(master_qe_calc)
            nspin2_tmpdir = f'{master_qe_calc.parameters.outdir}/{master_qe_calc.parameters.prefix}_{master_qe_calc.parameters.ndw}.save/K00001'

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
                nspin1_tmpdir = f'{qe_calc.parameters.outdir}/{qe_calc.parameters.prefix}_{qe_calc.parameters.ndw}.save/K00001'
                self.convert_wavefunction_2to1(nspin2_tmpdir, nspin1_tmpdir)

            # PBE with nspin=1
            qe_calc = copy.deepcopy(master_qe_calc)
            qe_calc.prefix += '_nspin1'
            qe_calc.parameters.nspin, qe_calc.parameters.nelup, qe_calc.parameters.neldw, qe_calc.parameters.tot_magnetization = 1, None, None, None
            qe_calc.parameters.ndw, qe_calc.parameters.ndr = 98, 98
            nspin1_tmpdir = f'{qe_calc.parameters.outdir}/{qe_calc.parameters.prefix}_{qe_calc.parameters.ndw}.save/K00001'
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
            nspin2_tmpdir = f'{qe_calc.parameters.outdir}/{qe_calc.parameters.prefix}_{qe_calc.parameters.ndw}.save/K00001'
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
        self.all_calcs.append(qe_calc)

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

            self.all_calcs.append(qe_calc)

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
            calc_str = f'{calc.name}_'

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

        # Ensure altering workflow.master_calcs won't affect self.master_calcs
        if workflow.master_calcs is self.master_calcs:
            workflow.master_calcs = copy.deepcopy(self.master_calcs)

        # Setting from_scratch to a non-None value will override the value of subworkflow.from_scratch...
        if from_scratch is None:
            workflow.parameters.from_scratch = self.parameters.from_scratch
        else:
            workflow.parameters.from_scratch = from_scratch

        # Link the list of calculations
        workflow.all_calcs = self.all_calcs

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
            for key in workflow.master_calcs.keys():
                calc = workflow.master_calcs[key]
                for setting in calc.parameters.are_paths:
                    path = getattr(calc.parameters, setting, None)
                    if path is not None and path.startswith(os.getcwd()):
                        new_path = os.path.abspath('./' + subdirectory + '/' + os.path.relpath(path))
                        setattr(calc.parameters, setting, new_path)

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

        # Removing keys we won't need to reconstruct the workflow
        del dct['valid_settings']

        # Adding information required by the json decoder
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    def fromdict(self, dct):
        for k, v in dct.items():
            setattr(self, k, v)

    @property
    def bands(self):
        if not hasattr(self, '_bands'):
            raise AttributeError('Bands have not been initialised')
        return self._bands

    @bands.setter
    def bands(self, value):
        assert isinstance(value, Bands)
        self._bands = value
