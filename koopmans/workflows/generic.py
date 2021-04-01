"""

Generic workflow object for python_KI

Written by Edward Linscott Oct 2020
Converted workflows from functions to objects Nov 2020

"""

import os
import sys
import copy
import numpy as np
from ase.calculators.calculator import CalculationFailed
from koopmans import io, utils
from koopmans.bands import Bands

valid_settings = [
    utils.Setting('task',
                  'Task to perform',
                  str, 'singlepoint', ('singlepoint', 'convergence', 'environ_dscf', 'ui')),
    utils.Setting('functional',
                  'Orbital-density-dependent-functional/density-functional to use',
                  str, 'ki', ('ki', 'kipz', 'pkipz', 'pbe', 'all')),
    utils.Setting('init_orbitals',
                  'which orbitals to use as an initial guess for the variational orbitals',
                  str, 'pz', ('pz', 'kohn-sham', 'mlwfs', 'projwfs', 'from old ki')),
    utils.Setting('init_empty_orbitals',
                  'which orbitals to use as an initial guess for the empty variational orbitals '
                  '(defaults to the same value as "init_orbitals")',
                  str, 'same', ('same', 'pz', 'kohn-sham', 'mlwfs', 'projwfs', 'from old ki')),
    utils.Setting('periodic',
                  'whether or not the system is periodic. If False, interaction between '
                  'periodic images will be corrected for',
                  bool, False, (True, False)),
    utils.Setting('mp_corrections',
                  'if True, the Makov-Payne corrections for charged systems will '
                  'be applied',
                  bool, False, (True, False)),
    utils.Setting('eps_inf',
                  'dielectric constant of the system; needed when mp_corrections is True',
                  float, None, None),
    utils.Setting('n_max_sc_steps',
                  'maximum number of self-consistency steps for calculating alpha',
                  int, 1, None),
    utils.Setting('alpha_conv_thr',
                  'convergence threshold for |delta E - lambda|; if below this '
                  'threshold, the corresponding alpha value is not updated',
                  (float, str), 1e-3, None),
    utils.Setting('calculate_alpha',
                  'if True, the screening parameters will be calculated; if False, '
                  'they will be read directly from file',
                  bool, True, (True, False)),
    utils.Setting('alpha_guess',
                  'starting guess for alpha (overridden if alpha_from_file is true)',
                  float, 0.6, None),
    utils.Setting('alpha_from_file',
                  'if True, uses the file_alpharef.txt from the base directory as a '
                  'starting guess',
                  bool, False, (True, False)),
    utils.Setting('print_qc',
                  'if True, prints out strings for the purposes of quality control',
                  bool, False, (True, False)),
    utils.Setting('from_scratch',
                  'if True, will delete any preexisting workflow and start again; '
                  'if False, will resume a workflow from where it was last up to',
                  bool, False, (True, False)),
    utils.Setting('orbital_groups',
                  'a list of integers the same length as the total number of bands, '
                  'denoting which bands to assign the same screening parameter to',
                  list, None, None),
    utils.Setting('enforce_spin_symmetry',
                  'if True, the spin-up and spin-down wavefunctions will be forced '
                  'to be the same',
                  bool, True, (True, False)),
    utils.Setting('convergence_observable',
                  'System observable of interest which we converge',
                  str, 'total energy', ('homo energy', 'lumo energy', 'total energy')),
    utils.Setting('convergence_threshold',
                  'Convergence threshold for the system observable of interest',
                  (str, float), None, None),
    utils.Setting('convergence_parameters',
                  'The observable of interest will be converged with respect to this/these '
                  'simulation parameter(s)',
                  (list, str), ['ecutwfc'], None),
    utils.Setting('eps_cavity',
                  'a list of epsilon_infinity values for the cavity in dscf calculations',
                  list, [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20], None)]


class Workflow(object):

    def __init__(self, workflow_settings=None, calcs_dct=None, dct={}):
        if dct:
            assert workflow_settings is None, f'If using the "dct" argument to initialise {self.__class__.__name__}, ' \
                'do not use any other arguments'
            assert calcs_dct is None, f'If using the "dct" argument to initialise {self.__class__.__name__}, do not ' \
                'use any other arguments'
            self.fromdict(dct)
        else:
            assert not dct, f'If using the "dct" argument to initialise {self.__class__.__name__}, do not use any ' \
                'other arguments'
            self.master_calcs = calcs_dct
            self.all_calcs = []
            self.silent = False
            self.valid_settings = valid_settings

            # Parsing workflow_settings
            checked_settings = utils.check_settings(workflow_settings, self.valid_settings, physicals=[
                'alpha_conv_thr', 'convergence_threshold'])
            self.list_of_settings = list(checked_settings.keys())
            for key, value in checked_settings.items():
                self.add_setting(key, value)

    @property
    def settings(self):
        return {k: getattr(self, k) for k in self.list_of_settings}

    def add_setting(self, key, val):
        if key not in self.list_of_settings:
            self.list_of_settings.append(key)
        setattr(self, key, val)

    def new_calculator(self, calc_type, **kwargs):
        if calc_type in self.master_calcs:
            calc = copy.deepcopy(self.master_calcs[calc_type])
            for keyword, value in kwargs.items():
                if keyword not in calc._valid_settings and not hasattr(calc, keyword):
                    raise ValueError(f'{keyword} is not a valid setting name')
                setattr(calc, keyword, value)
            return calc
        else:
            raise ValueError(f'Could not find a calculator of type {calc_type}')

    def run(self):
        raise NotImplementedError('This workflow class has not implemented the run() function')

    def convert_wavefunction_2to1(self, nspin2_tmpdir, nspin1_tmpdir):
        if self.from_scratch:
            utils.system_call(
                f'convert_nspin2_wavefunction_to_nspin1.sh {nspin2_tmpdir} {nspin1_tmpdir}')

    def convert_wavefunction_1to2(self, nspin1_tmpdir, nspin2_tmpdir):
        if self.from_scratch:
            utils.system_call(
                f'convert_nspin1_wavefunction_to_nspin2.sh {nspin1_tmpdir} {nspin2_tmpdir}')

    def run_calculator(self, master_qe_calc, enforce_ss=False):
        '''
        Wrapper for run_calculator_single that manages the optional enforcing of spin symmetry
        '''

        if enforce_ss:
            # Create a copy of the calculator object (to avoid modifying the input)
            qe_calc = copy.deepcopy(master_qe_calc)
            nspin2_tmpdir = f'{master_qe_calc.outdir}/{master_qe_calc.prefix}_{master_qe_calc.ndw}.save/K00001'

            if master_qe_calc.restart_mode == 'restart':
                # PBE with nspin=1 dummy
                qe_calc.name += '_nspin1_dummy'
                qe_calc.do_outerloop = False
                qe_calc.do_outerloop_empty = False
                qe_calc.nspin, qe_calc.nelup, qe_calc.neldw, qe_calc.tot_magnetization = 1, None, None, None
                qe_calc.ndw, qe_calc.ndr = 98, 98
                qe_calc.restart_mode = 'from_scratch'
                qe_calc.skip_qc = True
                self.run_calculator_single(qe_calc)
                # Copy over nspin=2 wavefunction to nspin=1 tmp directory (if it has not been done already)
                nspin1_tmpdir = f'{qe_calc.outdir}/{qe_calc.prefix}_{qe_calc.ndw}.save/K00001'
                self.convert_wavefunction_2to1(nspin2_tmpdir, nspin1_tmpdir)

            # PBE with nspin=1
            qe_calc = copy.deepcopy(master_qe_calc)
            qe_calc.name += '_nspin1'
            qe_calc.nspin, qe_calc.nelup, qe_calc.neldw, qe_calc.tot_magnetization = 1, None, None, None
            qe_calc.ndw, qe_calc.ndr = 98, 98
            nspin1_tmpdir = f'{qe_calc.outdir}/{qe_calc.prefix}_{qe_calc.ndw}.save/K00001'
            self.run_calculator_single(qe_calc)

            # PBE from scratch with nspin=2 (dummy run for creating files of appropriate size)
            qe_calc = copy.deepcopy(master_qe_calc)
            qe_calc.name += '_nspin2_dummy'
            qe_calc.restart_mode = 'from_scratch'
            qe_calc.do_outerloop = False
            qe_calc.do_outerloop_empty = False
            qe_calc.ndw = 99
            qe_calc.skip_qc = True
            self.run_calculator_single(qe_calc)

            # Copy over nspin=1 wavefunction to nspin=2 tmp directory (if it has not been done already)
            nspin2_tmpdir = f'{qe_calc.outdir}/{qe_calc.prefix}_{qe_calc.ndw}.save/K00001'
            self.convert_wavefunction_1to2(nspin1_tmpdir, nspin2_tmpdir)

            # PBE with nspin=2, reading in the spin-symmetric nspin=1 wavefunction
            master_qe_calc.name += '_nspin2'
            master_qe_calc.restart_mode = 'restart'
            master_qe_calc.ndr = 99
            self.run_calculator_single(master_qe_calc)

        else:

            self.run_calculator_single(master_qe_calc)

        return

    def run_calculator_single(self, qe_calc):
        # Runs qe_calc.calculate with additional checks

        # If an output file already exists, check if the run completed successfully
        verb = 'Running'
        if not self.from_scratch:

            calc_file = f'{qe_calc.directory}/{qe_calc.name}'

            if os.path.isfile(calc_file + qe_calc.ext_out):
                verb = 'Rerunning'

                self.load_old_calculator(qe_calc)

                old_calc = self.all_calcs[-1]
                if old_calc.is_complete():
                    # If it is complete, load the results, and exit
                    qe_calc.results = old_calc.results
                    if not self.silent:
                        print(f'Not running {calc_file} as it is already complete')
                    self.all_calcs.append(qe_calc)
                    return

        # Write out screening parameters to file
        if getattr(qe_calc, 'do_orbdep', False):
            qe_calc.write_alphas()

        if not self.silent:
            if qe_calc.directory == '.':
                dir_str = ''
            else:
                dir_str = qe_calc.directory + '/'
            print(f'{verb} {dir_str}{qe_calc.name}...', end='', flush=True)

        qe_calc.calculate()

        if not qe_calc.is_complete():
            print(' failed')
            raise CalculationFailed()

        if not self.silent:
            print(' done')

        # Check spin-up and spin-down eigenvalues match
        if 'eigenvalues' in qe_calc.results:

            if qe_calc.is_converged() and qe_calc.do_outerloop and qe_calc.nspin == 2 \
                    and qe_calc.tot_magnetization == 0 and not qe_calc.fixed_state \
                    and len(qe_calc.results['eigenvalues']) > 0:
                rms_eigenval_difference = np.sqrt(
                    np.mean(np.diff(qe_calc.results['eigenvalues'], axis=0)**2))
                if rms_eigenval_difference > 0.05:
                    utils.warn('Spin-up and spin-down eigenvalues differ substantially')

        # Store the calculator
        self.all_calcs.append(qe_calc)

        # Print quality control
        if self.print_qc and not qe_calc.skip_qc:
            for result in qe_calc.results_for_qc:
                val = qe_calc.results.get(result, None)
                if val:
                    self.print_qc_keyval(result, val, qe_calc)

        # If we reached here, all future calculations should be performed from scratch
        self.from_scratch = True

        return

    def load_old_calculator(self, qe_calc):
        # This is a separate function so that it can be monkeypatched by the test suite
        calc_file = f'{qe_calc.directory}/{qe_calc.name}'
        self.all_calcs.append(qe_calc.__class__(qe_files=calc_file))

    def print_qc_keyval(self, key, value, calc=None):
        '''
        Prints out a quality control message for testcode to evaluate, and if a calculator is provided, stores the
        result in calc.qc_results
        '''

        if calc is None:
            calc_str = ''
        else:
            calc_str = f'{calc.name}_'
        print(f'<QC> {calc_str}{key} {value}')

        if calc is not None:
            calc.qc_results[key] = value

    def run_subworkflow(self, workflow, subdirectory=None, from_scratch=None, **kwargs):
        '''
        Runs a workflow object, taking care of inheritance of several important properties

        '''

        # When testing, make sure the sub-workflow has access to the benchmark
        if hasattr(self, 'benchmark'):
            workflow.benchmark = self.benchmark

        # Ensure altering workflow.master_calcs won't affect self.master_calcs
        if workflow.master_calcs is self.master_calcs:
            workflow.master_calcs = copy.deepcopy(self.master_calcs)

        # Setting from_scratch to a non-None value will override the value of subworkflow.from_scratch...
        if from_scratch is None:
            workflow.from_scratch = self.from_scratch
        else:
            workflow.from_scratch = from_scratch

        # Link the list of calculations
        workflow.all_calcs = self.all_calcs

        if subdirectory is not None:
            # Update directories
            for key in workflow.master_calcs.keys():
                calc = workflow.master_calcs[key]
                for setting in calc._settings_that_are_paths:
                    path = getattr(calc, setting, None)
                    if path is not None and path.startswith(os.getcwd()):
                        new_path = os.path.abspath('./' + subdirectory + '/' + os.path.relpath(path))
                        setattr(calc, setting, new_path)

            # Run the workflow
            with utils.chdir(subdirectory):
                workflow.run(**kwargs)
        else:
            workflow.run(**kwargs)

        # ... and will prevent inheritance of from_scratch
        if from_scratch is None:
            self.from_scratch = workflow.from_scratch

    def todict(self):
        dct = self.__dict__
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
