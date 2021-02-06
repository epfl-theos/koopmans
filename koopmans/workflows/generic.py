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

valid_settings = [
    utils.Setting('task',
                  'Task to perform',
                  str, 'singlepoint', ('singlepoint', 'convergence', 'environ_dscf', 'ui')),
    utils.Setting('functional',
                  'Orbital-density-dependent-functional/density-functional to use',
                  str, 'ki', ('ki', 'kipz', 'pkipz', 'pbe', 'all')),
    utils.Setting('init_density',
                  'the functional to use to initialise the density',
                  str, 'pbe', ('pbe', 'pz', 'ki')),
    utils.Setting('init_variational_orbitals',
                  'which orbitals to use as an initial guess for the variational orbitals',
                  str, 'pz', ('pz', 'mlwfs', 'projw', 'ki', 'skip')),
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

    def __init__(self, workflow_settings, calcs_dct):
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
                if keyword not in calc._recognised_keywords:
                    raise ValueError(f'{keyword} is not a valid setting name')
                setattr(calc, keyword, value)
            return calc
        else:
            raise ValueError('Could not find a calculator of type {calc_type}')

    def run(self):
        raise NotImplementedError('This workflow class has not implemented the run() function')

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

                if self.from_scratch:
                    # Copy over nspin=2 wavefunction to nspin=1 tmp directory (if it has not been done already)
                    nspin1_tmpdir = f'{qe_calc.outdir}/{qe_calc.prefix}_{qe_calc.ndw}.save/K00001'
                    utils.system_call(
                        f'convert_nspin2_wavefunction_to_nspin1.sh {nspin2_tmpdir} {nspin1_tmpdir}')

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
            nspin2_tmpdir = f'{qe_calc.outdir}/{qe_calc.prefix}_{qe_calc.ndw}.save/K00001'

            # Copy over nspin=1 wavefunction to nspin=2 tmp directory
            if self.from_scratch:
                utils.system_call(
                    f'convert_nspin1_wavefunction_to_nspin2.sh {nspin1_tmpdir} {nspin2_tmpdir}')

            # PBE with nspin=2, reading in the spin-symmetric nspin=1 wavefunction
            master_qe_calc.name += '_nspin2'
            master_qe_calc.restart_mode = 'restart'
            master_qe_calc.ndr = 99
            self.run_calculator_single(master_qe_calc)

        else:

            self.run_calculator_single(master_qe_calc)

        return

    def run_calculator_single(self, qe_calc):
        # Runs qe_calc.calculate with additional options:

        ext_in = qe_calc.ext_in
        ext_out = qe_calc.ext_out

        # If an output file already exists, check if the run completed successfully
        verb = 'Running'
        if not self.from_scratch:

            calc_file = f'{qe_calc.directory}/{qe_calc.name}'

            if os.path.isfile(calc_file + ext_out):
                verb = 'Rerunning'

                # Load the old calc_file
                old_calc = qe_calc.__class__(qe_files=calc_file)

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
                    self.print_qc_keyval(qe_calc.name + '_' + result, val)

        # If we reached here, all future calculations should be performed from scratch
        self.from_scratch = True

        return

    def print_qc_keyval(self, key, value):
        '''
        Prints out a quality control message for testcode to evaluate
        '''
        print(f'<QC> {key} {value}')

    def run_subworkflow(self, workflow, from_scratch=None, **kwargs):
        '''
        Runs a workflow object, taking care of inheritance of several important properties

        Setting from_scratch to a non-None value will override the value of subworkflow.from_scratch
        and will prevent inheritance of from_scratch
        '''

        if from_scratch is None:
            workflow.from_scratch = self.from_scratch
        else:
            workflow.from_scratch = from_scratch

        workflow.run(**kwargs)

        self.all_calcs += workflow.all_calcs

        if from_scratch is None:
            self.from_scratch = workflow.from_scratch
