"""

Calculator module for python_KI

Written by Edward Linscott Jan 2020

Major modifications
May 2020: replaced Extended_Espresso_cp with QE_calc, a calculator class agnostic to the underlying ASE machinery
Sep 2020: moved individual calculators into calculators/

"""

import ase.io as ase_io
from ase.io import espresso_cp as cp_io
from ase.build import make_supercell
from koopmans.io import cpi_diff, read_alpharef, write_alpharef, warn
import os
import sys
import copy
import numpy as np
import koopmans.utils as utils


class QE_calc:

    '''
    A quantum espresso calculator object that...
     - stores cp.x/pw.x keywords in self._settings, but can be accessed like direct attributes
       e.g. self.<keyword> will return self._settings[<keyword>]
     - runs a cp.x/pw.x calculation upon self.calculate()
     - the calculation input/output files are 'self.directory/self.name.(cp/pw)(i/o)'
     - stores the results of this calculation in self.results

    Under the hood. it uses ASE to manage I/O and executing calculations.
      self.calc -> self._ase_calc
      self.directory -> self._ase_calc.directory
      self.name -> self._ase_calc.prefix
      self.results -> self._ase_calc.results
    This could be changed in the future without affecting the rest of the code

    From this generic class we will later define
      CP_calc for calculations using cp.x
      PW_calc for calculations using pw.x
      ... and others as listed in calculators/
    These are differentiated by self._io = cp_io/pw_io/...

    Arguments:
        calc       an ASE calculator object to initialize the QE calculation settings
        qe_files   an alternative to calc, initialize the settings using qe input file(s)
        kwargs     any valid quantum espresso keywords to alter
    '''

    def __init__(self, calc=None, qe_files=[], **kwargs):

        # If qe_input/output_files are provided, use them instead of calc
        if len(qe_files) > 0 and calc is not None:
            raise ValueError(
                f'Please only provide either the calc or qe_files argument to {self.__class__}')

        # Interpreting the qe_files argument
        if isinstance(qe_files, str):
            # If no extension is provided, automatically read both input and output files
            if '.' not in qe_files.split('/')[-1]:
                qe_files = [qe_files + self.ext_in, qe_files + self.ext_out]
            else:
                qe_files = [qe_files]

        # Read qe input file
        for qe_file in [f for f in qe_files if self.ext_in in f]:
            calc = ase_io.read(qe_file).calc

        # Read qe output file
        for qe_file in [f for f in qe_files if self.ext_out in f]:
            if calc is None:
                calc = ase_io.read(qe_file).calc
            else:
                # Copy over the results
                calc.results = ase_io.read(qe_file).calc.results

        # Initialise the calculator object
        if isinstance(calc, QE_calc):
            self._ase_calc = calc._ase_calc
        else:
            self._ase_calc = calc

        # Initialise a dictionary to store QE settings in
        self._settings = {}

        # Copy over settings from calc object into self._settings
        self._update_settings_dict()

        # Handle any recognised QE keywords passed as arguments
        for key, val in kwargs.items():
            if key not in self._recognised_keywords:
                raise ValueError(
                    f'{key} is not a recognised Quantum Espresso keyword')
            setattr(self, key, val)

    # By default, use cp.x
    _io = cp_io

    # extensions for i/o files
    ext_out = ''
    ext_in = ''

    @property
    def calc(self):
        # First, update the param block
        self._ase_calc.parameters['input_data'] = self.construct_namelist()

        return self._ase_calc

    @calc.setter
    def calc(self, value):
        self._ase_calc = value

    @property
    def directory(self):
        return self._ase_calc.directory

    @directory.setter
    def directory(self, value):
        self._ase_calc.directory = value

    @property
    def name(self):
        return self._ase_calc.prefix

    @name.setter
    def name(self, value):
        self._ase_calc.prefix = value

    @property
    def results(self):
        return self._ase_calc.results

    @results.setter
    def results(self, value):
        self._ase_calc.results = value

    @property
    def preprocessing_flags(self):
        return self._preprocessing_flags

    @preprocessing_flags.setter
    def preprocessing_flags(self, value):
        if not isinstance(value, list):
            value = [value]
        self._preprocessing_flags = value
        before, after = self._ase_calc.command.split(' ', 1)
        self._ase_calc.command = ' '.join([before] + value + [after])

    def calculate(self):
        # Generic function for running a calculation

        # If pseudo_dir is a relative path then make sure it accounts for self.directory
        if self.pseudo_dir is not None and self.pseudo_dir[0] != '/':
            directory_depth = self.directory.strip('./').count('/') + 1
            self.pseudo_dir = '../'*directory_depth + self.pseudo_dir

        self._ase_calculate()

    def _ase_calculate(self):
        # Update the keywords for the ASE calculator object
        self._ase_calc.parameters['input_data'] = self.construct_namelist()

        # Perform the calculation
        self._ase_calc.calculate()

    def _update_settings_dict(self):
        # Updates self._settings based on self._ase_calc
        self._settings = {}
        for namelist in self._ase_calc.parameters.get('input_data', {}).values():
            for key, val in namelist.items():
                self._settings[key] = val

    def construct_namelist(self):
        # Returns a namelist of settings, grouped by their Quantum Espresso headings
        return self._io.construct_namelist(**self._settings, warn=True)

    def parse_algebraic_setting(self, expr):
        # Checks if self._settings[expr] is defined algebraically, and evaluates them
        if not isinstance(expr, str):
            return expr
        if all([c.isalpha() or c in ['_'] for c in expr]):
            return expr

        expr = expr.replace('/', ' / ').replace('*', ' * ').split()
        for i, term in enumerate(expr):
            if term in ['*', '/']:
                continue
            elif all([c.isalpha() for c in term]):
                if self._settings.get(term, None) is None:
                    raise ValueError('Failed to parse ' + ''.join(expr))
                else:
                    expr[i] = self._settings[term]
            else:
                expr[i] = float(term)

        value = float(expr[0])
        for op, term in zip(expr[1::2], expr[2::2]):
            if op == '*':
                value *= float(term)
            elif op == '/':
                value /= float(term)
            else:
                raise ValueError('Failed to parse ' +
                                 ''.join([str(e) for e in expr]))

        return value

    def parse_algebraic_settings(self):
        # Checks self._settings for keywords defined algebraically, and evaluates them
        for key, value in self._settings.items():
            if key in ['pseudo_dir', 'outdir']:
                continue
            self._settings[key] = self.parse_algebraic_setting(value)

    def is_converged(self):
        raise ValueError(
            f'is_converged() function has not been implemented for {self.__class__}')

    def is_complete(self):
        raise ValueError(
            f'is_complete() function has not been implemented for {self.__class__}')

    def transform_to_supercell(self, matrix, **kwargs):
        # Converts to a supercell as given by a 3x3 transformation matrix
        assert np.shape(matrix) == (3, 3)
        self._ase_calc.atoms = make_supercell(self._ase_calc.atoms, matrix, **kwargs)
        self._ase_calc.atoms.calc = self._ase_calc

def run_qe(master_qe_calc, silent=True, from_scratch=False, enforce_ss=False):
    '''
    Wrapper for run_qe_single that manages the optional enforcing of spin symmetry
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
            if qe_calc.do_orbdep and qe_calc.odd_nkscalfact:
                # Rewrite alpha file for nspin=1
                alpha = read_alpharef(qe_calc)
                write_alpharef(alpha, qe_calc)
            from_scratch = run_qe_single(qe_calc, silent=silent, from_scratch=from_scratch)

            if from_scratch:
               # Copy over nspin=2 wavefunction to nspin=1 tmp directory (if it has not been done already)
               nspin1_tmpdir = f'{qe_calc.outdir}/{qe_calc.prefix}_{qe_calc.ndw}.save/K00001'
               utils.system_call(f'convert_nspin2_wavefunction_to_nspin1.sh {nspin2_tmpdir} {nspin1_tmpdir}')

        # PBE with nspin=1
        qe_calc = copy.deepcopy(master_qe_calc)
        qe_calc.name += '_nspin1'
        qe_calc.nspin, qe_calc.nelup, qe_calc.neldw, qe_calc.tot_magnetization = 1, None, None, None
        qe_calc.ndw, qe_calc.ndr = 98, 98
        nspin1_tmpdir = f'{qe_calc.outdir}/{qe_calc.prefix}_{qe_calc.ndw}.save/K00001'
        if qe_calc.do_orbdep and qe_calc.odd_nkscalfact and master_qe_calc.restart_mode == 'from_scratch':
            # Rewrite alpha file for nspin=1
            alpha = read_alpharef(qe_calc)
            write_alpharef(alpha, qe_calc)
        from_scratch = run_qe_single(qe_calc, silent=silent, from_scratch=from_scratch)

        # PBE from scratch with nspin=2 (dummy run for creating files of appropriate size)
        qe_calc = copy.deepcopy(master_qe_calc)
        qe_calc.name += '_nspin2_dummy'
        qe_calc.restart_mode = 'from_scratch'
        qe_calc.do_outerloop = False
        qe_calc.do_outerloop_empty = False
        qe_calc.ndw = 99
        if qe_calc.do_orbdep and qe_calc.odd_nkscalfact:
            # Rewrite alpha file for nspin=2
            alpha = read_alpharef(qe_calc, duplicated=False)
            write_alpharef(alpha, qe_calc)
        from_scratch = run_qe_single(qe_calc, silent=silent, from_scratch=from_scratch)
        nspin2_tmpdir = f'{qe_calc.outdir}/{qe_calc.prefix}_{qe_calc.ndw}.save/K00001'

        # Copy over nspin=1 wavefunction to nspin=2 tmp directory
        if from_scratch:
            utils.system_call(
                f'convert_nspin1_wavefunction_to_nspin2.sh {nspin1_tmpdir} {nspin2_tmpdir}')

        # PBE with nspin=2, reading in the spin-symmetric nspin=1 wavefunction
        master_qe_calc.name += '_nspin2'
        master_qe_calc.restart_mode = 'restart'
        master_qe_calc.ndr = 99
        return run_qe_single(master_qe_calc, silent=silent, from_scratch=from_scratch)

    else:
        return run_qe_single(master_qe_calc, silent, from_scratch)

def run_qe_single(qe_calc, silent=True, from_scratch=False):
    '''
    Runs qe_calc.calculate with additional options:

        silent :       if False, print status
        from_scratch : if False, check for a pre-existing calculation, and see
                       if it has completed. If so, skip this calculation.

    returns True if the calculation was run, and False if the calculation was skipped
    '''

    ext_in = qe_calc.ext_in
    ext_out = qe_calc.ext_out

    # If an output file already exists, check if the run completed successfully
    if not from_scratch:
        calc_file = f'{qe_calc.directory}/{qe_calc.name}'
        if os.path.isfile(calc_file + ext_out):
            old_calc = qe_calc.__class__(qe_files=calc_file)

            if old_calc.is_complete():
                # If it did, load the results, and exit
                qe_calc.results = old_calc.results
                if not silent:
                    print(
                        f'Not running {calc_file} as it is already complete')
                return False
            else:
                # If not, compare our settings with the settings of the preexisting .cpi file
                if ext_out == '.cpo':
                    diffs = cpi_diff([qe_calc, old_calc])
                else:
                    raise ValueError('pwi_diff needs to be implemented')
                if len(diffs) > 0:
                    for d in diffs:
                        old_value = getattr(old_calc, d, None)
                        setattr(qe_calc, d, old_value)
                        if not silent:
                            print(f'    Resetting {d}')

    if not silent:
        print('Running {}/{}...'.format(qe_calc.directory,
                                        qe_calc.name), end='', flush=True)

    qe_calc.calculate()

    if not qe_calc.is_complete():
        sys.exit(1)

    # Check spin-up and spin-down eigenvalues match
    if qe_calc.is_converged() and qe_calc.do_outerloop and qe_calc.nspin == 2 \
            and qe_calc.tot_magnetization == 0 and not qe_calc.fixed_state:
        rms_eigenval_difference = np.sqrt(
            np.mean(np.diff(qe_calc.results['eigenvalues'], axis=0)**2))
        if rms_eigenval_difference > 0.05:
            warn('Spin-up and spin-down eigenvalues differ substantially')

    if not silent:
        print(' done')

    return True


def calculate_alpha(calcs, filled=True, kipz=False):
    '''
    Calculates alpha via equation 10 of Nguyen et. al (2018) 10.1103/PhysRevX.8.021051
    If the band is filled, use s = 1; if the band is empty, use s = 0

    Arguments:
        calcs  -- a list of calculations
        filled -- True if the orbital for which we're calculating alpha is filled
        kipz   -- True if KIPZ; False if KI
    '''

    if kipz and filled:
        # KIPZ N
        [alpha_calc] = [c for c in calcs if c.which_orbdep == 'nkipz'
                        and c.do_orbdep and c.f_cutoff == 1.0]
        kipz = alpha_calc.results

        # KIPZ N-1
        [kipz_m1] = [c.results for c in calcs if c.which_orbdep == 'nkipz'
                     and c.do_orbdep and c.f_cutoff < 0.0001]

        # PBE N
        [pbe] = [c.results for c in calcs if not c.do_orbdep
                 and c.restart_mode == 'restart' and c.f_cutoff == 1.0]

        dE = kipz['energy'] - kipz_m1['energy']
        lambda_a = kipz['lambda_ii']
        lambda_0 = pbe['lambda_ii']

    elif kipz and not filled:
        # KIPZ N+1-1
        [alpha_calc] = [c for c in calcs if c.which_orbdep == 'nkipz'
                        and c.do_orbdep and c.f_cutoff < 0.0001]
        kipz = alpha_calc.results

        # KIPZ N+1
        [kipz_p1] = [c.results for c in calcs if c.which_orbdep == 'nkipz'
                     and c.do_orbdep and c.f_cutoff == 1.0]

        # PBE N+1-1
        [pbe] = [c.results for c in calcs if not c.do_orbdep
                 and c.restart_mode == 'restart' and c.f_cutoff < 0.0001]

        dE = kipz_p1['energy'] - kipz['energy']
        lambda_a = kipz['lambda_ii']
        lambda_0 = pbe['lambda_ii']

    elif not kipz and filled:
        # KI N
        [alpha_calc] = [c for c in calcs if c.which_orbdep == 'nki'
                        and c.do_orbdep and c.f_cutoff == 1.0]
        ki = alpha_calc.results

        # PBE N-1
        [pbe_m1] = [c.results for c in calcs if not c.do_orbdep
                    and c.restart_mode == 'restart' and c.f_cutoff < 0.0001]

        # PBE N
        [pbe] = [c.results for c in calcs if not c.do_orbdep
                 and c.restart_mode == 'restart' and c.f_cutoff == 1.0]

        dE = pbe['energy'] - pbe_m1['energy']
        lambda_a = ki['lambda_ii']
        lambda_0 = pbe['lambda_ii']

    else:
        # KI N+1-1
        [alpha_calc] = [c for c in calcs if c.which_orbdep == 'nki'
                        and c.do_orbdep and c.f_cutoff < 0.0001]
        ki = alpha_calc.results

        # PBE N+1
        [pbe_p1] = [c.results for c in calcs if not c.do_orbdep
                    and c.restart_mode == 'restart' and c.f_cutoff == 1.0]

        # PBE N+1-1
        [pbe] = [c.results for c in calcs if not c.do_orbdep
                 and c.restart_mode == 'restart' and c.f_cutoff < 0.0001]

        dE = pbe_p1['energy'] - pbe['energy']
        lambda_a = ki['lambda_ii']
        lambda_0 = pbe['lambda_ii']

    if not kipz:
        # Check total energy is unchanged
        dE_check = abs(ki['energy'] - pbe['energy'])
        if dE_check > 1e-5:
            warn('KI and PBE energies differ by {:.5f} eV'.format(dE_check))

    # Obtaining alpha
    if (alpha_calc.odd_nkscalfact and filled) or (alpha_calc.odd_nkscalfact_empty and not filled):
        alpha_guesses = read_alpharef(alpha_calc)
        alpha_guess = alpha_guesses[alpha_calc.fixed_band - 1]
    else:
        alpha_guess = alpha_calc.nkscalfact

    alpha = alpha_guess*(dE - lambda_0) / (lambda_a - lambda_0)

    # The error is lambda^alpha(1) - lambda^alpha_i(1)
    error = dE - lambda_a

    return alpha, error
