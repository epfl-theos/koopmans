"""

Calculator module for python_KI

Written by Edward Linscott Jan 2020

Major modifications
May 2020: replaced Extended_Espresso_cp with CP_calc, a calculator class agnostic to the underlying ASE machinery

"""

from ase.io import espresso_cp as cp_io
from ase.calculators.espresso_cp import Espresso_cp
from koopmans_cp.io import cpi_diff, read_alpharef, warn
import os
import sys
import copy

class CP_calc:

    '''
    A cp.x calculator object that...
     - stores cp.x keywords in self._settings, but can be accessed like direct attributes
       e.g. self.<keyword> will return self._settings[<keyword>]
     - runs a cp.x calculation upon self.calculate()
     - the calculation input/output files are 'self.directory/self.name.cp(i/o)'
     - stores the results of this calculation in self.results

    Under the hood. it uses ASE to manage I/O and executing calculations.
      self.calc -> self._ase_calc
      self.directory -> self._ase_calc.directory
      self.name -> self._ase_calc.prefix
      self.results -> self._ase_calc.results
    This could be changed in the future without affecting the rest of the code
    '''

    def __init__(self, calc, **kwargs):

        # Initialise the calculator object
        self._ase_calc = calc

        # Initialise a dictionary to store cp.x settings in
        self._settings = {}

        # Copy over settings from calc object into self._settings
        self._update_settings_dict()

        # Handle any recognised cp.x keywords passed as arguments
        for key, val in kwargs.items():
            if key not in self._recognised_cp_keywords:
                raise ValueError(f'{key} is not a recognised cp.x keyword')
            setattr(self, key, val)

    @property
    def calc(self):
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

    # Adding all cp.x keywords as decorated properties of the CP_calc class.
    # This means one can set and get cp.x keywords as self.<keyword> but
    # internally they are stored as self._settings['keyword'] rather than 
    # self.<keyword>
    _recognised_cp_keywords = []

    for keywords in cp_io.KEYS.values():
        for k in keywords: 
            _recognised_cp_keywords.append(k)

            # We need to use these make_get/set functions so that get/set_k are
            # evaluated immediately (otherwise we run into late binding and 'k'
            # is not defined when get/set_k are called)
            def make_get(key):
                def get_k(self):
                    # Return 'None' rather than an error if the keyword has not
                    # been defined
                    return self._settings.get(key, None)
                return get_k

            def make_set(key):
                def set_k(self, value):
                    self._settings[key] = value
                return set_k

            get_k = make_get(k)
            set_k = make_set(k)
            locals()[k] = property(get_k, set_k)     

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
        for namelist in self._ase_calc.parameters['input_data'].values():
            for key, val in namelist.items():
                self._settings[key] = val

    def construct_namelist(self):
        # Returns a namelist of settings, grouped by their cp.x headings
        return cp_io.construct_namelist(**self._settings)

    def parse_algebraic_setting(self, expr):
        # Checks if self._settings[expr] is defined algebraically, and evaluates them
        if not isinstance(expr, str):
            return expr
        if all([c.isalpha() for c in expr]):
            return expr

        expr = expr.replace('/', ' / ').replace('*', ' * ').split()
        for i, term in enumerate(expr):
            if term in ['*', '/']:
                continue
            elif all([c.isalpha() for c in term]):
                if self._settings.get(term, None) is None:
                    raise ValueError(f'Failed to parse ' + ''.join(expr))
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

def run_cp(cp_calc, silent=True, from_scratch=False):
    '''
    Runs cp_calc.calculate with additional options:

        silent :       if False, print status
        from_scratch : if False, check for a pre-existing calculation, and see 
                       if it has completed. If so, skip this calculation.

    returns True if the calculation was run, and False if the calculation was skipped
    '''

    # If an output file already exists, check if the run completed successfully
    if not from_scratch:
        calc_file = f'{cp_calc.directory}/{cp_calc.name}.cpo'
        if os.path.isfile(calc_file):
            old_cpo = next(cp_io.read_espresso_cp_out(calc_file)).calc

            if old_cpo.results['job_done']:
                # If it did, load the results, and exit
                cp_calc.results = old_cpo.results
                if not silent:
                    print(
                        f'Not running {calc_file} as it is already complete')
                return False
            else:
                # If not, compare our settings with the settings of the preexisting .cpi file
                old_cpi = CP_calc(cp_io.read_espresso_cp_in(calc_file.replace('cpo', 'cpi')).calc)
                diffs = cpi_diff([cp_calc, old_cpi])
                if not silent:
                    print(f'Rerunning {calc_file}')
                if len(diffs) > 0:
                    for d in diffs:
                        old_value = getattr(old_cpi, d)
                        setattr(cp_calc, d, old_value)
                        if not silent:
                            print(f'    Resetting {d}')

    if not silent:
        print('Running {}/{}...'.format(cp_calc.directory,
                                        cp_calc.name), end='', flush=True)

    cp_calc.calculate()

    if not cp_calc.results['job_done']:
        sys.exit(1)

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
