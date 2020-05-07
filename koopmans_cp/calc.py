"""

Calculator module for python_KI

Written by Edward Linscott Jan 2020

"""

from ase.io import espresso_cp as cp_io
from ase.calculators.espresso_cp import Espresso_cp
from koopmans_cp.io import cpi_diff, read_alpharef, read_pseudopotential, warn
import os
import sys

class Extended_Espresso_cp(Espresso_cp):

    '''
    An extension of the ASE Espresso_cp calculator class where we can access 
    the input parameters as attributes.

    e.g. 
    >>> calc.results
    will return calc.results
    >>> calc.nkscalfact
    will return calc.parameters['input_data']['nksic']['nkscalfact']

    I have not incorporated this directly into the ASE Espresso_cp module as 
    I imagine they won't be comfortable with overriding the __getattr__ and 
    __setattr__ methods
    '''

    def __init__(self, calc=None):
        # Optionally initialise using an existing ASE calculator
        if calc is None:
            super().__init__()
        else:
            super().__init__(label=getattr(calc, 'label', None), atoms=calc.atoms)
            self.results = calc.results
            self.parameters = calc.parameters
            if hasattr(calc, 'command'):
               self.command = calc.command
            if hasattr(calc, 'calc'):
               self.calc = calc.calc

    def __getattr__(self, name: str):
        '''
        When getting an attribute, also search self.parameters['input_data']
        Important notes:
          1. this routine is only called if self.__getattribute__ fails
          2. this routine will return 'None' rather than returning an AttributeError
             if the calculator does not possess the attribute AND if 'name' is a 
             recognised QE keyword
        '''

        parameters = self.__dict__.get('parameters', None)

        # Return calc.parameters['input_data'][block][name] if it exists
        if parameters is not None and 'input_data' in parameters:
            flattened_input_data = {k: v for block in parameters['input_data'].values()
                                    for k, v in block.items()}
            if name in flattened_input_data:
                return flattened_input_data[name]

        # Else return None if name is a recognised QE keyword
        for keywords in cp_io.KEYS.values():
            if name in keywords:
                return None

        # Else raise an AttributeError
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'")

    def __setattr__(self, name, value):

        # Set as self.parameters['input_data'] ahead of self.__dict__

        parameters = self.__dict__.get('parameters', None)

        if parameters is not None and 'input_data' in parameters:
            for block, keywords in cp_io.KEYS.items():
                if name in keywords:

                    # print(f'Setting {name} as a QE variable')
                    self.parameters['input_data'][block][name] = value

                    # If keyword is also an attribute (e.g. as in the case of 'prefix')
                    # set that too.
                    if name in self.__dict__:
                        self.__dict__[name] = value
                    return

        # print(f'Setting {name} as attribute')
        self.__dict__[name] = value

    def setattr_only(self, name, value):
        # If 'name' is both a QE parameter and an attribute of calc,
        # modify the attribute and leave the QE parameter unchanged

        parameters = self.__dict__.get('parameters', None)

        if parameters is not None and 'input_data' in parameters:
            for block, keywords in cp_io.KEYS.items():
                if name in keywords:
                    old_value = parameters['input_data'][block][name]
                    break

        setattr(self, name, value)
        parameters['input_data'][block][name] = old_value


def run_cp(calc, silent=True, from_scratch=False):
    '''
    Runs calc.calculate with additional options:

        silent :       if False, print status
        from_scratch : if False, check for a pre-existing calculation, and see 
                       if it has completed. If so, skip this calculation.

    returns True if the calculation was run, and False if the calculation was skipped
    '''

    # If an output file already exists, check if the run completed successfully
    if not from_scratch:
        calc_file = f'{calc.directory}/{calc.prefix}.cpo'
        if os.path.isfile(calc_file):
            old_cpo = next(cp_io.read_espresso_cp_out(calc_file)).calc

            if old_cpo.results['job_done']:
                # If it did, load the results, and exit
                calc.results = old_cpo.results
                if not silent:
                    print(
                        f'Not running {calc_file} as it is already complete')
                return False
            else:
                # If not, compare our settings with the settings of the preexisting .cpi file
                old_cpi = Extended_Espresso_cp(cp_io.read_espresso_cp_in(
                    calc_file.replace('cpo', 'cpi')).calc)
                diffs = cpi_diff([calc, old_cpi])
                if not silent:
                    print(f'Rerunning {calc_file}')
                if len(diffs) > 0:
                    for d in diffs:
                        old_value = getattr(old_cpi, d)
                        setattr(calc, d, old_value)
                        if not silent:
                            print(f'    Resetting {d}')

    if not silent:
        print('Running {}/{}...'.format(calc.directory,
                                        calc.prefix), end='', flush=True)

    calc.calculate()

    if not calc.results['job_done']:
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
