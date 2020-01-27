"""

Calculator module for python_KI

Written by Edward Linscott Jan 2020

"""

from ase.io import espresso_cp as cp_io
from ase.calculators.espresso_cp import Espresso_cp
import os
import warnings
import copy

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
            self.command = calc.command
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
                # If it did, exit immediately
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
            warnings.warn(
                'KI and PBE energies differ by {:.5f} eV'.format(dE_check))

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


def set_up_calculator(calc, calc_type='pbe_init', **kwargs):
    """
    Generates a new ASE calculator based on a template ASE calculator, modifying
    the appropriate settings to match the chosen calc_type, and altering any 
    Quantum Espresso keywords specified as kwargs

    Arguments:

        calc: an template ASE calculator. N.B. this object is not modified by
              this routine

        calc_type: the calculation type; must be one of

            Initialisation
            'pbe_init'  PBE calculation from scratch
            'pz'        PZ calculation starting from PBE restart
            'kipz_init' KIPZ starting from PBE restart

            For calculating alpha_i for filled orbitals
            'pbe'      PBE calculation starting from restart
            'pbe_n-1'  PBE calculation with N-1 electrons via fixed_state
            'kipz_n-1' KIPZ calculation with N-1 electrons via fixed_state
            'ki'       KI calculation with N electrons
            'kipz'     KIPZ calculation with N electrons

            For calculating alpha_i for empty orbitals
            'pz_print' PZ calculation that generates evcempty_fixed.dat file
            'kipz_print'    KIPZ calculation that generates evcempty_fixed.dat file
            'pbe_n+1_dummy' PBE dummy calculation that generates store files of
                            the correct dimensions
            'pbe_n+1'       PBE calculation with N+1 electrons
            'kipz_n+1'      KIPZ calculation with N+1 electrons
            'pbe_n+1-1'     PBE calculation with N electrons, starting from N+1 
                            but with f_cutoff = 0.00001
            'ki_n+1-1'      KI calculation with N electrons, starting from N+1 
                            but with f_cutoff = 0.00001
            'kipz_n+1-1'    KIPZ calculation with N electrons, starting from N+1 
                            but with f_cutoff = 0.00001

       **kwargs: accepts any Quantum Espresso keywords as an argument, and will
                 apply these options to the returned ASE calculator

    Returns: a new ASE calculator object

    """

    # Avoid modifying the master calculator
    calc = copy.deepcopy(calc)

    # Set up read/write indexes
    if calc_type == 'pbe_init':
        ndr = 50
        ndw = 50
    elif calc_type in ['pz', 'kipz_init']:
        ndr = 50
        ndw = 51
    elif calc_type == 'pbe':
        ndr = 51
        ndw = 52
    elif calc_type == 'pbe_n-1':
        ndr = 51
        ndw = 53
    elif calc_type in ['pz_print', 'kipz_print']:
        ndr = 51
        ndw = 54
    elif calc_type == 'pbe_n+1_dummy':
        ndr = 55
        ndw = 55
    elif calc_type == 'pbe_n+1-1':
        ndr = 55
        ndw = 56
    elif calc_type == 'pbe_n+1':
        ndr = 55
        ndw = 57
    elif calc_type in ['ki', 'kipz']:
        ndr = 51
        ndw = 60
    elif calc_type == 'kipz_n-1':
        ndr = 51
        ndw = 70
    elif calc_type == 'kipz_n+1':
        ndr = 55
        ndw = 80
    elif calc_type in ['ki_n+1-1', 'kipz_n+1-1']:
        ndr = 55
        ndw = 90
    else:
        raise ValueError('Invalid calc_type "{}"'.format(calc_type))

    # Pseudopotentials

    # This script will...
    #  1. try to locating the directory as currently specified by the calculator
    #  2. if that fails, it will check if $ESPRESSO_PSEUDO is set
    #  3. if that fails, it will raise an OS error

    # Therefore, if you want to use a pseudo directory that is different to your
    # $PSEUDO_DIR, reset pseudo_dir for the template calc object before entering
    # this routine.

    if calc.pseudo_dir is None or not os.path.isdir(calc.pseudo_dir):
        try:
            calc.pseudo_dir = os.environ.get('ESPRESSO_PSEUDO')
        except:
            raise OSError('Directory for pseudopotentials not found. Please define '
                          'the environment variable ESPRESSO_PSEUDO or provide the function run_cp '
                          'with a pseudo_dir argument')

    # CP options
    # control
    calc.setattr_only('prefix', calc_type)
    calc.ndw = ndw
    calc.ndr = ndr
    if calc.prefix in ['pbe_init', 'pbe_n+1_dummy']:
        calc.restart_mode = 'from_scratch'
    else:
        calc.restart_mode = 'restart'

    # system
    calc.nspin = 2
    if 'pz' in calc.prefix or 'ki' in calc.prefix:
        calc.do_orbdep = True
    else:
        calc.do_orbdep = False
    if calc.prefix in ['pbe_init', 'pz_print', 'pbe_n+1_dummy', 'pz', 'kipz_init']:
        calc.fixed_state = False
    else:
        calc.fixed_state = True
        if '-1' in calc.prefix:
            calc.f_cutoff = 1e-5
        else:
            calc.f_cutoff = 1.0
    if 'n+1' in calc.prefix:
        calc.nelec += 1
        calc.nelup += 1
    if calc.prefix in ['pbe_n+1', 'pbe_n+1-1', 'ki_n+1-1', 'kipz_n+1-1', 'kipz_n+1']:
        calc.restart_from_wannier_pwscf = True

    # electrons
    calc.conv_thr = calc.nelec*1e-8
    calc.empty_states_maxstep = 300
    if ('pz' in calc.prefix or 'ki' in calc.prefix) and 'kipz' not in calc.prefix:
        calc.maxiter = 2
    else:
        calc.maxiter = 200

    # nksic
    calc.esic_conv_thr = calc.nelec*1e-8
    calc.do_innerloop_cg = True
    if calc.prefix[:2] == 'pz':
        calc.do_innerloop = True
        calc.one_innerloop_only = True
    else:
        calc.do_innerloop = False
        calc.one_innerloop_only = False
    if 'kipz' in calc.prefix:
        calc.which_orbdep = 'nkipz'
    elif calc.prefix == 'pz':
        calc.which_orbdep = 'pz'
    elif 'ki' in calc.prefix:
        calc.which_orbdep = 'nki'
    if 'print' in calc.prefix:
        calc.print_wfc_anion = True

    # Handle any keywords provided by kwargs
    # Note that since this is performed after the above logic this can (deliberately
    # or accidentally) overwrite the above settings

    for keyword, value in kwargs.items():
        setattr(calc, keyword, value)

    # Sanity checking
    if calc.print_wfc_anion and calc.index_empty_to_save is None:
        raise ValueError('Error: print_wfc_anion is set to true but you have not selected '
                         ' an index_empty_to_save. Provide this as an argument to set_calculator')

    if calc.fixed_band is not None and calc.fixed_band > calc.nelup + 1:
        warnings.warn(
            'calc.fixed_band is higher than the LUMO; this should not happen')

    return calc


