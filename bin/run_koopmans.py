#!/usr/bin/env python3

import argparse
import textwrap
from koopmans import io, ase
from koopmans.workflows import workflow, pbe_dscf_with_pw
from koopmans import config
from collections import namedtuple

'''
Perform KI/KIPZ calculations
'''

Setting = namedtuple(
    'Setting', ['name', 'description', 'type', 'default', 'options'])

valid_settings = [
    Setting('task',
            'Task to perform',
            str, 'singlepoint', ('singlepoint', 'convergence', 'environ_dscf')),
    Setting('functional',
            'Orbital-density-dependent-functional/density-functional to use',
            str, 'ki', ('ki', 'kipz', 'pkipz', 'pbe', 'all')),
    Setting('init_density',
            'how to initialise the density',
            str, 'pbe', ('pbe', 'pz', 'ki')),
    Setting('init_manifold',
            'how to initialise the variational orbitals',
            str, 'pz', ('pz', 'ki', 'skip')),
    Setting('wannierize',
            'if True, a wannierization procedure precedes the rest of the workflow',
            bool, False, (True, False)),
    Setting('n_max_sc_steps',
            'maximum number of self-consistency steps for calculating alpha',
            int, 1, None),
    Setting('alpha_conv_thr',
            'convergence threshold for |delta E - lambda|; if below this '
            'threshold, the corresponding alpha value is not updated',
            (float, str), 1e-3, None),
    Setting('calculate_alpha',
            'if True, the screening parameters will be calculated; if False, '
            'they will be read directly from file',
            bool, True, (True, False)),
    Setting('alpha_guess',
            'starting guess for alpha (overridden if alpha_from_file is true)',
            float, 0.6, None),
    Setting('alpha_from_file',
            'if True, uses the file_alpharef.txt from the base directory as a '
            'starting guess',
            bool, False, (True, False)),
    Setting('print_qc',
            'if True, prints out strings for the purposes of quality control',
            bool, False, (True, False)),
    Setting('from_scratch',
            'if True, will delete any preexisting workflow and start again; '
            'if False, will resume a workflow from where it was last up to',
            bool, False, (True, False)),
    Setting('orbital_groups',
            'a list of integers the same length as the total number of bands, '
            'denoting which bands to assign the same screening parameter to',
            list, None, None),
    Setting('enforce_spin_symmetry',
            'if True, the spin-up and spin-down wavefunctions will be forced '
            'to be the same',
            bool, True, (True, False)),
    Setting('convergence_observable',
            'System observable of interest which we converge',
            str, 'total energy', ('homo energy', 'lumo energy', 'total energy')),
    Setting('convergence_threshold',
            'Convergence threshold for the system observable of interest',
            str, None, None),
    Setting('convergence_parameters',
            'The observable of interest will be converged with respect to this/these'
            'simulation parameter(s)',
            (list, str), ['ecutwfc'], None),
    Setting('eps_cavity',
            'a list of epsilon_infinity values for the cavity in dscf calculations',
            list, [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20], None)]

valid_settings_dict = {s.name: s for s in valid_settings}


def check_settings(settings):
    '''
    Checks workflow settings against the list of valid settings, populates
    missing keywords with their default values, and lowers any uppercases

    Arguments:
        settings -- a dictionary of workflow settings

    '''

    for key, value in settings.items():
        # Check key is a valid keyword
        if key in valid_settings_dict:
            valid_setting = valid_settings_dict[key]

            # Lowers any uppercase strings
            if isinstance(value, str):
                settings[key] = value.lower()
                value = value.lower()

            # Check value is the correct type
            if not isinstance(value, valid_setting.type) and value is not None:
                if isinstance(valid_setting.type, tuple):
                    raise ValueError(
                        f'{type(value).__name__} is an invalid type for "{key}" (must be '
                        'one of ' + '/'.join([t.__name__ for t in valid_setting.type]) + ')')
                else:
                    raise ValueError(
                        f'{type(value).__name__} is an invalid type for "{key}" (must be {valid_setting.type.__name__})')

            # Check value is among the valid options
            if valid_setting.options is not None and value not in valid_setting.options:
                raise ValueError(
                    f'"{value}" is an invalid value for "{key}" (options are {"/".join(valid_setting.options)})')
        else:
            raise ValueError(f'"{key}" is not a recognised workflow setting')

    # Populate missing settings with the default options
    for setting in valid_settings:
        if setting.name not in settings:
            settings[setting.name] = setting.default

    # Parse physicals
    for physical in ['alpha_conv_thr', 'convergence_threshold']:
        settings[physical] = io.parse_physical(settings[physical])

    # Set calculator.from_scratch
    config.init(from_scratch=workflow_settings['from_scratch'])


if __name__ == '__main__':

    # Automatically constructing a list of workflow keywords for 'run_koopmans.py --help'
    # from valid_settings
    epilog = ''
    maxlen = max([len(s.name) for s in valid_settings]) + 2
    for s in valid_settings:
        entry = f'  {s.name.ljust(maxlen)}{s.description} ('
        if isinstance(s.type, tuple):
            entry += '/'.join([t.__name__ for t in s.type])
        else:
            entry += s.type.__name__
        entry += f', default {s.default}'
        if s.options is not None and s.type is not bool:
            entry += ', must be ' + '/'.join([str(o) for o in s.options])
        entry += ')'
        for line in textwrap.wrap(entry, subsequent_indent=' ' * (maxlen + 2)):
            epilog += '\n' + line
        epilog += '\n'

    # Construct parser
    parser = argparse.ArgumentParser(
        description='Perform a KI/KIPZ calculation using Quantum Espresso',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="'workflow' block arguments:" + epilog)
    parser.add_argument('json', metavar='system.json', type=str,
                        help='a single JSON file containing the workflow and code settings')

    # Parse arguments
    args = parser.parse_args()

    # Reading in JSON file
    workflow_settings, calcs_dct = ase.read_json(args.json)

    # Check that the workflow settings are valid and populate missing
    # settings with default values
    check_settings(workflow_settings)

    if workflow_settings['task'] == 'convergence':
        workflow.run_convergence(workflow_settings, calcs_dct)
    elif workflow_settings['task'] == 'singlepoint':
        workflow.run_singlepoint(workflow_settings, calcs_dct)
    elif workflow_settings['task'] == 'environ_dscf':
        pbe_dscf_with_pw.run(workflow_settings, calcs_dct)
