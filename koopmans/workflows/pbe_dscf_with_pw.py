import os
from ase.calculators.calculator import CalculationFailed
from koopmans import utils, io
from koopmans.calculators.calculator import run_qe
from koopmans.calculators.environ import Environ_calc

'''
Workflow for performing delta SCF PBE calculations using pw.x --environ
'''


def run(workflow_settings, calcs_dct):

    from koopmans.config import from_scratch

    if 'pw' not in calcs_dct:
        raise ValueError(
            'You need to provide a pw block in your input .json file for task = environ_dscf')

    pw_calc = Environ_calc(calc=calcs_dct['pw'])

    # Run workflow
    calc_succeeded = True
    pw_calc.name = 'pbe'

    if from_scratch:
        utils.system_call('rm -r neutral charged 2> /dev/null', False)

    epsilons = sorted(workflow_settings['eps_cavity'], reverse=True)
    print('PBE ΔSCF WORKFLOW')

    for charge, label in zip([0, -1], ['neutral', 'charged']):
        print(f'\nPerforming {label} calculations...')

        # Initialize value of epsilon
        i_eps = 0
        epsilon = epsilons[i_eps]

        # Create working directories
        if not os.path.isdir(label):
            utils.system_call(f'mkdir {label}')
            utils.system_call(f'mkdir {label}/{epsilons[0]}')

        # Initialize system parameters
        pw_calc.restart_mode = 'from_scratch'
        pw_calc.disk_io = 'medium'  # checkpointing files will be required for later restarts
        pw_calc.environ_settings['ENVIRON']['environ_restart'] = False

        # Apply the desired charge
        pw_calc.tot_charge = charge
        pw_calc.tot_magnetization = -charge

        while calc_succeeded:
            pw_calc.directory = f'{label}/{epsilon}'
            pw_calc.environ_settings['ENVIRON']['env_static_permittivity'] = epsilon
            pw_calc.environ_settings['BOUNDARY']['solvent_mode'] = 'ionic'
            pw_calc.environ_settings['ELECTROSTATIC']['tol'] = 1e-8

            from_scratch = True
            # from_scratch = True means that run_qe won't try and skip this calculation
            # if it encounters pre-existing QE output files, and NOT that QE will use
            # restart_mode = 'from_scratch'

            try:
                run_qe(pw_calc, silent=False)
            except CalculationFailed:
                print(' failed to converge')
                print('\nWORKFLOW COMPLETE\n')
                return

            calc_succeeded = pw_calc.is_converged()

            # QC
            if workflow_settings['print_qc']:
                for key in ['energy', 'electrostatic embedding']:
                    value = pw_calc.results.get(key, None)
                    if value is not None:
                        io.print_qc(key, value)

            # Preparing for next loop
            i_eps += 1
            if epsilon == 1 or i_eps == len(epsilons):
                break
            new_epsilon = max(epsilons[i_eps], 1)
            utils.system_call(f'cp -r {label}/{epsilon} {label}/{new_epsilon}')
            epsilon = new_epsilon
            pw_calc.restart_mode = 'restart'

            pw_calc.environ_settings['ENVIRON']['environ_restart'] = True

    print('\nWORKFLOW COMPLETE\n')
