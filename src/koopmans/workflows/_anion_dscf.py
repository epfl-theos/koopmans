import os

from ase.calculators.calculator import CalculationFailed

from koopmans import utils
from koopmans.calculators import EnvironCalculator

from ._workflow import Workflow

'''
Workflow for performing delta SCF PBE calculations using pw.x --environ
'''


class DeltaSCFWorkflow(Workflow):

    def _run(self):
        # Run workflow
        if self.parameters.from_scratch:
            utils.system_call('rm -r neutral charged 2> /dev/null', False)

        epsilons = sorted(self.parameters.eps_cavity, reverse=True)
        self.print('PBE ΔSCF WORKFLOW', style='heading')

        # Remove settings from calculator_parameters that will be overwritten
        for key in ['tot_charge', 'tot_magnetization', 'disk_io', 'restart_mode']:
            self.calculator_parameters['pw'].pop(key, None)

        for charge, label in zip([0, -1], ['neutral', 'charged']):
            self.print(f'Performing {label} calculations', style='subheading')

            # Initialize variables
            i_eps = 0
            epsilon = epsilons[i_eps]
            environ_restart = False
            restart_mode = 'from_scratch'
            calc_succeeded = True

            # Create working directories
            if not os.path.isdir(label):
                utils.system_call(f'mkdir {label}')
                utils.system_call(f'mkdir {label}/{epsilons[0]}')

            while calc_succeeded:
                # Create a new Environ calculator object
                pw_calc = EnvironCalculator(atoms=self.atoms, restart_mode=restart_mode, disk_io='medium',
                                            tot_charge=charge, tot_magnetization=-charge,
                                            **self.calculator_parameters['pw'])
                pw_calc.prefix = 'pbe'
                pw_calc.directory = f'{label}/{epsilon}'
                pw_calc.parameters.outdir = 'TMP'

                # Update the environ settings
                pw_calc.environ_settings['ENVIRON']['environ_restart'] = environ_restart
                pw_calc.environ_settings['ENVIRON']['env_static_permittivity'] = epsilon
                pw_calc.environ_settings['BOUNDARY']['solvent_mode'] = 'ionic'
                pw_calc.environ_settings['ELECTROSTATIC']['tol'] = 1e-8

                # self.from_scratch = True means that run_qe won't try and skip this calculation
                # if it encounters pre-existing QE output files, and NOT that QE will use
                # restart_mode = 'from_scratch'
                self.parameters.from_scratch = True

                # Run the calculator
                try:
                    self.run_calculator(pw_calc)
                except CalculationFailed:
                    self.print(' failed to converge (which is expected for small values of epsilon)')
                    return
                calc_succeeded = pw_calc.is_converged()

                # Preparing for next loop
                i_eps += 1
                if epsilon == 1 or i_eps == len(epsilons):
                    break
                new_epsilon = max(epsilons[i_eps], 1)
                outdir = os.path.relpath(pw_calc.parameters.outdir, pw_calc.directory)
                utils.system_call(f'rsync -a {label}/{epsilon}/{outdir} {label}/{new_epsilon}/')
                epsilon = new_epsilon
                restart_mode = 'restart'
                environ_restart = True
