

"""

Workflow module for performing a single DFT calculation with either kcp.x or pw.x

Written by Edward Linscott Oct 2020

"""

import copy
import shutil
from pathlib import Path
from typing import TypeVar

from koopmans import calculators, pseudopotentials, utils

from ._workflow import Workflow

T = TypeVar('T', bound='calculators.CalculatorExt')


class DFTWorkflow(Workflow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.parameters.functional = 'dft'


class DFTCPWorkflow(DFTWorkflow):

    def _run(self):

        calc = self.new_calculator('kcp')

        # Removing old directories
        if self.parameters.from_scratch:
            utils.system_call(f'rm -r {calc.parameters.outdir} 2>/dev/null', False)

        calc.prefix = 'dft'
        calc.directory = '.'
        calc.parameters.ndr = 50
        calc.parameters.ndw = 51
        calc.parameters.restart_mode = 'from_scratch'
        calc.parameters.do_orbdep = False
        calc.parameters.fixed_state = False
        calc.parameters.do_outerloop = True
        calc.parameters.which_compensation = 'tcc'

        if calc.parameters.maxiter is None:
            calc.parameters.maxiter = 300
        if calc.has_empty_states():
            calc.parameters.do_outerloop_empty = True
            if calc.parameters.empty_states_maxstep is None:
                calc.parameters.empty_states_maxstep = 300

        self.run_calculator(calc, enforce_ss=self.parameters.fix_spin_contamination)

        return calc


class DFTPWWorkflow(DFTWorkflow):

    def _run(self):

        # Create the calculator
        calc = self.new_calculator('pw')

        # Update keywords
        calc.prefix = 'dft'
        calc.parameters.ndr = 50
        calc.parameters.ndw = 51
        calc.parameters.restart_mode = 'from_scratch'

        # Remove old directories
        if self.parameters.from_scratch:
            utils.system_call(f'rm -r {calc.parameters.outdir} 2>/dev/null', False)

        # Run the calculator
        self.run_calculator(calc)

        return


class DFTPhWorkflow(Workflow):

    def _run(self):

        self.print('Calculate the dielectric tensor', style='heading')

        calc_scf = self.new_calculator('pw', nbnd=None)
        calc_scf.prefix = 'scf'
        self.run_calculator(calc_scf)
        calc_ph = self.new_calculator('ph', epsil=True, fildyn=f'{self.name}.dynG')
        calc_ph.prefix = 'eps'
        self.run_calculator(calc_ph)


class DFTBandsWorkflow(DFTWorkflow):

    def _run(self):

        self.print('DFT bandstructure workflow', style='heading')

        if self.parameters.from_scratch:
            path = Path('dft_bands')
            if path.exists():
                shutil.rmtree(path)

        with utils.chdir('dft_bands'):
            # First, a scf calculation
            calc_scf = self.new_calculator('pw', nbnd=None)
            calc_scf.prefix = 'scf'
            self.run_calculator(calc_scf)

            # Second, a bands calculation
            calc_bands = self.new_calculator('pw', calculation='bands', kpts=self.kpoints.path)
            calc_bands.prefix = 'bands'
            self.run_calculator(calc_bands)

            # Prepare the band structure for plotting
            bs = calc_bands.results['band structure']

            # Third, a PDOS calculation
            pseudos = [pseudopotentials.read_pseudo_file(calc_scf.parameters.pseudo_dir / p) for p in
                       self.pseudopotentials.values()]
            if all([int(p['header'].get('number_of_wfc', 0)) > 0 for p in pseudos]):
                calc_dos = self.new_calculator('projwfc')
                calc_dos.pseudo_dir = calc_bands.parameters.pseudo_dir
                self.run_calculator(calc_dos)

                # Prepare the DOS for plotting
                dos = copy.deepcopy(calc_dos.results['dos'])
                dos._energies -= bs.reference
            else:
                # Skip if the pseudos don't have the requisite PP_PSWFC blocks
                utils.warn('Some of the pseudopotentials do not have PP_PSWFC blocks, which means a projected DOS '
                           'calculation is not possible. Skipping...')
                dos = None

            # Plot the band structure and DOS
            self.plot_bandstructure(bs.subtract_reference(), dos)

    def new_calculator(self,
                       calc_type: str,
                       *args,
                       **kwargs) -> T:   # type: ignore[type-var, misc]
        calc: T = super().new_calculator(calc_type, *args, **kwargs)
        if calc_type == 'projwfc':
            assert isinstance(calc, calculators.ProjwfcCalculator)
            calc.parameters.filpdos = self.name
            calc.pseudopotentials = self.pseudopotentials
            calc.spin_polarized = self.parameters.spin_polarized
        return calc
