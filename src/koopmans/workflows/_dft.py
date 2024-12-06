

"""

Workflow module for performing a single DFT calculation with either kcp.x or pw.x

Written by Edward Linscott Oct 2020

"""

import copy
import shutil
from pathlib import Path
from typing import Generator, List, TypeVar

from koopmans import calculators, pseudopotentials, utils
from koopmans.outputs import OutputModel
from koopmans.status import Status
from koopmans.step import Step

from ._workflow import Workflow

T = TypeVar('T', bound='calculators.CalculatorExt')


class DFTWorkflow(Workflow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.parameters.functional = 'dft'


class DFTCPOutput(OutputModel):
    pass


class DFTCPWorkflow(DFTWorkflow):
    output_model = DFTCPOutput  # type: ignore
    outputs: DFTCPOutput

    def _run(self) -> None:

        calc = self.new_calculator('kcp')
        assert isinstance(calc, calculators.KoopmansCPCalculator)

        calc.prefix = 'dft'
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

        status = self.run_steps(calc)
        if status == Status.COMPLETED:
            self.status = Status.COMPLETED

        return


class DFTPWOutput(OutputModel):
    pass


class DFTPWWorkflow(DFTWorkflow):

    output_model = DFTPWOutput  # type: ignore
    outputs: DFTPWOutput

    def _run(self) -> None:

        # Create the calculator
        calc = self.new_calculator('pw')

        # Update keywords
        calc.prefix = 'dft'
        calc.parameters.ndr = 50
        calc.parameters.ndw = 51
        calc.parameters.restart_mode = 'from_scratch'

        status = self.run_steps(calc)
        if status == Status.COMPLETED:
            self.status = Status.COMPLETED

        return


class DFTPhOutput(OutputModel):
    pass


class DFTPhWorkflow(Workflow):

    output_model = DFTPhOutput
    outputs: DFTPhOutput

    def _run(self) -> None:

        self.print('Calculate the dielectric tensor', style='heading')

        calc_scf = self.new_calculator('pw', nbnd=None)
        calc_scf.prefix = 'scf'
        status = self.run_steps(calc_scf)
        if status != Status.COMPLETED:
            return

        calc_ph = self.new_calculator('ph', epsil=True, fildyn=f'{self.name}.dynG')
        calc_ph.prefix = 'eps'
        self.link(calc_scf, calc_scf.parameters.outdir, calc_ph, calc_ph.parameters.outdir, symlink=True)
        status = self.run_steps(calc_ph)

        if status == Status.COMPLETED:
            self.status = Status.COMPLETED

        return


class DFTBandsOutput(OutputModel):
    pass


class DFTBandsWorkflow(DFTWorkflow):

    output_model = DFTBandsOutput
    outputs: DFTBandsOutput

    def _run(self) -> None:

        self.print('DFT bandstructure workflow', style='heading')

        # First, a scf calculation
        calc_scf: calculators.PWCalculator = self.new_calculator('pw', nbnd=None)
        calc_scf.prefix = 'scf'
        status = self.run_steps(calc_scf)
        if status != Status.COMPLETED:
            return

        # Second, a bands calculation
        if self.parameters.calculate_bands in (True, None):
            calc_bands: calculators.PWCalculator = self.new_calculator(
                'pw', calculation='bands', kpts=self.kpoints.path)
        else:
            calc_bands = self.new_calculator('pw', calculation='nscf')
        calc_bands.prefix = 'bands'
        assert isinstance(calc_scf.parameters.outdir, Path)
        assert isinstance(calc_bands.parameters.outdir, Path)
        self.link(calc_scf, calc_scf.parameters.outdir, calc_bands, calc_bands.parameters.outdir, symlink=True)
        status = self.run_steps(calc_bands)
        if status != Status.COMPLETED:
            return

        # Prepare the band structure for plotting
        if self.parameters.calculate_bands in (True, None):
            bs = calc_bands.results['band structure']

        # Third, a PDOS calculation
        if all([p['header'].get('number_of_wfc', 0) for p in self.pseudopotentials.values()]):
            calc_dos: calculators.ProjwfcCalculator = self.new_calculator('projwfc')
            assert isinstance(calc_dos.parameters.outdir, Path)
            self.link(calc_bands, calc_bands.parameters.outdir, calc_dos, calc_dos.parameters.outdir, symlink=True)
            status = self.run_steps(calc_dos)
            if status != Status.COMPLETED:
                return

            # Prepare the DOS for plotting
            dos = copy.deepcopy(calc_dos.results['dos'])
            if self.parameters.calculate_bands in (True, None):
                dos._energies -= bs.reference
        else:
            # Skip if the pseudos don't have the requisite PP_PSWFC blocks
            utils.warn('Some of the pseudopotentials do not have `PP_PSWFC` blocks, which means a projected DOS '
                       'calculation is not possible. Skipping...')
            dos = None

        # Plot the band structure and DOS
        if self.parameters.calculate_bands in (True, None):
            self.plot_bandstructure(bs.subtract_reference(), dos)
        elif dos is not None:
            workflow_name = self.__class__.__name__.lower()
            filename = f'{self.name}_{workflow_name}_dos'
            dos.plot(filename=filename)

        self.status = Status.COMPLETED

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
