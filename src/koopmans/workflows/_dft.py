"""Workflows for performing simple DFT calculations with either kcp.x or pw.x."""

import copy
from pathlib import Path
from typing import TypeVar

from ase_koopmans.spectrum.band_structure import BandStructure
from ase_koopmans.spectrum.doscollection import DOSCollection

from koopmans import calculators
from koopmans.files import File
from koopmans.process_io import IOModel
from koopmans.status import Status
from koopmans.utils.warnings import warn

from ._workflow import Workflow, spin_symmetrize

T = TypeVar('T', bound='calculators.CalculatorExt')
OutputModel = TypeVar('OutputModel', bound=IOModel)


class DFTWorkflow(Workflow[OutputModel]):
    """Basic template for a DFT workflow."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.parameters.functional = 'dft'


class DFTCPOutput(IOModel):
    """Output model for the DFTCPWorkflow."""

    pass


class DFTCPWorkflow(DFTWorkflow[DFTCPOutput]):
    """Basic workflow for performing an SCF DFT calculation with kcp.x."""

    output_model = DFTCPOutput

    def _run(self) -> None:
        """Run the workflow."""
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

        if self.parameters.fix_spin_contamination:
            status = spin_symmetrize(self, calc)
            if status == Status.COMPLETED:
                self.status = status
        else:
            status = self.run_steps(calc)
            if status == Status.COMPLETED:
                self.status = status

        return


class DFTPWOutput(IOModel):
    """Output model for the DFTPWWorkflow."""

    pass


class DFTPWWorkflow(DFTWorkflow[DFTPWOutput]):
    """Basic workflow for performing an SCF DFT calculation with pw.x."""

    output_model = DFTPWOutput

    def _run(self) -> None:
        """Run the workflow."""
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


class DFTPhOutput(IOModel):
    """Output model for the DFTPhWorkflow."""

    pass


class DFTPhWorkflow(Workflow[DFTPhOutput]):
    """Workflow for calculating the dielectric tensor using ph.x."""

    output_model = DFTPhOutput

    def _run(self) -> None:
        """Run the workflow."""
        calc_scf = self.new_calculator('pw', nbnd=None)
        calc_scf.prefix = 'scf'
        status = self.run_steps(calc_scf)
        if status != Status.COMPLETED:
            return

        calc_ph = self.new_calculator('ph', epsil=True, fildyn=f'{self.name}.dynG')
        calc_ph.prefix = 'eps'
        calc_ph.link(File(calc_scf, calc_scf.parameters.outdir), calc_ph.parameters.outdir, symlink=True)
        status = self.run_steps(calc_ph)

        if status == Status.COMPLETED:
            self.status = Status.COMPLETED

        return


class DFTBandsOutput(IOModel):
    """Output model for the DFTBandsWorkflow."""

    band_structure: BandStructure
    dos: DOSCollection | None


class DFTBandsWorkflow(DFTWorkflow[DFTBandsOutput]):
    """Workflow for calculating the band structure and projected density of states using pw.x."""

    output_model = DFTBandsOutput

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
        calc_bands.link(File(calc_scf, calc_scf.parameters.outdir), calc_bands.parameters.outdir, symlink=True)
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
            calc_dos.link(File(calc_bands, calc_bands.parameters.outdir), calc_dos.parameters.outdir, symlink=True)
            status = self.run_steps(calc_dos)
            if status != Status.COMPLETED:
                return

            # Prepare the DOS for plotting
            dos = copy.deepcopy(calc_dos.results['dos'])
            if self.parameters.calculate_bands in (True, None):
                dos._energies -= bs.reference
        else:
            # Skip if the pseudos don't have the requisite PP_PSWFC blocks
            warn('Some of the pseudopotentials do not have `PP_PSWFC` blocks, which means a projected DOS '
                 'calculation is not possible. Skipping...')
            dos = None

        # Plot the band structure and DOS
        if self.parameters.calculate_bands in (True, None):
            self.plot_bandstructure(bs.subtract_reference(), dos)
        elif dos is not None:
            workflow_name = self.__class__.__name__.lower()
            filename = f'{self.name}_{workflow_name}_dos'
            dos.plot(filename=filename)

        self.outputs = self.output_model(band_structure=bs, dos=dos)

        self.status = Status.COMPLETED

    def new_calculator(self,
                       calc_type: str,
                       *args,
                       **kwargs) -> T:   # type: ignore[type-var, misc]
        """Create a new calculator."""
        calc: T = super().new_calculator(calc_type, *args, **kwargs)
        if calc_type == 'projwfc':
            assert isinstance(calc, calculators.ProjwfcCalculator)
            calc.parameters.filpdos = self.name
            calc.pseudopotentials = self.pseudopotentials
            calc.spin_polarized = self.parameters.spin_polarized
        return calc
