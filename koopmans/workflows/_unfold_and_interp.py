"""

Workflow module for koopmans, containing the workflow for unfolding and interpolating
results from the supercell to generate a bandstructure

For the moment the code works only with cubic, tetragonal and orthorhombic systems.

Originally written by Riccardo De Gennaro as the standalone 'unfolding and interpolate' code
Integrated within koopmans by Edward Linscott Jan 2021
"""

import numpy as np
from pathlib import Path
from typing import Optional
from ._workflow import Workflow
from koopmans import calculators
from ase.spectrum.band_structure import BandStructure


class UnfoldAndInterpolateWorkflow(Workflow):

    def __init__(self, *args, redo_smooth_dft: Optional[bool] = None, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._redo_smooth_dft = redo_smooth_dft

    def _run(self):
        '''

        Wrapper for the whole unfolding and interpolation workflow, consisting of:
        - a smooth WannierizeWorkflow (if required)
        - an UnfoldAndInterpolateCalculator for occ states
        - an UnfoldAndInterpolateCalculator for emp states
        - an UnfoldAndInterpolateCalculator to merge occ and emp results

        '''
        # Import these here so that if these have been monkey-patched, we get the monkey-patched version
        from koopmans.workflows import WannierizeWorkflow

        # Transform self.atoms back to the primitive cell
        if self.parameters.method == 'dscf':
            self.supercell_to_primitive()

        # Store the original w90 calculations
        w90_calcs = [c for c in self.calculations if isinstance(c, calculators.Wannier90Calculator)
                     and c.command.flags == ''][-len(self.projections):]

        if self.master_calc_params['ui'].do_smooth_interpolation:
            wf_kwargs = self.wf_kwargs
            wf_kwargs['kgrid'] = [x * y for x,
                                  y in zip(wf_kwargs['kgrid'], self.master_calc_params['ui'].smooth_int_factor)]
            wannier_workflow = WannierizeWorkflow(**wf_kwargs)
            # For the moment not enabling this change
            # wannier_workflow = WannierizeWorkflow(scf_kgrid=self.kgrid, **wf_kwargs)

            # Here, we allow for skipping of the smooth dft calcs (assuming they have been already run)
            # This is achieved via the optional argument of from_scratch in run_subworkflow(), which
            # overrides the value of wannier_workflow.from_scratch, as well as preventing the inheritance of
            # self.from_scratch to wannier_workflow.from_scratch and back again after the subworkflow finishes
            self.run_subworkflow(wannier_workflow, from_scratch=self._redo_smooth_dft)

        calc: calculators.UnfoldAndInterpolateCalculator
        if self.parameters.spin_polarised:
            spins = ['up', 'down']
        else:
            spins = [None]

        for spin in spins:
            for filling in ['occ', 'emp']:
                calc_presets = filling
                if spin:
                    calc_presets += '_' + spin
                calc = self.new_ui_calculator(calc_presets)
                calc.centers = np.array([center for c in w90_calcs for center in c.results['centers']
                                        if calc_presets in c.directory.name])
                calc.spreads = [spread for c in w90_calcs for spread in c.results['spreads']
                                if calc_presets in c.directory.name]
                self.run_calculator(calc, enforce_ss=False)

        # Merge the two calculations to print out the DOS and bands
        calc = self.new_ui_calculator('merge')

        # Merge the bands
        if self.parameters.spin_polarised:
            energies = [[c.results['band structure'].energies for c in subset]
                        for subset in [self.calculations[-4:-2], self.calculations[-2:]]]
            reference = np.max([e[0] for e in energies])
            energies_np = np.concatenate([np.concatenate(e, axis=2) for e in energies], axis=0)
        else:
            energies = [c.results['band structure'].energies for c in self.calculations[-2:]]
            reference = np.max(energies[0])
            energies_np = np.concatenate(energies, axis=2)
        calc.results['band structure'] = BandStructure(self.kpath, energies_np - reference)

        if calc.parameters.do_dos:
            # Generate the DOS
            calc.calc_dos()

        # Print out the merged bands and DOS
        if self.parameters.from_scratch:
            calc.write_results()

        # Store the calculator in the workflow's list of all the calculators
        self.calculations.append(calc)

    def new_ui_calculator(self, calc_presets: str, **kwargs) -> calculators.UnfoldAndInterpolateCalculator:
        valid_calc_presets = ['occ', 'occ_up', 'occ_down', 'emp', 'emp_up', 'emp_down', 'merge']
        assert calc_presets in valid_calc_presets, \
            'In UnfoldAndInterpolateWorkflow.new_calculator() calc_presets must be ' \
            + '/'.join([f'"{s}"' for s in valid_calc_presets]) + \
            f', but you have tried to set it equal to {calc_presets}'

        if calc_presets == 'merge':
            # Dummy calculator for merging bands and dos
            kwargs['directory'] = Path('./')
            pass
        else:
            # Automatically generating UI calculator settings
            kwargs['directory'] = Path(f'{calc_presets}')
            kwargs['plot_params'] = self.plot_params
            if self.parameters.method == 'dscf':
                # DSCF case
                if '_' in calc_presets:
                    ham_prefix = calc_presets.replace('up', '1').replace('down', '2')
                else:
                    ham_prefix = calc_presets + '_1'
                kwargs['kc_ham_file'] = Path(f'../final/ham_{ham_prefix}.dat').resolve()
                kwargs['w90_seedname'] = Path(f'../init/wannier/{calc_presets}/wann').resolve()
                if self.master_calc_params['ui'].do_smooth_interpolation:
                    kwargs['dft_smooth_ham_file'] = Path(f'wannier/{calc_presets}/wann_hr.dat').resolve()
                    kwargs['dft_ham_file'] = Path(f'../init/wannier/{calc_presets}/wann_hr.dat').resolve()
            else:
                # DFPT case
                if self.parameters.spin_polarised:
                    raise NotImplementedError()
                kwargs['kc_ham_file'] = Path(f'../hamiltonian/kc.kcw_hr_{calc_presets}.dat').resolve()
                kwargs['w90_seedname'] = Path(f'../wannier/{calc_presets}/wann').resolve()
                if self.master_calc_params['ui'].do_smooth_interpolation:
                    kwargs['dft_smooth_ham_file'] = Path(f'wannier/{calc_presets}/wann_hr.dat').resolve()
                    kwargs['dft_ham_file'] = Path(f'../wannier/{calc_presets}/wann_hr.dat').resolve()

        calc: calculators.UnfoldAndInterpolateCalculator = super().new_calculator('ui', **kwargs)
        calc.prefix = self.parameters.functional

        return calc


class SingleUnfoldAndInterpolateWorkflow(Workflow):
    '''
    A workflow that runs a single UnfoldAndInterpolateCalculator

    This exists to make it possible for .uii (i.e. Unfold and Interpolate input) files to be run using the koopmans
    command
    '''

    def _run(self):
        '''
        '''
        ui_calc = self.new_calculator('ui')
        ui_calc.prefix = self.name

        ui_calc.calculate()
        self.calculations = [ui_calc]

        # Print quality control
        if self.parameters.print_qc and not ui_calc.skip_qc:
            for key in ui_calc.results_for_qc:
                val = ui_calc.results.get(key, None)
                if val is not None:
                    ui_calc.qc_results[key] = val
