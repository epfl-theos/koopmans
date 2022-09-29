"""

Workflow module for koopmans, containing the workflow for performing KI and KIPZ calculations with KC_WANN

Written by Edward Linscott Feb 2021

"""

import os
import shutil
from pathlib import Path

import numpy as np

from koopmans import utils
from koopmans.bands import Bands
from koopmans.calculators import (KoopmansHamCalculator, PWCalculator,
                                  Wann2KCCalculator)

from ._workflow import Workflow


class KoopmansDFPTWorkflow(Workflow):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Check the consistency of keywords
        if self.parameters.spin_polarized:
            raise NotImplementedError(
                'Calculating screening parameters with DFPT is not yet possible for spin-polarized systems')
        if self.parameters.functional != 'ki':
            raise NotImplementedError(
                'Calculating screening parameters with DFPT is not yet possible with functionals other than KI')
        if all(self.atoms.pbc):
            if self.parameters.init_orbitals not in ['mlwfs', 'projwfs']:
                raise ValueError(
                    'Calculating screening parameters with DFPT for a periodic system is only possible with MLWFs '
                    'or projected WFs as the variational orbitals')
        else:
            if self.parameters.init_orbitals != 'kohn-sham':
                raise ValueError(
                    'Calculating screening parameters with DFPT for a non-periodic system is only possible '
                    'with Kohn-Sham orbitals as the variational orbitals')
        for params in self.calculator_parameters.values():
            if all(self.atoms.pbc):
                # Gygi-Baldereschi
                if self.parameters.gb_correction:
                    if hasattr(params, 'l_vcut'):
                        if params.l_vcut is None:
                            params.l_vcut = True
                        if not params.l_vcut:
                            raise ValueError('Gygi-Baldereschi corrections require l_vcut = True')

                # Makov-Payne
                if self.parameters.mp_correction:
                    if hasattr(params, 'eps_inf'):
                        params.eps_inf = self.parameters.eps_inf

                # Martyna-Tuckerman (not used for periodic systems)
                if hasattr(params, 'assume_isolated'):
                    if params.assume_isolated not in ['none', None]:
                        raise ValueError('Periodic systems must have "assume_isolated" = "none", but it is set to '
                                         + params.assume_isolated)

            else:
                # Gygi-Baldereschi (not used for aperiodic systems)
                if hasattr(params, 'l_vcut'):
                    if params.l_vcut is None:
                        params.l_vcut = False
                    if params.l_vcut:
                        raise ValueError('Aperiodic systems require l_vcut = False')

                # Makov-Payne (not used for aperiodic systems)
                if getattr(params, 'eps_inf', None) not in [None, 1.0]:
                    raise ValueError('Aperiodic systems should not have a non-zero eps_inf')

                # Martyna-Tuckerman
                if self.parameters.mt_correction:
                    if hasattr(params, 'assume_isolated'):
                        if params.assume_isolated is None:
                            params.assume_isolated = 'm-t'
                        if params.assume_isolated not in ['martyna-tuckerman', 'm-t', 'mt']:
                            raise ValueError(
                                f'assume_isolated = {params.assume_isolated} is incompatible with mt_correction = True')

        # Initialize the bands
        if all(self.atoms.pbc):
            nocc = self.projections.num_wann(occ=True)
            nemp = self.projections.num_wann(occ=False)
        else:
            nocc = self.calculator_parameters['kcp'].nelec // 2
            nemp = self.calculator_parameters['pw'].nbnd - nocc
        if self.parameters.orbital_groups is None:
            self.parameters.orbital_groups = [list(range(nocc + nemp))]
        self.bands = Bands(n_bands=nocc + nemp, filling=[[True] * nocc + [False] * nemp],
                           groups=self.parameters.orbital_groups)

        # Populating kpoints if absent
        if not all(self.atoms.pbc):
            for key in ['pw', 'kc_screen']:
                self.calculator_parameters[key].kpts = [1, 1, 1]

    def _run(self):
        '''
        This function runs the workflow from start to finish
        '''

        # Import these here so that if they have been monkey-patched, we get the monkey-patched version
        from koopmans.workflows import DFTPWWorkflow, WannierizeWorkflow

        if self.parameters.from_scratch:
            for output_directory in [self.calculator_parameters['pw'].outdir, 'wannier', 'init', 'screening',
                                     'hamiltonian', 'postproc']:
                output_directory = Path(output_directory)
                if output_directory.exists():
                    shutil.rmtree(output_directory)

        if all(self.atoms.pbc):
            # Run PW and Wannierization
            for key in self.calculator_parameters.keys():
                if key.startswith('w90'):
                    self.calculator_parameters[key].write_u_matrices = True
                    self.calculator_parameters[key].write_xyz = True
            wf_workflow = WannierizeWorkflow.fromparent(self, force_nspin2=True)
            wf_workflow.run()

        else:
            # Run PW
            self.print('Initialization of density and variational orbitals', style='heading')

            # Create the workflow
            pw_workflow = DFTPWWorkflow.fromparent(self)

            # Update settings
            pw_params = pw_workflow.calculator_parameters['pw']
            pw_params.nspin = 2
            pw_params.tot_magnetization = 0

            # Run the subworkflow
            with utils.chdir('init'):
                pw_workflow.run()

        # Copy the outdir to the base directory
        base_outdir = self.calculator_parameters['pw'].outdir
        base_outdir.mkdir(exist_ok=True)
        scf_calcs = [c for c in self.calculations if isinstance(c, PWCalculator) and c.parameters.calculation == 'scf']
        init_outdir = scf_calcs[-1].parameters.outdir
        if self.parameters.from_scratch and init_outdir != base_outdir:
            utils.symlink(f'{init_outdir}/*', base_outdir)

        # Convert from wannier to KC
        self.print('Conversion to Koopmans format', style='subheading')
        wann2kc_calc = self.new_calculator('wann2kc')
        self.run_calculator(wann2kc_calc)

        # Calculate screening parameters
        if self.parameters.calculate_alpha:
            self.print('Calculation of screening parameters', style='heading')
            all_groups = [g for gspin in self.parameters.orbital_groups for g in gspin]
            if len(set(all_groups)) == len(all_groups):
                # If there is no orbital grouping, do all orbitals in one calculation
                # 1) Create the calculator
                kc_screen_calc = self.new_calculator('kc_screen')

                # 2) Run the calculator
                self.run_calculator(kc_screen_calc)

                # 3) Store the computed screening parameters
                self.bands.alphas = kc_screen_calc.results['alphas']
            else:
                # If there is orbital grouping, do the orbitals one-by-one
                for band in self.bands.to_solve:
                    # 1) Create the calculator (in a subdirectory)
                    kc_screen_calc = self.new_calculator('kc_screen', i_orb=band.index)
                    kc_screen_calc.directory /= f'band_{band.index}'

                    # 2) Run the calculator
                    self.run_calculator(kc_screen_calc)

                    # 3) Store the computed screening parameter (accounting for band groupings)
                    for b in self.bands:
                        if b.group == band.group:
                            alpha = kc_screen_calc.results['alphas'][band.spin]
                            b.alpha = alpha[band.spin]
        else:
            # Load the alphas
            if self.parameters.alpha_from_file:
                self.bands.alphas = [utils.read_alpha_file(Path())]
            else:
                self.bands.alphas = self.parameters.alpha_guess

        # Calculate the Hamiltonian
        self.print('Construction of the Hamiltonian', style='heading')
        kc_ham_calc = self.new_calculator('kc_ham', kpts=self.kpoints.path)

        if self.parameters.calculate_alpha and kc_ham_calc.parameters.lrpa != kc_screen_calc.parameters.lrpa:
            raise ValueError('Do not set "lrpa" to different values in the "screen" and "ham" blocks')
        self.run_calculator(kc_ham_calc)

        # Postprocessing
        if all(self.atoms.pbc) and self.projections and self.kpoints.path is not None \
                and self.calculator_parameters['ui'].do_smooth_interpolation:
            from koopmans.workflows import UnfoldAndInterpolateWorkflow
            self.print(f'\nPostprocessing', style='heading')
            ui_workflow = UnfoldAndInterpolateWorkflow.fromparent(self)
            ui_workflow.run(subdirectory='postproc')

        # Plotting
        self.plot_bandstructure()

    def plot_bandstructure(self):
        if not all(self.atoms.pbc):
            return

        # Identify the relevant calculators
        [wann2kc_calc] = [c for c in self.calculations if isinstance(c, Wann2KCCalculator)]
        [kc_ham_calc] = [c for c in self.calculations if isinstance(c, KoopmansHamCalculator)]

        # Plot the bandstructure if the band path has been specified
        bs = kc_ham_calc.results['band structure']
        if bs.path.path:
            super().plot_bandstructure(bs.subtract_reference())

    def new_calculator(self, calc_presets, **kwargs):
        if calc_presets not in ['kc_ham', 'kc_screen', 'wann2kc']:
            raise ValueError(
                f'Invalid choice calc_presets={calc_presets} in {self.__class__.__name__}.new_calculator()')

        calc = super().new_calculator(calc_presets)

        calc.prefix = 'kc'
        calc.parameters.prefix = 'kc'
        calc.parameters.outdir = 'TMP'
        calc.parameters.seedname = 'wann'
        calc.parameters.spin_component = 1
        calc.parameters.kcw_at_ks = not all(self.atoms.pbc)
        calc.parameters.read_unitary_matrix = all(self.atoms.pbc)

        if calc_presets == 'wann2kc':
            if all(self.atoms.pbc):
                calc.directory = 'wannier'
            else:
                calc.directory = 'init'
        elif calc_presets == 'kc_screen':
            calc.directory = 'screening'
            # If eps_inf is not provided in the kc_wann:screen subdictionary but there is a value provided in the
            # workflow parameters, adopt that value
            if self.parameters.eps_inf is not None and calc.parameters.eps_inf is None and all(self.atoms.pbc):
                calc.parameters.eps_inf = self.parameters.eps_inf
        else:
            calc.directory = 'hamiltonian'
            calc.parameters.do_bands = all(self.atoms.pbc)
            calc.alphas = self.bands.alphas

        if all(self.atoms.pbc):
            nocc = self.projections.num_wann(occ=True)
            nemp = self.projections.num_wann(occ=False)
            have_empty = (self.projections.num_wann(occ=False) > 0)
            has_disentangle = (self.projections.num_bands() != self.projections.num_wann())
        else:
            nocc = self.calculator_parameters['kcp'].nelec // 2
            nemp = self.calculator_parameters['pw'].nbnd - nocc
            have_empty = nemp > 0
            has_disentangle = False
        calc.parameters.num_wann_occ = nocc
        calc.parameters.num_wann_emp = nemp
        calc.parameters.have_empty = have_empty
        calc.parameters.has_disentangle = has_disentangle

        # Apply any additional calculator keywords passed as kwargs
        for k, v in kwargs.items():
            setattr(calc.parameters, k, v)

        return calc

    def run_calculator(self, calc):
        # Create this (possibly nested) directory
        calc.directory.mkdir(parents=True, exist_ok=True)

        # Provide the rotation matrices and the wannier centers
        if all(self.atoms.pbc):
            utils.symlink(f'wannier/occ/wann_u.mat', f'{calc.directory}/', exist_ok=True)
            utils.symlink(f'wannier/emp/wann_u.mat', f'{calc.directory}/wann_emp_u.mat', exist_ok=True)
            if Path('wannier/emp/wann_u_dis.mat').exists():
                utils.symlink(f'wannier/emp/wann_u_dis.mat',
                              f'{calc.directory}/wann_emp_u_dis.mat', exist_ok=True)
            utils.symlink(f'wannier/occ/wann_centres.xyz', f'{calc.directory}/', exist_ok=True)
            utils.symlink(f'wannier/emp/wann_centres.xyz',
                          f'{calc.directory}/wann_emp_centres.xyz', exist_ok=True)

        super().run_calculator(calc)
