"""

Workflow module for koopmans, containing the workflow for performing KI and KIPZ calculations with KC_WANN

Written by Edward Linscott Feb 2021

"""

import os
import numpy as np
import copy
import matplotlib
from koopmans import utils, io
from koopmans.bands import Bands
from koopmans.calculators import PWCalculator, Wann2KCCalculator, KoopmansScreenCalculator, KoopmansHamCalculator
from koopmans.workflows.generic import Workflow
from koopmans.workflows.wf_with_w90 import WannierizeWorkflow
from koopmans.workflows.dft_with_pw import DFTPWWorkflow
matplotlib.use('Agg')


class KoopmansPWWorkflow(Workflow):

    def __init__(self, workflow_settings, calcs_dct):
        super().__init__(workflow_settings, calcs_dct)

        # Check the consistency of keywords
        if self.functional != 'ki':
            raise ValueError('Calculating screening parameters with DFPT is only possible with the KI functional')
        if self.periodic:
            if self.init_orbitals not in ['mlwfs', 'projwfs']:
                raise ValueError(
                    'Calculating screening parameters with DFPT for a periodic system is only possible with MLWFs as '
                    'the variational orbitals')
        else:
            if self.init_orbitals != 'kohn-sham':
                raise ValueError(
                    'Calculating screening parameters with DFPT for a periodic system is only possible with Kohn-Sham '
                    'orbitals as the variational orbitals')
        for calc in self.master_calcs.values():
            if self.periodic:
                # Gygi-Baldereschi
                if self.gb_correction:
                    if hasattr(calc, 'l_vcut'):
                        if calc.l_vcut is None:
                            calc.l_vcut = True
                        if not calc.l_vcut:
                            raise ValueError('Gygi-Baldereschi corrections require l_vcut = True')

                # Makov-Payne
                if self.mp_correction:
                    if hasattr(calc, 'eps_inf'):
                        calc.eps_inf = self.eps_inf

                # Martyna-Tuckerman (not used for periodic systems)
                if hasattr(calc, 'assume_isolated'):
                    if calc.assume_isolated not in ['none', None]:
                        raise ValueError(
                            f'Periodic systems must have "assume_isolated" = "none", but it is set to {val}')

            else:
                # Gygi-Baldereschi (not used for aperiodic systems)
                if hasattr(calc, 'l_vcut'):
                    if calc.l_vcut is None:
                        calc.l_vcut = False
                    if calc.l_vcut:
                        raise ValueError('Aperiodic systems require l_vcut = False')

                # Makov-Payne (not used for aperiodic systems)
                if getattr(calc, 'eps_inf', None) not in [None, 1.0]:
                    raise ValueError('Aperiodic systems should not have a non-zero eps_inf')

                # Martyna-Tuckerman
                if self.mt_correction:
                    if hasattr(calc, 'assume_isolated'):
                        if calc.assume_isolated is None:
                            calc.assume_isolated = 'm-t'
                        if calc.assume_isolated not in ['martyna-tuckerman', 'm-t', 'mt']:
                            raise ValueError(
                                f'assume_isolated = {calc.assume_isolated} is incompatible with mt_correction = True')

        # Delete any pre-existing directories if running from scratch
        if self.from_scratch:
            for directory in ['init', 'wannier', 'screening', 'hamiltonian', 'TMP']:
                if os.path.isdir(directory):
                    utils.system_call(f'rm -r {directory}')

        # Initialise the bands
        if self.periodic:
            nocc = self.master_calcs['w90_occ'].num_wann
            nemp = self.master_calcs['w90_emp'].num_wann
        else:
            pw_calc = self.master_calcs['pw']
            nocc = pw_calc.nelec // 2
            nemp = pw_calc.nbnd - nocc
        if self.orbital_groups is None:
            self.orbital_groups = list(range(nocc + nemp))
        self.bands = Bands(n_bands=nocc + nemp, filling=[True] * nocc + [False] * nemp, groups=self.orbital_groups)

        # Populating kpoints if absent
        if not self.periodic:
            for calc in self.master_calcs.values():
                if isinstance(calc, (PWCalculator, KoopmansScreenCalculator)):
                    calc.calc.parameters['kpts'] = [1, 1, 1]

    def run(self):
        '''
        This function runs the workflow from start to finish
        '''

        if self.periodic:
            # Run PW and Wannierisation
            for key in self.master_calcs.keys():
                if key.startswith('w90'):
                    self.master_calcs[key].write_u_matrices = True
                    self.master_calcs[key].write_xyz = True
            wf_workflow = WannierizeWorkflow(self.settings, self.master_calcs, nspin=2)
            self.run_subworkflow(wf_workflow)

        else:
            # Run PW
            self.print('Initialisation of density and variational orbitals', style='heading')

            # Create the workflow
            pw_workflow = DFTPWWorkflow(self.settings, self.master_calcs)

            # Update settings
            pw_calc = pw_workflow.master_calcs['pw']
            pw_calc.directory = 'init'
            pw_calc.outdir = '../TMP'
            pw_calc.nspin = 2
            pw_calc.tot_magnetization = 0

            # Run the subworkflow
            self.run_subworkflow(pw_workflow)

        # Copy the outdir to the base directory
        base_outdir = self.master_calcs['pw'].outdir
        if not os.path.isdir(base_outdir):
            utils.system_call(f'mkdir {base_outdir}')
        init_outdir = self.all_calcs[0].outdir
        if self.from_scratch and init_outdir != base_outdir:
            utils.system_call(f'ln -sf {init_outdir}/* {base_outdir}/')

        # Convert from wannier to KC
        self.print('Conversion to Koopmans format', style='subheading')
        wann2kc_calc = self.new_calculator('wann2kc')
        self.run_calculator(wann2kc_calc)

        # Calculate screening parameters
        if self.calculate_alpha:
            self.print('Calculation of screening parameters', style='heading')
            if len(set(self.orbital_groups)) == len(self.orbital_groups):
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
                    kc_screen_calc.directory += f'/band_{band.index}'

                    # 2) Run the calculator
                    self.run_calculator(kc_screen_calc)

                    # 3) Store the computed screening parameter (accounting for band groupings)
                    for b in self.bands:
                        if b.group == band.group:
                            b.alpha = kc_screen_calc.results['alphas'][0]
        else:
            # Load the alphas
            if self.alpha_from_file:
                self.bands.alphas = io.read_alpha_file()
            else:
                self.bands.alphas = self.alpha_guess

        # Calculate the Hamiltonian
        self.print('Construction of the Hamiltonian', style='heading')
        kc_ham_calc = self.new_calculator('kc_ham')

        if self.calculate_alpha and kc_ham_calc.lrpa != kc_screen_calc.lrpa:
            raise ValueError('Do not set "lrpa" to different values in the "screen" and "ham" blocks')
        self.run_calculator(kc_ham_calc)

        # Plotting
        self.plot_bandstructure()

    def plot_bandstructure(self):
        # Identify the relevant calculators
        [wann2kc_calc] = [c for c in self.all_calcs if isinstance(c, Wann2KCCalculator)]
        [kc_ham_calc] = [c for c in self.all_calcs if isinstance(c, KoopmansHamCalculator)]

        if self.periodic:
            # Align the bandstructure to the VBM
            bs = kc_ham_calc.results['band structure']
            reference = None
            n_occ = wann2kc_calc.num_wann_occ
            if n_occ is not None:
                reference = np.max(bs.energies, axis=1)[:, n_occ - 1][0]
                bs._reference = reference
                bs = bs.subtract_reference()

            # Plot the bandstructure if the band path has been specified
            if bs.path.path:
                bs.plot(emin=-20, emax=20, filename=f'{self.name}_bandstructure.png')

    def new_calculator(self, calc_presets, **kwargs):
        if calc_presets not in ['kc_ham', 'kc_screen', 'wann2kc']:
            raise ValueError(
                f'Invalid choice calc_presets={calc_presets} in {self.__class__.__name__}.new_calculator()')

        calc = copy.deepcopy(self.master_calcs[calc_presets])
        calc.name = 'kc'
        calc.outdir = 'TMP'
        calc.seedname = 'wann'
        calc.spin_component = 1
        calc.kc_at_ks = not self.periodic
        calc.read_unitary_matrix = self.periodic

        if calc_presets == 'wann2kc':
            if self.periodic:
                calc.directory = 'wannier'
            else:
                calc.directory = 'init'
            calc.calc.parameters.pop('kpath', None)
        elif calc_presets == 'kc_screen':
            calc.directory = 'screening'
            # If eps_inf is not provided in the kc_wann:screen subdictionary but there is a value provided in the
            # workflow parameters, adopt that value
            if self.eps_inf is not None and calc.eps_inf is None and self.periodic:
                calc.eps_inf = self.eps_inf
        else:
            calc.directory = 'hamiltonian'
            calc.results['alphas'] = self.bands.alphas

        calc.num_wann_occ = self.master_calcs['w90_occ'].num_wann
        calc.num_wann_emp = self.master_calcs['w90_emp'].num_wann

        # Apply any additional calculator keywords passed as kwargs
        for k, v in kwargs.items():
            setattr(calc, k, v)

        return calc

    def run_calculator(self, calc):
        # Create this (possibly nested) directory
        utils.mkdir(calc.directory)

        # Provide the rotation matrices and the wannier centres
        if self.periodic:
            utils.system_call(f'ln -srf wannier/occ/wann_u.mat {calc.directory}/')
            utils.system_call(f'ln -srf wannier/emp/wann_u.mat {calc.directory}/wann_emp_u.mat')
            utils.system_call(f'ln -srf wannier/emp/wann_u_dis.mat {calc.directory}/wann_emp_u_dis.mat')
            utils.system_call(f'ln -srf wannier/occ/wann_centres.xyz {calc.directory}/wann_centres.xyz')
            utils.system_call(f'ln -srf wannier/emp/wann_centres.xyz {calc.directory}/wann_emp_centres.xyz')

        super().run_calculator(calc)
