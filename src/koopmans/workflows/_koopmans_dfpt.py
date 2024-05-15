"""

Workflow module for koopmans, containing the workflow for performing KI and KIPZ calculations with KC_WANN

Written by Edward Linscott Feb 2021

"""

import shutil
from pathlib import Path
from typing import Dict

from koopmans import utils, pseudopotentials
from koopmans.bands import Bands
from koopmans.calculators import (KoopmansHamCalculator, PWCalculator,
                                  Wann2KCCalculator)

from ._workflow import Workflow


class KoopmansDFPTWorkflow(Workflow):

    def __init__(self, scf_kgrid=None, *args, **kwargs):
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
        if self.parameters.dfpt_coarse_grid is not None and not self.parameters.calculate_alpha:
            raise ValueError('Using a DFPT coarse grid to calculate the screening parameters is not possible '
                             'when calculate_alpha = False')
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
        nocc = pseudopotentials.nelec_from_pseudos(
            self.atoms, self.pseudopotentials, self.parameters.pseudo_directory) // 2
        if all(self.atoms.pbc):
            exclude_bands = self.calculator_parameters['w90'].get('exclude_bands', [])
            nocc -= len(exclude_bands)
            ntot = self.projections.num_wann()
        else:
            ntot = self.calculator_parameters['pw'].nbnd
        nemp = ntot - nocc
        if self.parameters.orbital_groups is None:
            self.parameters.orbital_groups = [list(range(nocc + nemp))]
        tols: Dict[str, float] = {}
        for key in ['self_hartree', 'spread']:
            val = self.parameters.get(f'orbital_groups_{key}_tol', None)
            if val is not None:
                tols[key] = val
        self.bands = Bands(n_bands=nocc + nemp, filling=[[True] * nocc + [False] * nemp],
                           groups=self.parameters.orbital_groups, tolerances=tols)

        # Populating kpoints if absent
        if not all(self.atoms.pbc):
            for key in ['pw', 'kc_screen']:
                self.calculator_parameters[key].kpts = [1, 1, 1]

        self._perform_ham_calc: bool = True

        self._scf_kgrid = scf_kgrid

    def _run(self):
        '''
        This function runs the workflow from start to finish
        '''

        # Import these here so that if they have been monkey-patched, we get the monkey-patched version
        from koopmans.workflows import DFTPWWorkflow, WannierizeWorkflow

        # MB mod
        if self.parameters.from_scratch and self.parameters.mode == "ase":
            for output_directory in [self.calculator_parameters['pw'].outdir, 'wannier', 'init', 'screening',
                                     'hamiltonian', 'postproc']:
                output_directory = Path(output_directory)
                if output_directory.exists():
                    shutil.rmtree(output_directory)

        if self.parameters.dfpt_coarse_grid is not None:
            self.print('Coarse grid calculations', style='heading')
            coarse_wf = self.__class__.fromparent(self, scf_kgrid=self.kpoints.grid)
            coarse_wf.parameters.dfpt_coarse_grid = None
            coarse_wf.kpoints.grid = self.parameters.dfpt_coarse_grid
            coarse_wf._perform_ham_calc = False
            coarse_wf.run(subdirectory='coarse_grid')
            self.print('Regular grid calculations', style='heading')

        if all(self.atoms.pbc):
            # Run PW and Wannierization
            for key in self.calculator_parameters.keys():
                if key.startswith('w90'):
                    self.calculator_parameters[key].write_u_matrices = True
                    self.calculator_parameters[key].write_xyz = True
            wf_workflow = WannierizeWorkflow.fromparent(self, force_nspin2=True, scf_kgrid = self._scf_kgrid)
            wf_workflow.run()
            
            # MB mod
            if hasattr(wf_workflow,"dft_wchains"): self.dft_wchains = wf_workflow.dft_wchains
            if hasattr(wf_workflow,"w90_wchains"): self.w90_wchains = wf_workflow.w90_wchains
            if hasattr(wf_workflow,"wannier90_files"): self.wannier90_files = wf_workflow.wannier90_files

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
            
            # MB mod
            if pw_workflow.parameters.mode == "ase":
                with utils.chdir('init'):
                    pw_workflow.run()
            else:
                pw_workflow.run()
                # MB mod
                if hasattr(pw_workflow,"dft_wchains"): self.dft_wchains = pw_workflow.dft_wchains
                
        # MB mod
        if self.parameters.mode == "ase":
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
        
        # MB mod
        if not self.parameters.mode == "ase":
            if hasattr(self,"wannier90_files"): wann2kc_calc.wannier90_files = self.wannier90_files           
        
        self.run_calculator(wann2kc_calc)
        
        # MB mod
        if not self.parameters.mode == "ase":
            self.wann2kc_calculation = wann2kc_calc.calculation
        
        # Calculate screening parameters
        if self.parameters.calculate_alpha:
            if self.parameters.dfpt_coarse_grid is None:
                self.print('Calculation of screening parameters', style='heading')

                # Group the bands by spread
                self.bands.assign_groups(sort_by='spread', allow_reassignment=True)

                if len(self.bands.to_solve) == len(self.bands):
                    # If there is no orbital grouping, do all orbitals in one calculation
                    # 1) Create the calculator
                    kc_screen_calc = self.new_calculator('kc_screen')
                    
                    # MB mod
                    if not self.parameters.mode == "ase":
                        kc_screen_calc.parent_folder = self.wann2kc_calculation.outputs.remote_folder
                        if hasattr(self,"wannier90_files"): kc_screen_calc.wannier90_files = self.wannier90_files  

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

                        # MB mod
                        if not self.parameters.mode == "ase":
                            kc_screen_calc.parent_folder = self.wann2kc_calculation.outputs.remote_folder
                            if hasattr(self,"wannier90_files"): kc_screen_calc.wannier90_files = self.wannier90_files
                        
                        # 2) Run the calculator
                        self.run_calculator(kc_screen_calc)

                        # 3) Store the computed screening parameter (accounting for band groupings)
                        for b in self.bands:
                            if b.group == band.group:
                                alpha = kc_screen_calc.results['alphas'][band.spin]
                                b.alpha = alpha[band.spin]
            else:
                self.bands.alphas = coarse_wf.bands.alphas
            
            # MB mod
            if not self.parameters.mode == "ase":
                self.kc_screen_calculation = kc_screen_calc.calculation
        
        else:
            # Load the alphas
            if self.parameters.alpha_from_file:
                self.bands.alphas = [utils.read_alpha_file(Path())]
            else:
                self.bands.alphas = self.parameters.alpha_guess

        # Calculate the Hamiltonian
        if self._perform_ham_calc:
            self.print('Construction of the Hamiltonian', style='heading')
            kc_ham_calc = self.new_calculator('kc_ham', kpts=self.kpoints.path)

            if self.parameters.calculate_alpha and self.parameters.dfpt_coarse_grid is None:
                if kc_ham_calc.parameters.lrpa != kc_screen_calc.parameters.lrpa:
                    raise ValueError('Do not set "lrpa" to different values in the "screen" and "ham" blocks')
            
            # MB mod
            if not self.parameters.mode == "ase":
                if hasattr(self, "kc_screen_calculation"):
                    kc_ham_calc.parent_folder = self.kc_screen_calculation.outputs.remote_folder
                else:
                    kc_ham_calc.parent_folder = self.wann2kc_calculation.outputs.remote_folder
                if hasattr(self,"wannier90_files"): 
                    kc_ham_calc.wannier90_files = self.wannier90_files
                    
                    #the kcw.x internal interpolation is done on the same path of the DFT wannierized
                    kc_ham_calc.kpoints = self.w90_wchains["occ"][0].outputs.band_structure.get_kpoints() # just an array for now.
            
            self.run_calculator(kc_ham_calc)
            
            # MB mod
            if not self.parameters.mode == "ase":
                self.kc_ham_calculation = kc_ham_calc.calculation

            # MB mod
            if self.parameters.mode == "ase":
                # Postprocessing
                if all(self.atoms.pbc) and self.projections and self.kpoints.path is not None \
                        and self.calculator_parameters['ui'].do_smooth_interpolation:
                    from koopmans.workflows import UnfoldAndInterpolateWorkflow
                    self.print(f'\nPostprocessing', style='heading')
                    ui_workflow = UnfoldAndInterpolateWorkflow.fromparent(self)
                    ui_workflow.run(subdirectory='postproc')

                # Plotting
                self.plot_bandstructure()
            else:
                """
                # Postprocessing
                if all(self.atoms.pbc) and self.projections and self.kpoints.path is not None \
                        and self.calculator_parameters['ui'].do_smooth_interpolation:
                    from koopmans.workflows import UnfoldAndInterpolateWorkflow
                    self.print(f'\nPostprocessing', style='heading')
                    ui_workflow = UnfoldAndInterpolateWorkflow.fromparent(self)
                    ui_workflow.run(subdirectory='postproc')
                    if hasattr(ui_workflow,"w90_wchains"): self.w90_wchains_kcw = ui_workflow.w90_wchains
                """
                pass    

                # Plotting
                #self.plot_bandstructure()
                
                if hasattr(self,'dft_wchains'): 
                    self.dft_wchains_pk = [] 
                    for label,wchain in self.dft_wchains.items():
                        self.dft_wchains_pk.append(wchain.pk)
                    del self.dft_wchains
                if hasattr(self, 'w90_wchains'):
                    self.w90_wchains_pk = [] 
                    for label,wchains in self.w90_wchains.items():
                        for wchain in wchains:
                            self.w90_wchains_pk.append(wchain.pk)
                    del self.w90_wchains
                if hasattr(self, 'w90_wchains_kcw'):
                    self.w90_wchains_kcw_pk = [] 
                    for label,wchains in self.w90_wchains_kcw.items():
                        for wchain in wchains:
                            self.w90_wchains_kcw_pk.append(wchain.pk)
                    del self.w90_wchains_kcw
                if hasattr(self,'wann2kc_calculation'): del self.wann2kc_calculation
                if hasattr(self,'kc_screen_calculation'): del self.kc_screen_calculation
                if hasattr(self,'kc_ham_calculation'): 
                    self.kc_ham_calculation = self.kc_ham_calculation.pk
                    print(f"The last kcw AiiDA calculation was the one with pk <{self.kc_ham_calculation}>")

    def plot_bandstructure(self):
        if not all(self.atoms.pbc):
            return

        # Identify the relevant calculators
        wann2kc_calc = [c for c in self.calculations if isinstance(c, Wann2KCCalculator)][-1]
        kc_ham_calc = [c for c in self.calculations if isinstance(c, KoopmansHamCalculator)][-1]

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
            nocc = self.bands.num(filled=True)
            nemp = self.bands.num(filled=False)
            have_empty = (nemp > 0)
            has_disentangle = (self.projections.num_bands() != nocc + nemp)
        else:
            nocc = self.calculator_parameters['kcp'].nelec // 2
            nemp = self.calculator_parameters['pw'].nbnd - nocc
            have_empty = (nemp > 0)
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
        
        # MB mod
        if self.parameters.mode == "ase":
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
