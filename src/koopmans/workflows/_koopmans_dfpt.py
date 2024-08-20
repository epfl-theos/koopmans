"""

Workflow module for koopmans, containing the workflow for performing KI and KIPZ calculations with KC_WANN

Written by Edward Linscott Feb 2021

"""

import shutil
from pathlib import Path
from typing import Dict

from koopmans import pseudopotentials, utils
from koopmans.bands import Bands
from koopmans.calculators import (KoopmansHamCalculator, PWCalculator,
                                  Wann2KCCalculator, Wannier90Calculator)
from koopmans.files import FilePointer
from koopmans.outputs import OutputModel

from ._workflow import Workflow


class KoopmansDFPTOutputs(OutputModel):
    pass


class KoopmansDFPTWorkflow(Workflow):

    output_model = KoopmansDFPTOutputs  # type: ignore
    outputs: KoopmansDFPTOutputs

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
                             'when `calculate_alpha = False`')
        for params in self.calculator_parameters.values():
            if all(self.atoms.pbc):
                # Gygi-Baldereschi
                if self.parameters.gb_correction:
                    if hasattr(params, 'l_vcut'):
                        if params.l_vcut is None:
                            params.l_vcut = True
                        if not params.l_vcut:
                            raise ValueError('Gygi-Baldereschi corrections require `l_vcut = True`')

                # Makov-Payne
                if self.parameters.mp_correction:
                    if hasattr(params, 'eps_inf'):
                        params.eps_inf = self.parameters.eps_inf

                # Martyna-Tuckerman (not used for periodic systems)
                if hasattr(params, 'assume_isolated'):
                    if params.assume_isolated not in ['none', None]:
                        raise ValueError('Periodic systems must have `assume_isolated = "none"`, but it is set to `'
                                         + params.assume_isolated + '`')

            else:
                # Gygi-Baldereschi (not used for aperiodic systems)
                if hasattr(params, 'l_vcut'):
                    if params.l_vcut is None:
                        params.l_vcut = False
                    if params.l_vcut:
                        raise ValueError('Aperiodic systems require `l_vcut = False`')

                # Makov-Payne (not used for aperiodic systems)
                if getattr(params, 'eps_inf', None) not in [None, 1.0]:
                    raise ValueError('Aperiodic systems should not have a non-zero `eps_inf`')

                # Martyna-Tuckerman
                if self.parameters.mt_correction:
                    if hasattr(params, 'assume_isolated'):
                        if params.assume_isolated is None:
                            params.assume_isolated = 'm-t'
                        if params.assume_isolated not in ['martyna-tuckerman', 'm-t', 'mt']:
                            raise ValueError(
                                f'`assume_isolated = {params.assume_isolated}` is incompatible with `mt_correction = True`')

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

        if self.parameters.from_scratch:
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

        wannier_files_to_link = {}
        dft_ham_files = {}
        if all(self.atoms.pbc):
            # Run PW and Wannierization
            for key in self.calculator_parameters.keys():
                if key.startswith('w90'):
                    self.calculator_parameters[key].write_u_matrices = True
                    self.calculator_parameters[key].write_xyz = True
            wf_workflow = WannierizeWorkflow.fromparent(self, force_nspin2=True, scf_kgrid=self._scf_kgrid)
            wf_workflow.run(subdirectory='wannier')

            # Store the NSCF calculation because it uses nosym, so we can use its outdir for subsequent calculations
            init_pw = [c for c in self.calculations if isinstance(
                c, PWCalculator) and c.parameters.calculation == 'nscf'][-1]

            # Populate a list of files to link to subsequent calculations
            for f in [wf_workflow.outputs.u_matrices_files["occ"],
                      wf_workflow.outputs.hr_files["occ"],
                      wf_workflow.outputs.u_dis_file,
                      wf_workflow.outputs.centers_files["occ"]]:
                assert f is not None
                wannier_files_to_link[f.name] = f

            for f in [wf_workflow.outputs.u_matrices_files["emp"],
                      wf_workflow.outputs.hr_files["emp"],
                      wf_workflow.outputs.centers_files["emp"]]:
                assert f is not None
                wannier_files_to_link[f.parent.prefix + '_emp' + str(f.name)[len(f.parent.prefix):]] = f

            if self.parameters.spin_polarized:
                raise NotImplementedError('Need to adapt the following code for spin-polarized calculations')
            for label in ['occ', 'emp']:
                dft_ham_files[(label, None)] = wf_workflow.outputs.hr_files[label]

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
            pw_workflow.run()

            init_pw = pw_workflow.calculations[0]

        # Convert from wannier to KC
        wann2kc_calc = self.new_calculator('wann2kc')
        self.link(init_pw, init_pw.parameters.outdir, wann2kc_calc, wann2kc_calc.parameters.outdir)
        for dst, f in wannier_files_to_link.items():
            self.link(f.parent, f.name, wann2kc_calc, dst, symlink=True)
        self.run_calculator(wann2kc_calc)

        # Calculate screening parameters
        if self.parameters.calculate_alpha:
            if self.parameters.dfpt_coarse_grid is None:
                screen_wf = ComputeScreeningViaDFPTWorkflow.fromparent(
                    self, wannier_files_to_link=wannier_files_to_link)
                screen_wf.run(subdirectory='screening')
            else:
                self.bands.alphas = coarse_wf.bands.alphas
        else:
            # Load the alphas
            if self.parameters.alpha_from_file:
                self.bands.alphas = [utils.read_alpha_file(Path())]
            else:
                self.bands.alphas = self.parameters.alpha_guess

        # Calculate the Hamiltonian
        if self._perform_ham_calc:
            kc_ham_calc = self.new_calculator('kc_ham', kpts=self.kpoints.path)
            self.link(wann2kc_calc, wann2kc_calc.parameters.outdir / (wann2kc_calc.parameters.prefix + '.save'),
                      kc_ham_calc, kc_ham_calc.parameters.outdir / (kc_ham_calc.parameters.prefix + '.save'), symlink=True)
            self.link(wann2kc_calc, wann2kc_calc.parameters.outdir / (wann2kc_calc.parameters.prefix + '.xml'),
                      kc_ham_calc, kc_ham_calc.parameters.outdir / (kc_ham_calc.parameters.prefix + '.xml'), symlink=True)
            self.link(wann2kc_calc, wann2kc_calc.parameters.outdir / 'kcw', kc_ham_calc,
                      kc_ham_calc.parameters.outdir / 'kcw', recursive_symlink=True)
            for dst, f in wannier_files_to_link.items():
                self.link(f.parent, f.name, kc_ham_calc, dst, symlink=True)
            self.run_calculator(kc_ham_calc)

            # Postprocessing
            if all(self.atoms.pbc) and self.projections and self.kpoints.path is not None \
                    and self.calculator_parameters['ui'].do_smooth_interpolation:
                from koopmans.workflows import UnfoldAndInterpolateWorkflow

                # Assemble the Hamiltonian files required by the UI workflow
                koopmans_ham_files = {("occ", None): FilePointer(kc_ham_calc, f'{kc_ham_calc.parameters.prefix}.kcw_hr_occ.dat'),
                                      ("emp", None): FilePointer(kc_ham_calc, f'{kc_ham_calc.parameters.prefix}.kcw_hr_emp.dat')}
                if not dft_ham_files:
                    raise ValueError(
                        'The DFT Hamiltonian files have not been generated but are required for the UI workflow')
                ui_workflow = UnfoldAndInterpolateWorkflow.fromparent(
                    self, dft_ham_files=dft_ham_files, koopmans_ham_files=koopmans_ham_files)
                ui_workflow.run(subdirectory='postproc')

            # Plotting
            self.plot_bandstructure()

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
        return internal_new_calculator(self, calc_presets, **kwargs)


class ComputeScreeningViaDFPTOutputs(OutputModel):
    pass


class ComputeScreeningViaDFPTWorkflow(Workflow):
    output_model = ComputeScreeningViaDFPTOutputs
    outputs: ComputeScreeningViaDFPTOutputs

    def __init__(self, *args, wannier_files_to_link: Dict[str, FilePointer], **kwargs):
        super().__init__(*args, **kwargs)

        self._wannier_files_to_link = wannier_files_to_link

    def _run(self):
        # Group the bands by spread
        self.bands.assign_groups(sort_by='spread', allow_reassignment=True)

        if len(self.bands.to_solve) == len(self.bands):
            # If there is no orbital grouping, do all orbitals in one calculation

            # 1) Create the calculator
            kc_screen_calc = self.new_calculator('kc_screen')

            # 2) Run the calculator
            self.run_calculator(kc_screen_calc)

            # 3) Store the computed screening parameters
            self.bands.alphas = kc_screen_calc.results['alphas']
        else:
            # If there is orbital grouping, do the orbitals one-by-one

            kc_screen_calcs = []
            # 1) Create the calculators
            for band in self.bands.to_solve:
                kc_screen_calc = self.new_calculator('kc_screen', i_orb=band.index)
                kc_screen_calc.prefix += f'_orbital_{band.index}'
                kc_screen_calcs.append(kc_screen_calc)

            # 2) Run the calculators (possibly in parallel)
            self.run_calculators(kc_screen_calcs)

            # 3) Store the computed screening parameters (accounting for band groupings)
            for band, kc_screen_calc in zip(self.bands.to_solve, kc_screen_calcs):
                for b in self.bands:
                    if b.group == band.group:
                        alpha = kc_screen_calc.results['alphas'][band.spin]
                        b.alpha = alpha[band.spin]

    def new_calculator(self, calc_type: str, *args, **kwargs):
        assert calc_type == 'kc_screen', 'Only the "kc_screen" calculator is supported in DFPTScreeningWorkflow'

        calc = internal_new_calculator(self, calc_type, *args, **kwargs)

        # Link to the most recent wann2kc calculation
        wann2kc_calc = [c for c in self.calculations if isinstance(c, Wann2KCCalculator)][-1]
        self.link(wann2kc_calc, wann2kc_calc.parameters.outdir / (wann2kc_calc.parameters.prefix + '.save'),
                  calc, calc.parameters.outdir / (calc.parameters.prefix + '.save'), symlink=True)
        self.link(wann2kc_calc, wann2kc_calc.parameters.outdir / (wann2kc_calc.parameters.prefix + '.xml'),
                  calc, calc.parameters.outdir / (calc.parameters.prefix + '.xml'), symlink=True)
        self.link(wann2kc_calc, wann2kc_calc.parameters.outdir / 'kcw', calc,
                  calc.parameters.outdir / 'kcw', recursive_symlink=True)

        for dst, f in self._wannier_files_to_link.items():
            self.link(f.parent, f.name, calc, dst, symlink=True)

        return calc


def internal_new_calculator(workflow, calc_presets, **kwargs):
    if calc_presets not in ['kc_ham', 'kc_screen', 'wann2kc']:
        raise ValueError(
            f'Invalid choice `calc_presets={calc_presets}` in `{workflow.__class__.__name__}.new_calculator()`')

    calc = super(workflow.__class__, workflow).new_calculator(calc_presets)

    calc.prefix = calc_presets
    calc.parameters.prefix = 'kc'
    calc.parameters.outdir = 'TMP'
    if all(workflow.atoms.pbc):
        calc.parameters.seedname = [c for c in workflow.calculations if isinstance(c, Wannier90Calculator)][-1].prefix
    calc.parameters.spin_component = 1
    calc.parameters.kcw_at_ks = not all(workflow.atoms.pbc)
    calc.parameters.read_unitary_matrix = all(workflow.atoms.pbc)

    if calc_presets == 'wann2kc':
        pass
    elif calc_presets == 'kc_screen':
        # If eps_inf is not provided in the kc_wann:screen subdictionary but there is a value provided in the
        # workflow parameters, adopt that value
        if workflow.parameters.eps_inf is not None and calc.parameters.eps_inf is None and all(workflow.atoms.pbc):
            calc.parameters.eps_inf = workflow.parameters.eps_inf
    else:
        calc.parameters.do_bands = all(workflow.atoms.pbc)
        calc.alphas = workflow.bands.alphas

    if all(workflow.atoms.pbc):
        nocc = workflow.bands.num(filled=True)
        nemp = workflow.bands.num(filled=False)
        have_empty = (nemp > 0)
        has_disentangle = (workflow.projections.num_bands() != nocc + nemp)
    else:
        nocc = workflow.calculator_parameters['kcp'].nelec // 2
        nemp = workflow.calculator_parameters['pw'].nbnd - nocc
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
