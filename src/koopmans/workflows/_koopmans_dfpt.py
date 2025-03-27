"""The workflow for performing KI and KIPZ calculations with kcw.x."""

from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from koopmans import pseudopotentials, utils
from koopmans.bands import Bands
from koopmans.calculators import (KoopmansHamCalculator, PWCalculator,
                                  Wann2KCCalculator, Wannier90Calculator)
from koopmans.files import File
from koopmans.process_io import IOModel
from koopmans.projections import BlockID
from koopmans.status import Status

from ._dft import DFTPWWorkflow
from ._unfold_and_interp import UnfoldAndInterpolateWorkflow
from ._wannierize import WannierizeWorkflow
from ._workflow import Workflow


class KoopmansDFPTOutputs(IOModel):
    """Pydantic model for the outputs of a `KoopmansDFPTWorkflow`."""

    pass


class KoopmansDFPTWorkflow(Workflow[KoopmansDFPTOutputs]):
    """Workflow for calculating the screening parameters of a system using DFPT."""

    output_model = KoopmansDFPTOutputs  # type: ignore

    def __init__(self, scf_kgrid: Optional[List[int]] = None, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

        # Check the consistency of keywords
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
                                f'`assume_isolated = {params.assume_isolated}` is incompatible with '
                                '`mt_correction = True`')

        # Initialize the bands
        tols: Dict[str, float] = {}
        if self.parameters.spin_polarized:
            nelec = pseudopotentials.nelec_from_pseudos(self.atoms, self.pseudopotentials)
            tot_mag = self.calculator_parameters['pw'].tot_magnetization
            nocc_up = (nelec + tot_mag) // 2
            nocc_dw = (nelec - tot_mag) // 2
            if all(self.atoms.pbc):
                # Using Wannier functions
                ntot_up = self.projections.num_wann(spin='up')
                ntot_dw = self.projections.num_wann(spin='down')
            else:
                # Using KS bands
                ntot_up = self.calculator_parameters['pw'].nbnd
                ntot_dw = self.calculator_parameters['pw'].nbnd
            nemp_up = ntot_up - nocc_up
            nemp_dw = ntot_dw - nocc_dw
            filling = [[True for _ in range(nocc_up)] + [False for _ in range(nemp_up)],
                       [True for _ in range(nocc_dw)] + [False for _ in range(nemp_dw)]]
            if self.parameters.orbital_groups is None:
                self.parameters.orbital_groups = [
                    list(range(nocc_up + nemp_up)), list(range(nocc_dw + nemp_dw))]
            for key in ['self_hartree', 'spread']:
                val = self.parameters.get(f'orbital_groups_{key}_tol', None)
                if val is not None:
                    tols[key] = val
            self.bands = Bands(n_bands=[len(f) for f in filling], n_spin=2,
                               spin_polarized=self.parameters.spin_polarized,
                               filling=filling, groups=self.parameters.orbital_groups, tolerances=tols)
        else:
            nocc = pseudopotentials.nelec_from_pseudos(self.atoms, self.pseudopotentials) // 2
            if all(self.atoms.pbc):
                exclude_bands = self.calculator_parameters['w90'].get(
                    'exclude_bands', [])
                nocc -= len(exclude_bands)
                ntot = self.projections.num_wann()
            else:
                ntot = self.calculator_parameters['pw'].nbnd
            nemp = ntot - nocc
            filling = [[True] * nocc + [False] * nemp]
            if self.parameters.orbital_groups is None:
                self.parameters.orbital_groups = [list(range(nocc + nemp))]
            for key in ['self_hartree', 'spread']:
                val = self.parameters.get(f'orbital_groups_{key}_tol', None)
                if val is not None:
                    tols[key] = val
            self.bands = Bands(n_bands=nocc + nemp, filling=filling,
                               groups=self.parameters.orbital_groups, tolerances=tols)

        # Populating kpoints if absent
        if not all(self.atoms.pbc):
            for key in ['pw', 'kcw_screen']:
                self.calculator_parameters[key].kpts = [1, 1, 1]

        self._perform_ham_calc: bool = True

        self._scf_kgrid = scf_kgrid

    def _run(self):
        """Run the workflow."""
        if self.parameters.dfpt_coarse_grid is not None:
            coarse_wf = self.__class__.fromparent(
                self, scf_kgrid=self.kpoints.grid)
            coarse_wf.parameters.dfpt_coarse_grid = None
            coarse_wf.kpoints.grid = self.parameters.dfpt_coarse_grid
            coarse_wf._perform_ham_calc = False
            coarse_wf.name += '-coarse'
            coarse_wf.proceed()
            if coarse_wf.status != Status.COMPLETED:
                return

        wannier_files_to_link_by_spin = []
        dft_ham_files = {}
        if all(self.atoms.pbc):
            # Run PW and Wannierization
            for key in self.calculator_parameters.keys():
                if key.startswith('w90'):
                    self.calculator_parameters[key].write_u_matrices = True
                    self.calculator_parameters[key].write_xyz = True
            wf_workflow = WannierizeWorkflow.fromparent(self, force_nspin2=True, scf_kgrid=self._scf_kgrid)
            wf_workflow.proceed()
            if wf_workflow.status != Status.COMPLETED:
                return

            # Store the NSCF calculation because it uses nosym, so we can use its outdir for subsequent calculations
            init_pw = [c for c in self.calculations if isinstance(
                c, PWCalculator) and c.parameters.calculation == 'nscf'][-1]

            # Populate a list of files to link to subsequent calculations
            spins = ['up', 'down'] if self.parameters.spin_polarized else [None]
            for spin in spins:
                wannier_files_to_link_by_spin.append({})
                for filled in [True, False]:
                    block_id = BlockID(filled=filled, spin=spin)
                    for f in [wf_workflow.outputs.u_matrices_files[block_id],
                              wf_workflow.outputs.hr_files[block_id],
                              wf_workflow.outputs.centers_files[block_id]]:
                        assert f is not None

                        if filled:
                            wannier_files_to_link_by_spin[-1][f.name] = f
                        else:
                            wannier_files_to_link_by_spin[-1]['wannier90_emp' + str(f.name)[9:]] = f

                    dft_ham_files[block_id] = wf_workflow.outputs.hr_files[block_id]

                # Empty blocks might also have a disentanglement file than we need to copy
                if wf_workflow.outputs.u_dis_files[block_id] is not None:
                    wannier_files_to_link_by_spin[-1]['wannier90_emp_u_dis.mat'] = \
                        wf_workflow.outputs.u_dis_files[block_id]

        else:
            # Run PW
            self.print('Initialization of density and variational orbitals', style='heading')

            # Create the workflow
            pw_workflow = DFTPWWorkflow.fromparent(self)

            # Update settings
            pw_params = pw_workflow.calculator_parameters['pw']
            pw_params.nspin = 2

            # Run the subworkflow
            pw_workflow.proceed()
            if pw_workflow.status != Status.COMPLETED:
                return

            init_pw = pw_workflow.calculations[0]
            wannier_files_to_link_by_spin = [{}, {}] if self.parameters.spin_polarized else [{}]

        spin_components = [1, 2] if self.parameters.spin_polarized else [1]

        # Convert from wannier to KC
        wann2kc_calcs = []
        for spin_component, wannier_files_to_link in zip(spin_components, wannier_files_to_link_by_spin):
            spin_suffix = f'_spin_{spin_component}' if self.parameters.spin_polarized else ''

            wann2kc_calc = self.new_calculator('kcw_wannier', spin_component=spin_component)
            wann2kc_calc.link(File(init_pw, init_pw.parameters.outdir), wann2kc_calc.parameters.outdir)

            for dst, f in wannier_files_to_link.items():
                wann2kc_calc.link(f, dst, symlink=True)

            wann2kc_calc.prefix += spin_suffix
            wann2kc_calcs.append(wann2kc_calc)
        status = self.run_steps(wann2kc_calcs)
        if status != Status.COMPLETED:
            return

        # Calculate screening parameters
        screen_wfs = []
        for spin_component, wannier_files_to_link, wann2kc_calc in zip(spin_components,
                                                                       wannier_files_to_link_by_spin,
                                                                       wann2kc_calcs):
            spin_suffix = f'_spin_{spin_component}' if self.parameters.spin_polarized else ''
            if self.parameters.calculate_alpha:
                if self.parameters.dfpt_coarse_grid is None:
                    screen_wf = ComputeScreeningViaDFPTWorkflow.fromparent(
                        self, wannier_files_to_link=wannier_files_to_link, spin_component=spin_component)
                    screen_wf.name += spin_suffix
                    screen_wfs.append(screen_wf)
            else:
                # Load the alphas
                if self.parameters.alpha_from_file:
                    self.bands.alphas = [utils.read_alpha_file(Path())]
                else:
                    self.bands.alphas = self.parameters.alpha_guess

        for wf in screen_wfs:
            wf.proceed()

        if any([wf.status != Status.COMPLETED for wf in screen_wfs]):
            return

        # Calculate the Hamiltonian
        kc_ham_calcs = []
        for spin_component, wannier_files_to_link, wann2kc_calc in zip(spin_components, wannier_files_to_link_by_spin,
                                                                       wann2kc_calcs):
            spin_suffix = f'_spin_{spin_component}' if self.parameters.spin_polarized else ''
            if self._perform_ham_calc:
                kc_ham_calc = self.new_calculator('kcw_ham', kpts=self.kpoints.path, spin_component=spin_component)
                kc_ham_calc.link(
                    File(wann2kc_calc, wann2kc_calc.parameters.outdir / (wann2kc_calc.parameters.prefix + '.save')),
                    kc_ham_calc.parameters.outdir / (kc_ham_calc.parameters.prefix + '.save'),
                    symlink=True
                )
                kc_ham_calc.link(
                    File(wann2kc_calc, wann2kc_calc.parameters.outdir / (wann2kc_calc.parameters.prefix + '.xml')),
                    kc_ham_calc.parameters.outdir / (kc_ham_calc.parameters.prefix + '.xml'),
                    symlink=True
                )
                kc_ham_calc.link(File(wann2kc_calc, wann2kc_calc.parameters.outdir / 'kcw'),
                                 kc_ham_calc.parameters.outdir / 'kcw', recursive_symlink=True)
                for dst, f in wannier_files_to_link.items():
                    kc_ham_calc.link(f, dst, symlink=True)
                kc_ham_calc.prefix += spin_suffix
                kc_ham_calcs.append(kc_ham_calc)

        status = self.run_steps(kc_ham_calcs)
        if status != Status.COMPLETED:
            return

        # Postprocessing
        if self._perform_ham_calc:
            if all(self.atoms.pbc) and self.projections and self.kpoints.path is not None \
                    and self.calculator_parameters['ui'].do_smooth_interpolation:

                # Assemble the Hamiltonian files required by the UI workflow
                if self.parameters.spin_polarized:
                    koopmans_ham_files = {
                        BlockID(filled=True, spin="up"): File(kc_ham_calc,
                                                              f'{kc_ham_calcs[0].parameters.prefix}.kcw_hr_occ.dat'),
                        BlockID(filled=False, spin="up"): File(kc_ham_calc,
                                                               f'{kc_ham_calcs[0].parameters.prefix}.kcw_hr_emp.dat'),
                        BlockID(filled=True, spin="down"): File(kc_ham_calc,
                                                                f'{kc_ham_calcs[1].parameters.prefix}.kcw_hr_occ.dat'),
                        BlockID(filled=False, spin="down"): File(kc_ham_calc,
                                                                 f'{kc_ham_calcs[1].parameters.prefix}.kcw_hr_emp.dat')
                    }
                else:
                    koopmans_ham_files = {
                        BlockID(filled=True): File(kc_ham_calc,
                                                   f'{kc_ham_calc.parameters.prefix}.kcw_hr_occ.dat'),
                        BlockID(filled=False): File(kc_ham_calc,
                                                    f'{kc_ham_calc.parameters.prefix}.kcw_hr_emp.dat')
                    }
                if not dft_ham_files:
                    raise ValueError(
                        'The DFT Hamiltonian files have not been generated but are required for the UI workflow')
                ui_workflow = UnfoldAndInterpolateWorkflow.fromparent(
                    self, dft_ham_files=dft_ham_files, koopmans_ham_files=koopmans_ham_files)

                ui_workflow.proceed()
                if ui_workflow.status != Status.COMPLETED:
                    return

        # Plotting
        if self._perform_ham_calc:
            self.plot_bandstructure()

        self.status = Status.COMPLETED
        return

    def plot_bandstructure(self):
        """Plot the band structure."""
        if not all(self.atoms.pbc):
            return

        # Identify the relevant calculators
        n_calc = 2 if self.parameters.spin_polarized else 1
        kc_ham_calcs = [c for c in self.calculations if isinstance(c, KoopmansHamCalculator)][-n_calc:]

        # Plot the bandstructure if the band path has been specified
        bandstructures = [c.results['band structure'] for c in kc_ham_calcs]
        ref = max([b.reference for b in bandstructures])
        bs_to_plot = bandstructures[0].subtract_reference(ref)
        if not bs_to_plot.path.path:
            return

        if self.parameters.spin_polarized:
            bs_to_plot._energies = np.append(bs_to_plot._energies, bandstructures[1]._energies - ref, axis=0)

        super().plot_bandstructure(bs_to_plot)

    def new_calculator(self, calc_presets, **kwargs):
        """Create a new calculator for the workflow."""
        return internal_new_calculator(self, calc_presets, **kwargs)


class ComputeScreeningViaDFPTOutputs(IOModel):
    """Pydantic model for the outputs of the `ComputeScreeningViaDFPTWorkflow`."""

    pass


class ComputeScreeningViaDFPTWorkflow(Workflow[ComputeScreeningViaDFPTOutputs]):
    """Workflow that computes the screening parameters of a system using DFPT."""

    output_model = ComputeScreeningViaDFPTOutputs

    def __init__(self, *args, spin_component: int, wannier_files_to_link: Dict[str, File], **kwargs):
        super().__init__(*args, **kwargs)

        self._wannier_files_to_link = wannier_files_to_link
        self._spin_component = spin_component

    def _run(self):
        # Group the bands by spread
        self.bands.assign_groups(sort_by='spread', allow_reassignment=True)

        if len(self.bands.to_solve) == len(self.bands):
            # If there is no orbital grouping, do all orbitals in one calculation

            # 1) Create the calculator
            kc_screen_calc = self.new_calculator('kcw_screen', spin_component=self._spin_component)

            # 2) Run the calculator
            status = self.run_steps(kc_screen_calc)
            if status != Status.COMPLETED:
                return

            # 3) Store the computed screening parameters
            self.bands.alphas = kc_screen_calc.results['alphas']
        else:
            # If there is orbital grouping, do the orbitals one-by-one
            assert self.bands is not None
            bands_to_solve = [b for b in self.bands.to_solve if b.spin == self._spin_component - 1]

            # 1) Create the calculators
            kc_screen_calcs = []
            for band in bands_to_solve:
                kc_screen_calc = self.new_calculator('kcw_screen', i_orb=band.index,
                                                     spin_component=self._spin_component)
                kc_screen_calc.prefix += f'_orbital_{band.index}_spin_{self._spin_component}'
                kc_screen_calcs.append(kc_screen_calc)

            # 2) Run the calculators (possibly in parallel)
            status = self.run_steps(kc_screen_calcs)
            if status != Status.COMPLETED:
                return

            # 3) Store the computed screening parameters (accounting for band groupings)
            for band, kc_screen_calc in zip(bands_to_solve, kc_screen_calcs):
                for b in self.bands:
                    if b.group == band.group:
                        [[alpha]] = kc_screen_calc.results['alphas']
                        b.alpha = alpha

        self.status = Status.COMPLETED

        return

    def new_calculator(self, calc_type: str, *args, **kwargs):
        """Create a new calculator for the workflow."""
        assert calc_type == 'kcw_screen', 'Only the "kcw_screen" calculator is supported in DFPTScreeningWorkflow'

        calc = internal_new_calculator(self, calc_type, *args, **kwargs)

        # Link to the most recent wann2kc calculation
        wann2kc_calc = [c for c in self.calculations if isinstance(c, Wann2KCCalculator)][-1]
        calc.link(File(wann2kc_calc, wann2kc_calc.parameters.outdir / (wann2kc_calc.parameters.prefix + '.save')),
                  calc.parameters.outdir / (calc.parameters.prefix + '.save'), symlink=True)
        calc.link(File(wann2kc_calc, wann2kc_calc.parameters.outdir / (wann2kc_calc.parameters.prefix + '.xml')),
                  calc.parameters.outdir / (calc.parameters.prefix + '.xml'), symlink=True)
        calc.link(File(wann2kc_calc, wann2kc_calc.parameters.outdir / 'kcw'),
                  calc.parameters.outdir / 'kcw', recursive_symlink=True)

        for dst, f in self._wannier_files_to_link.items():
            calc.link(f, dst, symlink=True)

        return calc


def internal_new_calculator(workflow, calc_presets, **kwargs):
    """Create a new calculator for the workflow."""
    if calc_presets not in ['kcw_ham', 'kcw_screen', 'kcw_wannier']:
        raise ValueError(
            f'Invalid choice `calc_presets={calc_presets}` in `{workflow.__class__.__name__}.new_calculator()`')

    calc = super(workflow.__class__, workflow).new_calculator(calc_presets)

    calc.prefix = calc_presets
    calc.parameters.prefix = 'kc'
    calc.parameters.outdir = 'TMP'
    if all(workflow.atoms.pbc):
        calc.parameters.seedname = [c for c in workflow.calculations if isinstance(c, Wannier90Calculator)][-1].prefix
    calc.parameters.spin_component = kwargs['spin_component'] if 'spin_component' in kwargs else 1
    calc.parameters.kcw_at_ks = not all(workflow.atoms.pbc)
    calc.parameters.read_unitary_matrix = all(workflow.atoms.pbc)

    if calc_presets == 'kcw_wannier':
        pass
    elif calc_presets == 'kcw_screen':
        # If eps_inf is not provided in the kc_wann:screen subdictionary but there is a value provided in the
        # workflow parameters, adopt that value
        if workflow.parameters.eps_inf is not None and calc.parameters.eps_inf is None and all(workflow.atoms.pbc):
            calc.parameters.eps_inf = workflow.parameters.eps_inf
    else:
        calc.parameters.do_bands = all(workflow.atoms.pbc)
        if not workflow.parameters.spin_polarized and len(workflow.bands.alphas) != 1:
            raise ValueError('The list of screening parameters should be length 1 for spin-unpolarized calculations')
        calc.alphas = workflow.bands.alphas[calc.parameters.spin_component - 1]

    if all(workflow.atoms.pbc):
        if workflow.parameters.spin_polarized:
            if calc.parameters.spin_component == 1:
                nocc = workflow.calculator_parameters['kcp'].nelup
                nemp = workflow.projections.num_wann(spin='up') - nocc
            else:
                nocc = workflow.calculator_parameters['kcp'].neldw
                nemp = workflow.projections.num_wann(spin='down') - nocc
        else:
            nocc = workflow.bands.num(filled=True)
            nemp = workflow.bands.num(filled=False)
        nemp_pw = workflow.calculator_parameters['pw'].nbnd - nocc
        have_empty = (nemp > 0)
        has_disentangle = (nemp != nemp_pw)
    else:
        if workflow.parameters.spin_polarized:
            if calc.parameters.spin_component == 1:
                nocc = workflow.calculator_parameters['kcp'].nelup
            else:
                nocc = workflow.calculator_parameters['kcp'].neldw
        else:
            nocc = workflow.calculator_parameters['kcp'].nelec // 2
        nemp = workflow.calculator_parameters['pw'].nbnd - nocc
        have_empty = (nemp > 0)
        has_disentangle = False
    calc.parameters.num_wann_occ = nocc
    calc.parameters.num_wann_emp = nemp
    calc.parameters.have_empty = have_empty
    calc.parameters.has_disentangle = has_disentangle

    # Ensure that any additional calculator keywords passed as kwargs are still that value
    for k, v in kwargs.items():
        setattr(calc.parameters, k, v)

    return calc
