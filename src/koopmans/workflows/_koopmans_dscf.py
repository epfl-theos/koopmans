"""

Workflow module for koopmans, containing the workflow for performing KI and KIPZ calculations

Written by Edward Linscott Jan 2020
Split off from workflow.py Oct 2020

"""

import logging
import shutil
from pathlib import Path
from typing import Dict, Generator, List, Mapping, Optional, Tuple

import numpy as np
from ase_koopmans.dft import DOS
from pydantic import ConfigDict

from koopmans import calculators, utils
from koopmans.bands import Band, Bands
from koopmans.files import File
from koopmans.process_io import IOModel
from koopmans.processes.koopmans_cp import (ConvertFilesFromSpin1To2,
                                            ConvertFilesFromSpin2To1)
from koopmans.projections import BlockID
from koopmans.settings import KoopmansCPSettingsDict
from koopmans.status import Status

from ._folding import FoldToSupercellWorkflow
from ._koopmans_cp_with_spin_swap import KoopmansCPWithSpinSwapWorkflow
from ._ml import PowerSpectrumDecompositionWorkflow
from ._unfold_and_interp import UnfoldAndInterpolateWorkflow
from ._wannierize import WannierizeWorkflow
from ._workflow import Workflow, spin_symmetrize

logger = logging.getLogger(__name__)


class KoopmansDSCFOutputs(IOModel):
    '''
    Outputs for the KoopmansDSCFWorkflow
    '''
    variational_orbital_files: Dict[str, File]
    final_calc: calculators.KoopmansCPCalculator
    wannier_hamiltonian_files: Dict[BlockID, File] | None = None
    smooth_dft_ham_files: Dict[BlockID, File] | None = None
    model_config = ConfigDict(arbitrary_types_allowed=True)


class KoopmansDSCFWorkflow(Workflow[KoopmansDSCFOutputs]):

    output_model = KoopmansDSCFOutputs

    def __init__(self, *args,
                 initial_variational_orbital_files: Dict[str, File] | None = None,
                 previous_cp_calc: calculators.KoopmansCPCalculator | None = None,
                 smooth_dft_ham_files: Dict[BlockID, File] | None = None,
                 precomputed_descriptors: List[File] | None = None, **kwargs) -> None:
        super().__init__(*args, **kwargs)

        # The following two additional keywords allow for some tweaking of the workflow when running a singlepoint
        # workflow with functional == 'all'

        # If the higher-resolution DFT calculations have already been computed, the smooth DFT Hamiltonian can be provided
        # via this keyword to avoid redoing these calculations unnecessarily
        self._smooth_dft_ham_files = smooth_dft_ham_files

        # For the KIPZ calculation we restart from the old KI calculation
        self._initial_variational_orbital_files = initial_variational_orbital_files

        # For workflows where we have already performed initialization elsewhere, we can restart using that directory
        self._previous_cp_calc = previous_cp_calc

        # For workflows where we have previously computed descriptors, we can provide them here
        self._precomputed_descriptors = precomputed_descriptors

        # If periodic, convert the kcp calculation into a Γ-only supercell calculation
        kcp_params = self.calculator_parameters['kcp']
        if all(self.atoms.pbc):
            spins: List[Optional[str]]
            if self.parameters.spin_polarized:
                spins = ['up', 'down']
                nelecs = [kcp_params.nelup, kcp_params.neldw]
            else:
                spins = [None]
                nelecs = [kcp_params.nelec // 2]

            for spin, nelec in zip(spins, nelecs):
                # Check that we have wannierized every filled orbital
                nbands_occ = nelec
                if self.projections:
                    label = 'w90'
                    if spin:
                        label += f'_{spin}'
                    nbands_excl = len(self.calculator_parameters[label].get('exclude_bands', []))
                    if nbands_excl > 0:
                        raise ValueError('Excluding bands is incompatible with `method == "dscf"`. Please provide '
                                         'projections for every band and remove the `exclude_bands` `Wannier90` keyword.')

                    nwann = self.projections.num_wann(spin=spin)

                    if nwann < nelec:
                        raise ValueError('You have configured this calculation to only wannierize a subset of the '
                                         'occupied bands:\n'
                                         f' number of occupied bands = {nbands_occ}\n'
                                         f' number of Wannier functions = {nwann}\n'
                                         'This is incompatible with the subsequent Koopmans '
                                         'calculation.\nPlease modify the `Wannier90` settings in order to wannierize '
                                         'all of the occupied bands.')

                    nbands_emp = nwann - nbands_occ
                else:
                    nbands_emp = self.calculator_parameters['pw'].nbnd - nbands_occ

                # Check the number of empty states has been correctly configured
                spin_info = f'spin {spin} ' if self.parameters.spin_polarized else ''
                if kcp_params.nbnd is None:
                    if nbands_emp != 0:
                        kcp_params.nbnd = nbands_occ + nbands_emp
                elif nbands_occ > kcp_params.nbnd:
                    raise ValueError(f'The value you have provided for `nbnd` is less than the number of {spin_info}'
                                     f'electrons. Please increase `nbnd` to at least {nbands_occ}')
                elif kcp_params.nbnd != nbands_occ + nbands_emp:
                    raise ValueError(f'The number of {spin_info}empty states are inconsistent:\n'
                                     f' number of empty bands = {kcp_params.nbnd - nbands_occ}\n'
                                     f' number of empty Wannier functions = {nbands_emp}\n'
                                     'If you have provided `nbnd` explicitly to the `kcp` calculator, check that it '
                                     'matches with the number of empty projections/bands in your system.')

            # Populating self.parameters.orbital_groups if needed
            # N.B. self.bands.groups is guaranteed to be 2 x num_wann, but self.parameters.orbital_groups
            # is either 1- or 2- long, depending on if we are spin-polarized or not
            if self.parameters.orbital_groups is None:
                orbital_groups: List[List[int]] = []
                i_start = 0
                for nelec in nelecs:
                    i_end = i_start + kcp_params.get('nbnd', nelec) - 1
                    orbital_groups.append(list(range(i_start, i_end + 1)))
                    i_start = i_end + 1
                self.parameters.orbital_groups = orbital_groups

            if not self.kpoints.gamma_only:
                # Update the KCP settings to correspond to a supercell (leaving self.atoms unchanged for the moment)
                self.convert_kcp_to_supercell()

                # Expanding self.parameters.orbital_groups to account for the supercell, grouping equivalent wannier
                # functions together
                for i_spin, nelec in enumerate(nelecs):
                    assert self.kpoints.grid is not None
                    self.parameters.orbital_groups[i_spin] = [i for _ in range(np.prod(self.kpoints.grid))
                                                              for i in self.parameters.orbital_groups[i_spin][:nelec]] \
                        + [i for _ in range(np.prod(self.kpoints.grid)) for i in
                           self.parameters.orbital_groups[i_spin][nelec:]]

        # Check the shape of self.parameters.orbital_groups is as expected
        if self.parameters.spin_polarized:
            target_length = 2
        else:
            target_length = 1
        if self.parameters.orbital_groups is not None:
            assert len(self.parameters.orbital_groups) == target_length

        # Constructing the arrays required to initialize a Bands object
        if self.parameters.spin_polarized:
            if 'nbnd' in kcp_params:
                n_emp_up = kcp_params.nbnd - kcp_params.nelup
                n_emp_dw = kcp_params.nbnd - kcp_params.neldw
            else:
                n_emp_up = 0
                n_emp_dw = 0
            filling = [[True for _ in range(kcp_params.nelup)] + [False for _ in range(n_emp_up)],
                       [True for _ in range(kcp_params.neldw)] + [False for _ in range(n_emp_dw)]]
            groups = self.parameters.orbital_groups
        else:
            if 'nbnd' in kcp_params:
                n_emp = kcp_params.nbnd - kcp_params.nelec // 2
            else:
                n_emp = 0
            filling = [[True for _ in range(kcp_params.nelec // 2)]
                       + [False for _ in range(n_emp)] for _ in range(2)]
            # self.parameters.orbital_groups does not have a spin index
            if self.parameters.orbital_groups is None:
                groups = None
            else:
                groups = [self.parameters.orbital_groups[0] for _ in range(2)]

        # Checking groups and filling are the same dimensions
        if groups is not None:
            for g, f in zip(groups, filling):
                assert len(g) == len(f), 'orbital_groups is the wrong dimension; its length should match the number ' \
                    'of bands'

        # Initialize the bands object
        tols: Dict[str, float] = {}
        for key in ['self_hartree', 'spread']:
            val = self.parameters.get(f'orbital_groups_{key}_tol', None)
            if val is not None:
                tols[key] = val
        self.bands = Bands(n_bands=[len(f) for f in filling], n_spin=2, spin_polarized=self.parameters.spin_polarized,
                           filling=filling, groups=groups, tolerances=tols)

        if self.parameters.alpha_from_file:
            # Reading alpha values from file
            self.bands.alphas = self.read_alphas_from_file()
        else:
            # Initializing alpha with a guess
            self.bands.alphas = self.parameters.alpha_guess

        # Raise errors if any UI keywords are provided but will be overwritten by the workflow
        for ui_keyword in ['kc_ham_file', 'w90_seedname', 'dft_ham_file', 'dft_smooth_ham_file']:
            for ui_kind in ['occ', 'emp']:
                value = getattr(self.calculator_parameters[f'ui_{ui_kind}'], ui_keyword)
                [default_value] = [s.default for s in self.calculator_parameters['ui'].settings if s.name == ui_keyword]
                if value != default_value:
                    raise ValueError(f'UI keyword `{ui_keyword}` has been set in the input file, but this will be '
                                     'automatically set by the Koopmans workflow. Remove this keyword from the input '
                                     'file')

        # Check self.init_empty_orbitals
        if self.parameters.init_empty_orbitals != self.parameters.init_orbitals:
            raise NotImplementedError(f'The combination `init_orbitals` = {self.parameters.init_orbitals} '
                                      f'and `init_empty_orbitals` = {self.parameters.init_empty_orbitals} '
                                      'has not yet been implemented')

    def convert_kcp_to_supercell(self):
        # Multiply all extensive KCP settings by the appropriate prefactor
        assert self.kpoints.grid is not None
        prefactor = np.prod(self.kpoints.grid)
        for attr in ['nelec', 'nelup', 'neldw', 'nbnd', 'conv_thr', 'esic_conv_thr', 'tot_charge', 'tot_magnetization']:
            value = getattr(self.calculator_parameters['kcp'], attr, None)
            if value is not None:
                setattr(self.calculator_parameters['kcp'], attr, prefactor * value)

    def read_alphas_from_file(self, directory: Path = Path()):
        '''
        This routine reads in the contents of file_alpharef.txt and file_alpharef_empty.txt

        Since utils.read_alpha_file provides a flattened list of alphas so we must convert this
        to a nested list using convert_flat_alphas_for_kcp()
        '''

        flat_alphas = utils.read_alpha_file(self)
        params = self.calculator_parameters['kcp']
        assert isinstance(params, KoopmansCPSettingsDict)
        alphas = calculators.convert_flat_alphas_for_kcp(flat_alphas, params)

        if self.parameters.spin_polarized:
            raise NotImplementedError('Need to check implementation')

        return alphas

    def _run(self) -> None:
        '''
        This function runs a KI/pKIPZ/KIPZ workflow from start to finish
        '''

        init_wf: Optional[InitializationWorkflow] = None
        if self._initial_variational_orbital_files is None:
            init_wf = InitializationWorkflow.fromparent(self)
            init_wf.run()
            if init_wf.status != Status.COMPLETED:
                return
            self._initial_variational_orbital_files = init_wf.outputs.variational_orbital_files
            initial_cp_calculation = init_wf.outputs.final_calc
        else:
            assert self._previous_cp_calc is not None
            initial_cp_calculation = self._previous_cp_calc

        self.primitive_to_supercell()

        if self.parameters.calculate_alpha:
            screening_wf = CalculateScreeningViaDSCF.fromparent(self, initial_variational_orbital_files=self._initial_variational_orbital_files,
                                                                initial_cp_calculation=initial_cp_calculation,
                                                                precomputed_descriptors=self._precomputed_descriptors)
            screening_wf.run()
            if screening_wf.status != Status.COMPLETED:
                return

            # Store the files which will be needed for the final calculation
            n_electron_restart_dir = screening_wf.outputs.n_electron_restart_dir

        else:
            self.print('Skipping calculation of screening parameters', end='')
            assert self.bands is not None
            if len(self.bands.alpha_history()) == 0:
                self.print('; reading values from file')
                self.bands.alphas = self.read_alphas_from_file()
            print_alpha_history(self)

            # In this case the final calculation will restart from the initialization calculations
            if self._previous_cp_calc is None:
                assert init_wf is not None
                n_electron_restart_dir = init_wf.outputs.final_calc.write_directory
            else:
                n_electron_restart_dir = self._previous_cp_calc.write_directory

        # Final calculations
        if self.ml.test:
            # If we are testing the model, we want to run both with and without the ML model
            use_mls = [True, False]
        elif self.ml.predict:
            # Use the ML model
            use_mls = [True]
        else:
            # Don't use the ML model
            use_mls = [False]

        for use_ml in use_mls:
            if self.parameters.functional == 'pkipz':
                if self._previous_cp_calc is None:
                    final_calc_types = ['ki', 'pkipz']
                else:
                    final_calc_types = ['pkipz']
            else:
                assert isinstance(self.parameters.functional, str)
                final_calc_types = [self.parameters.functional]

            for final_calc_type in final_calc_types:

                final_calc_type += '_final'

                assert self.bands is not None
                if use_ml:
                    alphas = self.bands.predicted_alphas
                else:
                    alphas = self.bands.alphas

                calc = internal_new_kcp_calculator(self, final_calc_type, write_hr=True, alphas=alphas)

                if self.parameters.functional == 'ki' and self.parameters.init_orbitals in ['mlwfs', 'projwfs'] \
                        and not self.parameters.calculate_alpha:
                    calc.parameters.restart_from_wannier_pwscf = True

                if use_ml:
                    calc.prefix += '_ml'

                calc.link(n_electron_restart_dir, calc.read_directory, recursive_symlink=True)
                status = self.run_steps(calc)
                if status != Status.COMPLETED:
                    return

        final_calc = calc
        variational_orbital_files = {f: final_calc.read_directory / 'K00001' / f
                                     for f in ['evc01.dat', 'evc02.dat', 'evc0_empty1.dat', 'evc0_empty2.dat']}

        # Postprocessing
        smooth_dft_ham_files: Dict[BlockID, File] | None = None
        if all(self.atoms.pbc):
            if self.parameters.calculate_bands in [None, True] and self.projections and self.kpoints.path is not None:
                # Calculate interpolated band structure and DOS with UI
                final_koopmans_calc = self.calculations[-1]
                koopmans_ham_files: Dict[BlockID, File]
                if self.parameters.spin_polarized:
                    koopmans_ham_files = {BlockID(filled=True, spin="up"): File(final_koopmans_calc, Path('ham_occ_1.dat')),
                                          BlockID(filled=False, spin="up"): File(final_koopmans_calc, Path('ham_emp_1.dat')),
                                          BlockID(filled=True, spin="down"): File(final_koopmans_calc, Path('ham_occ_2.dat')),
                                          BlockID(filled=False, spin="down"): File(final_koopmans_calc, Path('ham_emp_2.dat'))}
                else:
                    koopmans_ham_files = {BlockID(filled=True): File(final_koopmans_calc, Path('ham_occ_1.dat')),
                                          BlockID(filled=False): File(final_koopmans_calc, Path('ham_emp_1.dat'))}
                assert init_wf is not None
                dft_ham_files = init_wf.outputs.wannier_hamiltonian_files
                ui_workflow = UnfoldAndInterpolateWorkflow.fromparent(
                    self,
                    koopmans_ham_files=koopmans_ham_files,
                    dft_ham_files=dft_ham_files,
                    smooth_dft_ham_files=self._smooth_dft_ham_files)
                ui_workflow.run()
                if ui_workflow.status != Status.COMPLETED:
                    return
                smooth_dft_ham_files = ui_workflow.outputs.smooth_dft_ham_files
            else:
                # Generate the DOS only
                dos = DOS(self.calculations[-1], width=self.plotting.degauss, npts=self.plotting.nstep + 1)
                self.calculations[-1].results['dos'] = dos

        self.outputs = self.output_model(variational_orbital_files=variational_orbital_files, final_calc=final_calc,
                                         smooth_dft_ham_files=smooth_dft_ham_files)

        self.status = Status.COMPLETED

        return

    def _overwrite_canonical_with_variational_orbitals(self, calc: calculators.KoopmansCPCalculator) -> None:
        self.print('Overwriting the variational orbitals with Kohn-Sham orbitals')
        savedir = calc.parameters.write_directory / f'{calc.parameters.prefix}_{calc.parameters.ndw}.save/K00001'
        for ispin in range(2):
            raise NotImplementedError('Need to replace this shutil call')
            shutil.copy(savedir / f'evc{ispin + 1}.dat', savedir / f'evc0{ispin + 1}.dat')
            if calc.has_empty_states(ispin):
                shutil.copy(savedir / f'evc_empty{ispin + 1}.dat', savedir / f'evc0_empty{ispin + 1}.dat')


class CalculateScreeningViaDSCFOutput(IOModel):
    n_electron_restart_dir: File
    model_config = ConfigDict(arbitrary_types_allowed=True)


class CalculateScreeningViaDSCF(Workflow[CalculateScreeningViaDSCFOutput]):

    output_model = CalculateScreeningViaDSCFOutput

    def __init__(self, *args, initial_variational_orbital_files: Dict[str, Tuple[Workflow, str]],
                 initial_cp_calculation: calculators.KoopmansCPCalculator,
                 precomputed_descriptors: List[File] | None = None, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._initial_variational_orbital_files = initial_variational_orbital_files
        self._initial_cp_calculation = initial_cp_calculation
        self._precomputed_descriptors = precomputed_descriptors

    def _run(self) -> None:
        converged = False
        i_sc = 0

        alpha_indep_calcs: List[calculators.KoopmansCPCalculator] = []

        variational_orbital_files = self._initial_variational_orbital_files
        n_electron_calc = self._initial_cp_calculation
        dummy_outdirs: Dict[Tuple[int, int], File | None] = {}

        assert isinstance(self.parameters.alpha_numsteps, int)
        while not converged and i_sc < self.parameters.alpha_numsteps:
            i_sc += 1

            iteration_wf = DeltaSCFIterationWorkflow.fromparent(self, variational_orbital_files=variational_orbital_files,
                                                                previous_n_electron_calculation=n_electron_calc,
                                                                precomputed_descriptors=self._precomputed_descriptors,
                                                                dummy_outdirs=dummy_outdirs,
                                                                i_sc=i_sc, alpha_indep_calcs=alpha_indep_calcs)
            iteration_wf.name = f'Iteration {i_sc}'

            if i_sc == 1:

                # For the first iteration, the spin contamination has already been addressed during the initialization
                iteration_wf.parameters.fix_spin_contamination = False

            iteration_wf.run()
            if iteration_wf.status != Status.COMPLETED:
                return

            converged = iteration_wf.outputs.converged or self.ml.predict

            assert self.bands is not None
            if self.parameters.functional == 'ki' and self.bands.num(filled=False) == 0:
                # For this case the screening parameters are guaranteed to converge instantly
                if self.parameters.alpha_numsteps == 1:
                    # Print the "converged" message rather than the "determined but not necessarily converged" message
                    converged = True
                else:
                    # Do the subsequent loop
                    utils.warn('The screening parameters for a KI calculation with no empty states will converge '
                               'instantly; to save computational time set `alpha_numsteps == 1`')

            parent = iteration_wf.outputs.n_electron_restart_dir.parent_process
            assert isinstance(parent, calculators.KoopmansCPCalculator)
            n_electron_calc = parent
            dummy_outdirs = iteration_wf.outputs.dummy_outdirs
            variational_orbital_files = {}

        if not converged:
            utils.warn('The screening parameters have been calculated but are not necessarily self-consistent. '
                       'You may want to increase `alpha_numsteps` to obtain a more accurate result.')

        self.outputs = CalculateScreeningViaDSCFOutput(
            n_electron_restart_dir=iteration_wf.outputs.n_electron_restart_dir)

        self.status = Status.COMPLETED

        return


class DeltaSCFIterationOutputs(IOModel):
    converged: bool
    n_electron_restart_dir: File
    dummy_outdirs: Dict[Tuple[int, int], File | None]
    model_config = ConfigDict(arbitrary_types_allowed=True)


class DeltaSCFIterationWorkflow(Workflow[DeltaSCFIterationOutputs]):

    output_model = DeltaSCFIterationOutputs

    def __init__(self, *args, variational_orbital_files: Dict[str, File],
                 previous_n_electron_calculation=calculators.KoopmansCPCalculator,
                 dummy_outdirs: Dict[Tuple[int, int], File | None], i_sc: int,
                 alpha_indep_calcs: List[calculators.KoopmansCPCalculator], precomputed_descriptors: List[File] | None, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._variational_orbital_files = variational_orbital_files
        self._previous_n_electron_calculation = previous_n_electron_calculation
        self._precomputed_descriptors = precomputed_descriptors
        self._dummy_outdirs = dummy_outdirs
        self._i_sc = i_sc
        self._alpha_indep_calcs = alpha_indep_calcs

        # Set a more instructive name
        self.name = 'Iteration_' + str(i_sc)

    def _run(self) -> None:
        # Do a KI/KIPZ calculation with the updated alpha values
        restart_from_wannier_pwscf = 'evc_occupied1.dat' in self._variational_orbital_files
        if self.parameters.task in ['singlepoint', 'trajectory'] and self.ml.descriptor == 'orbital_density':
            print_real_space_density = True
        else:
            print_real_space_density = False
        assert self.bands is not None
        trial_calc = internal_new_kcp_calculator(self, calc_presets=self.parameters.functional.replace('pkipz', 'ki'),
                                                 print_real_space_density=print_real_space_density,
                                                 alphas=self.bands.alphas,
                                                 restart_from_wannier_pwscf=restart_from_wannier_pwscf)

        # Link the temporary files from the previous calculation
        previous_calc = self._previous_n_electron_calculation
        assert isinstance(previous_calc, calculators.KoopmansCPCalculator)
        trial_calc.link(previous_calc.write_directory, trial_calc.read_directory, recursive_symlink=True)

        for filename, src_file in self._variational_orbital_files.items():
            trial_calc.link(src_file, trial_calc.read_directory / 'K00001' / filename, symlink=True, overwrite=True)

        if self.parameters.functional == 'kipz' and not all(self.atoms.pbc):
            # For the first KIPZ trial calculation, do the innerloop
            trial_calc.parameters.do_innerloop = True

        # Run the calculation
        if self.parameters.fix_spin_contamination:
            status = spin_symmetrize(self, trial_calc)
            if status != Status.COMPLETED:
                return
        else:
            status = self.run_steps(trial_calc)
            if status != Status.COMPLETED:
                return
        alpha_dep_calcs = [trial_calc]

        # Update the bands' self-Hartree and energies (assuming spin-symmetry)
        self.bands.self_hartrees = trial_calc.results['orbital_data']['self-Hartree']

        # Group the bands
        self.bands.assign_groups(allow_reassignment=True)

        skipped_orbitals = []
        first_band_of_each_channel = [self.bands.get(spin=spin)[0] for spin in range(2)]

        # Calculate the power spectrum if required
        if self.ml.descriptor == 'orbital_density' and (self.ml.train or self.ml.predict or self.ml.test) \
                and self.ml.estimator != 'mean':
            if self._precomputed_descriptors is None:
                if self.ml.descriptor == 'orbital_density':
                    psfit_workflow = PowerSpectrumDecompositionWorkflow.fromparent(
                        self, calc_that_produced_orbital_densities=trial_calc)
                    psfit_workflow.run()
                    if psfit_workflow.status != Status.COMPLETED:
                        return

                    descriptors = psfit_workflow.outputs.descriptors
            else:
                descriptors = self._precomputed_descriptors
            for band, power_spectrum in zip(self.bands.to_solve, descriptors):
                band.power_spectrum = power_spectrum

        # Loop over removing/adding an electron from/to each orbital
        assert self.bands is not None
        for band in self.bands:
            # For a KI calculation with only filled bands, we don't have any further calculations to
            # do, so in this case don't print any headings
            print_headings = self.parameters.functional != 'ki' or not band.filled

            # Working out what to print for the orbital heading (grouping skipped bands together)
            if band in self.bands.to_solve or band == self.bands.get(spin=band.spin)[-1]:
                if band not in self.bands.to_solve and (self.parameters.spin_polarized or band.spin == 0):
                    skipped_orbitals.append(band.index)
                if len(skipped_orbitals) > 0:
                    skipped_orbitals = []
                if band not in self.bands.to_solve:
                    continue
            elif not self.parameters.spin_polarized and band.spin == 1:
                # In this case, skip over the bands entirely and don't include it in the printout about which
                # bands we've skipped
                continue
            else:
                # Skip the bands which can copy the screening parameter from another
                # calculation in the same orbital group
                skipped_orbitals.append(band.index)
                continue

            # Use the ML model to predict the screening parameters
            if self.ml.predict or self.ml.test:
                assert self.ml_model is not None
                alpha_pred = self.ml_model.predict(band)
            else:
                alpha_pred = None

            if self.ml.predict:
                alpha = None
                error = None
            else:
                # Calculate the screening parameters ab initio
                assert isinstance(band.index, int)
                dummy_outdir = self._dummy_outdirs.get((band.index, band.spin), None)
                subwf = OrbitalDeltaSCFWorkflow.fromparent(
                    self, band=band, trial_calc=trial_calc, dummy_outdir=dummy_outdir, i_sc=self._i_sc, alpha_indep_calcs=self._alpha_indep_calcs)
                subwf.run()
                if subwf.status != Status.COMPLETED:
                    return
                alpha = subwf.outputs.alpha
                error = subwf.outputs.error
                self._dummy_outdirs[(band.index, band.spin)] = subwf.outputs.dummy_outdir

            for b in self.bands:
                if b == band or (b.group is not None and b.group == band.group):
                    if alpha:
                        b.alpha = alpha
                    if error:
                        b.error = error
                    if alpha_pred:
                        b.predicted_alpha = alpha_pred

            # add alpha to training data
            if self.ml.train:
                assert self.ml_model is not None
                self.ml_model.add_training_data([band])
                # if the user wants to train on the fly, train the model after the calculation of each orbital
                if self.ml.train_on_the_fly:
                    assert self.ml_model is not None
                    self.ml_model.train()

        if not self.steps_are_running():
            # We only want to print the history if we are not waiting for the results of a calculation
            print_alpha_history(self)

        assert isinstance(self.ml.predict, bool)
        converged = self.ml.predict or all([abs(b.error) < 1e-3 for b in self.bands])

        if self.ml.train:
            # if the user doesn't want to train on the fly, train the model at the end of each snapshot
            if not self.ml.train_on_the_fly:
                assert self.ml_model is not None
                self.ml_model.train()

        self.outputs = DeltaSCFIterationOutputs(converged=converged,
                                                n_electron_restart_dir=trial_calc.write_directory,
                                                dummy_outdirs=self._dummy_outdirs)

        self.status = Status.COMPLETED


class OrbitalDeltaSCFOutputs(IOModel):
    alpha: float
    error: float
    dummy_outdir: File | None
    model_config = ConfigDict(arbitrary_types_allowed=True)


class OrbitalDeltaSCFWorkflow(Workflow[OrbitalDeltaSCFOutputs]):

    output_model = OrbitalDeltaSCFOutputs

    def __init__(self, band: Band, trial_calc: calculators.KoopmansCPCalculator,
                 dummy_outdir: File | None, i_sc: int,
                 alpha_indep_calcs: List[calculators.KoopmansCPCalculator],
                 **kwargs):
        super().__init__(**kwargs)
        self.band = band
        self._trial_calc = trial_calc
        self._dummy_outdir = dummy_outdir
        self._i_sc = i_sc
        self._alpha_indep_calcs = alpha_indep_calcs

        # Set a more instructive name
        self.name = 'Orbital ' + str(self.band.index)
        if self.parameters.spin_polarized:
            self.name += ' Spin ' + str(self.band.spin + 1)

    def _run(self) -> None:

        assert self.bands is not None

        alpha_dep_calcs = [self._trial_calc]

        # Don't repeat if this particular alpha_i was converged
        if hasattr(self.band, 'error') and abs(self.band.error) < self.parameters.alpha_conv_thr:
            assert self.band.alpha is not None
            self.outputs = self.output_model(alpha=self.band.alpha, error=self.band.error,
                                             dummy_outdir=self._dummy_outdir)
            self.status = Status.COMPLETED
            return

        # When we write/update the alpharef files in the work directory
        # make sure to include the fixed band alpha in file_alpharef.txt
        # rather than file_alpharef_empty.txt
        if self.band.filled:
            index_empty_to_save = None
        else:
            index_empty_to_save = self.band.index - self.bands.num(filled=True, spin=self.band.spin)
            if self.parameters.spin_polarized and self.band.spin == 1:
                index_empty_to_save += self.bands.num(filled=False, spin=0)

        # Perform the fixed-band-dependent calculations
        if self.parameters.functional in ['ki', 'pkipz']:
            if self.band.filled:
                calc_types = ['dft_n-1']
            else:
                if self._i_sc == 1:
                    calc_types = ['dft_n+1_dummy', 'pz_print', 'dft_n+1']
                else:
                    assert self._dummy_outdir is not None
                    calc_types = ['pz_print', 'dft_n+1']
        else:
            if self.band.filled:
                calc_types = ['kipz_n-1']
            else:
                if self._i_sc == 1:
                    calc_types = ['dft_n+1_dummy', 'kipz_print', 'kipz_n+1']
                else:
                    assert self._dummy_outdir is not None
                    calc_types = ['kipz_print', 'kipz_n+1']

        dummy_outdir = self._dummy_outdir
        print_calc = None

        for calc_type in calc_types:
            if self.parameters.functional in ['ki', 'pkipz']:
                # The calculations whose results change with alpha are...
                #  - the KI calculations
                #  - DFT calculations on empty variational orbitals
                # We don't need to redo any of the others
                if not self._trial_calc.has_empty_states() or self.band.filled:
                    if self._i_sc > 1 and 'ki' not in calc_type:
                        self.print('No further calculations are required to calculate this screening parameter')
                        continue

            if 'print' in calc_type:
                # Note that the 'print' calculations for empty bands do not
                # in fact involve the fixing of that band (and thus for the
                # 'fixed' band the corresponding alpha should be in
                # file_alpharef_empty.txt)
                alphas = self.bands.alphas
                filling = self.bands.filling
            elif not self.band.filled:
                # In the case of empty orbitals, we gain an extra orbital in
                # the spin-up channel, so we explicitly construct both spin
                # channels for "alphas" and "filling"
                alphas = self.bands.alphas
                alphas[self.band.spin].append(alphas[self.band.spin][-1])
                filling = self.bands.filling
                assert self.band.index is not None
                filling[self.band.spin][self.band.index - 1] = True
                filling[self.band.spin].append(False)
            else:
                alphas = self.bands.alphas
                filling = self.bands.filling

            # Work out the index of the band that is fixed (noting that we will be throwing away all empty
            # bands)
            fixed_band = min(self.band.index, self.bands.num(filled=True, spin=self.band.spin) + 1)
            if self.parameters.spin_polarized and self.band.spin == 1:
                fixed_band += self.bands.num(filled=True, spin=0)

            # Set up calculator
            calc = internal_new_kcp_calculator(self, calc_type, alphas=alphas, filling=filling, fixed_band=fixed_band,
                                               index_empty_to_save=index_empty_to_save,
                                               add_to_spin_up=(self.band.spin == 0))

            if calc.parameters.ndr == self._trial_calc.parameters.ndw:
                calc.link(self._trial_calc.write_directory, calc.read_directory,
                          recursive_symlink=True)

            if calc_type in ['dft_n+1', 'kipz_n+1']:
                assert dummy_outdir is not None
                calc.link(dummy_outdir, calc.parameters.outdir, recursive_symlink=True)
                # Copying of evcfixed_empty.dat to evc_occupied.dat
                assert print_calc is not None
                for ispin in range(1, 3):
                    calc.link(print_calc.write_directory / f'K00001/evcfixed_empty{ispin}.dat',
                              calc.read_directory / f'K00001/evc_occupied{ispin}.dat', symlink=True, overwrite=True)

            # Run kcp.x
            if calc.parameters.nelup < calc.parameters.neldw:
                subwf = KoopmansCPWithSpinSwapWorkflow.fromparent(self, calc=calc)
                subwf.run()
                if subwf.status != Status.COMPLETED:
                    return
                if 'dummy' in calc_type:
                    dummy_outdir = subwf.outputs.outdir
            else:
                self.run_steps(calc)
                if 'dummy' in calc_type:
                    dummy_outdir = File(calc, calc.parameters.outdir)

            # Store the band that we've perturbed as calc.fixed_band. Note that we can't use
            # calc.parameters.fixed_band to keep track of which band we held fixed, because for empty
            # orbitals, calc.parameters.fixed_band is always set to the LUMO but in reality we're fixing
            # the band corresponding # to index_empty_to_save from an earlier calculation
            calc.fixed_band = self.band

            # Store the result
            # We store the results in one of two lists: alpha_indep_calcs and
            # alpha_dep_calcs. The latter is overwritten at each new self-
            # consistency loop.
            if 'ki' in calc_type and 'print' not in calc_type:
                alpha_dep_calcs.append(calc)
            elif 'dft' in calc_type and 'dummy' not in calc_type:
                if self.parameters.functional in ['ki', 'pkipz']:
                    # For KI, the results of the DFT calculations are typically independent of alpha so we
                    # store these in a list that is never overwritten

                    # The exception to this are KI calculations on empty states. When we update alpha, the
                    # empty manifold changes, which in turn affects the lambda values
                    if self._trial_calc.has_empty_states() and not self.band.filled:
                        alpha_dep_calcs.append(calc)
                    else:
                        self._alpha_indep_calcs.append(calc)
                else:
                    # For KIPZ, the DFT calculations are dependent on alpha via
                    # the definition of the variational orbitals. We only want to
                    # store the calculations that used the most recent value of alpha

                    alpha_dep_calcs.append(calc)

            # Storing the calculators to allow for the copying of evcfixed_empty.dat to evc_occupied.dat
            if calc_type in ['pz_print', 'kipz_print']:
                print_calc = calc

        # Calculate an updated alpha and a measure of the error
        # E(N) - E_i(N - 1) - lambda^alpha_ii(1)     (filled)
        # E_i(N + 1) - E(N) - lambda^alpha_ii(0)     (empty)

        calcs = [c for c in alpha_dep_calcs + self._alpha_indep_calcs if c.fixed_band == self.band]

        alpha, error = self.calculate_alpha_from_list_of_calcs(
            calcs, self._trial_calc, self.band, filled=self.band.filled)

        # Mixing
        alpha = self.parameters.alpha_mixing * alpha + (1 - self.parameters.alpha_mixing) * self.band.alpha

        warning_message = 'The computed screening parameter is {0}. Proceed with caution.'
        failure_message = 'The computed screening parameter is significantly {0}. This should not ' \
            'happen. Decrease `alpha_mixing` and/or change `alpha_guess`.'

        if alpha < -0.1:
            raise ValueError(failure_message.format('less than 0'))
        elif alpha < 0:
            utils.warn(warning_message.format('less than 0'))
        elif alpha > 1.1:
            raise ValueError(failure_message.format('greater than 1'))
        elif alpha > 1:
            utils.warn(warning_message.format('greater than 1'))

        self.outputs = self.output_model(alpha=alpha, error=error, dummy_outdir=dummy_outdir)

        self.status = Status.COMPLETED

    def calculate_alpha_from_list_of_calcs(self,
                                           calcs: List[calculators.KoopmansCPCalculator],
                                           trial_calc: calculators.KoopmansCPCalculator,
                                           band: Band,
                                           filled: bool = True) -> Tuple[float, float]:
        '''

        Calculates alpha via equation 10 of Nguyen et. al (2018) 10.1103/PhysRevX.8.021051
        If the band is filled, use s = 1; if the band is empty, use s = 0

        Arguments:
            calcs          -- a list of selected calculations from which to calculate alpha
            trial_calc     -- the N-electron Koopmans calculation
            filled         -- True if the orbital for which we're calculating alpha is filled

        '''

        # Extract the energy difference Delta E
        if self.parameters.functional == 'kipz':
            if filled:
                # KIPZ N-1
                [kipz_m1_calc] = [c for c in calcs if c.parameters.which_orbdep == 'nkipz'
                                  and c.parameters.do_orbdep and c.parameters.f_cutoff < 0.0001]
                kipz_m1 = kipz_m1_calc.results
                charge = 1 - kipz_m1_calc.parameters.f_cutoff

                dE = trial_calc.results['energy'] - kipz_m1['energy']
                mp1 = kipz_m1['mp1_energy']
                mp2 = kipz_m1['mp2_energy']

            else:
                # KIPZ N+1
                [kipz_p1_calc] = [c for c in calcs if c.parameters.which_orbdep == 'nkipz'
                                  and c.parameters.do_orbdep and c.parameters.f_cutoff == 1.0]
                kipz_p1 = kipz_p1_calc.results
                charge = - kipz_p1_calc.parameters.f_cutoff

                dE = kipz_p1['energy'] - trial_calc.results['energy']
                mp1 = kipz_p1['mp1_energy']
                mp2 = kipz_p1['mp2_energy']

        else:
            # self.functional in ['ki', 'pkipz']
            if filled:
                # DFT N-1
                [dft_m1_calc] = [c for c in calcs if not c.parameters.do_orbdep
                                 and c.parameters.restart_mode == 'restart' and c.parameters.f_cutoff < 0.0001]
                dft_m1 = dft_m1_calc.results
                charge = 1 - dft_m1_calc.parameters.f_cutoff

                dE = trial_calc.results['energy'] - dft_m1['energy']
                mp1 = dft_m1['mp1_energy']
                mp2 = dft_m1['mp2_energy']

            else:
                # DFT N+1
                [dft_p1_calc] = [c for c in calcs if not c.parameters.do_orbdep
                                 and c.parameters.restart_mode == 'restart' and c.parameters.f_cutoff == 1.0]
                dft_p1 = dft_p1_calc.results
                charge = - dft_p1_calc.parameters.f_cutoff

                dE = dft_p1['energy'] - trial_calc.results['energy']
                mp1 = dft_p1['mp1_energy']
                mp2 = dft_p1['mp2_energy']

        # Extract lambda from the base calculator
        assert band.index is not None
        iband = band.index - 1  # converting from 1-indexing to 0-indexing
        lambda_a = trial_calc.results['lambda'][band.spin][iband, iband].real
        lambda_0 = trial_calc.results['bare lambda'][band.spin][iband, iband].real

        # Obtaining alpha
        if (trial_calc.parameters.odd_nkscalfact and filled) \
                or (trial_calc.parameters.odd_nkscalfact_empty and not filled):
            alpha_guess = trial_calc.alphas[band.spin][iband]
        else:
            alpha_guess = trial_calc.parameters.nkscalfact

        # Checking Makov-Payne correction energies and applying them (if needed)
        if self.parameters.mp_correction:
            if mp1 is None:
                raise ValueError('Could not find 1st order Makov-Payne energy')
            if mp2 is None:
                # utils.warn('Could not find 2nd order Makov-Payne energy; applying first order only')
                mp_energy = mp1
            else:
                mp_energy = mp1 + mp2

            dE -= np.sign(charge) * mp_energy / self.parameters.eps_inf

        alpha = alpha_guess * (dE - lambda_0) / (lambda_a - lambda_0)

        # The error is lambda^alpha(1) - lambda^alpha_i(1)
        error = dE - lambda_a

        return alpha, error


def internal_new_kcp_calculator(workflow,
                                calc_presets: str = 'dft_init',
                                alphas: Optional[List[List[float]]] = None,
                                filling: Optional[List[List[bool]]] = None,
                                add_to_spin_up: bool = True,
                                **kwargs) -> calculators.KoopmansCPCalculator:
    """

    Generates a new KCP calculator based on the self.calculator_parameters["kcp"]
    parameters, modifying the appropriate settings to match the
    chosen calc_presets, and altering any Quantum Espresso keywords
    specified as kwargs

    Arguments:

        calc_presets
            The set of preset values to use; must be one of the following strings:

            Initialization
            'dft_init'            DFT calculation from scratch
            'pz_init'             PZ calculation starting from DFT restart
            'pz_innerloop_init'   PZ calculation starting from DFT restart (innerloop only)
            'dft_dummy'           DFT dummy calculation that generate store files
                                  for a periodic calculation restarting from Wannier functions

            Trial calculations
            'ki'     KI calculation with N electrons and empty bands if specified
            'kipz'   As above, but for KIPZ

            For calculating alpha_i for filled orbitals.
            'dft_n-1'       DFT calculation with N-1 electrons via fixed_state
            'kipz_n-1'      KIPZ calculation with N-1 electrons via fixed_state

            For calculating alpha_i for empty orbitals
            'pz_print'         PZ calculation that generates evcfixed_empty.dat file
            'kipz_print'       KIPZ calculation that generates evcfixed_empty.dat file
            'dft_n+1_dummy'    DFT dummy calculation that generates store files of
                               the correct dimensions
            'dft_n+1'          DFT calculation with N+1 electrons
            'kipz_n+1'         KIPZ calculation with N+1 electrons

            Note that when calculating alpha_i, all empty states (bar orbital_i if it is empty) are removed
            and the convergence criteria are loosened

            Final calculation
            'ki_final'    Final KI calculation with N electrons and empty bands if specified
            'kipz_final'  As above, but for KIPZ
            'pkipz_final' A KIPZ calculation that leaves the manifold unchanged (for the
                          purposes of performing KIPZ on top of the KI manifold)

        alphas
            an array of screening parameters

       **kwargs accepts any Quantum Espresso keywords as an argument, and will
                apply these options to the returned ASE calculator

    Returns: a new KCP calculator object

    """

    # By default, use the last row in the alpha table for the screening parameters
    if alphas is None:
        alphas = workflow.bands.alphas

    # Generate a new kcp calculator copied from the master calculator
    calc: calculators.KoopmansCPCalculator = workflow.new_calculator('kcp', alphas=alphas, filling=filling)

    # Set up read/write indexes
    if calc_presets in ['dft_init', 'dft_dummy']:
        ndr = 50
        ndw = 50
    elif calc_presets in ['pz_init', 'pz_innerloop_init']:
        ndr = 50
        ndw = 51
    elif calc_presets in ['ki', 'kipz']:
        ndr = 51
        ndw = 60
    elif calc_presets in ['dft_n-1', 'kipz_n-1']:
        ndr = 60
        ndw = 63
    elif calc_presets in ['pz_print', 'kipz_print']:
        ndr = 60
        ndw = 64
    elif calc_presets == 'dft_n+1_dummy':
        ndr = 65
        ndw = 65
    elif calc_presets in ['dft_n+1', 'kipz_n+1']:
        ndr = 65
        ndw = 68
    elif calc_presets in ['ki_final', 'kipz_final']:
        if workflow.parameters.calculate_alpha:
            ndr = 60
        else:
            ndr = 51
        ndw = 70
    elif calc_presets in ['pkipz_final']:
        ndr = 70
        ndw = 71
    else:
        raise ValueError('Invalid calc_presets "{}"'.format(calc_presets))

    # KCP options
    # control
    calc.prefix = calc_presets
    calc.parameters.ndw = ndw
    calc.parameters.ndr = ndr
    if calc.prefix in ['dft_init', 'dft_dummy', 'dft_n+1_dummy']:
        calc.parameters.restart_mode = 'from_scratch'
    else:
        calc.parameters.restart_mode = 'restart'

    # system
    if 'pz' in calc.prefix or 'ki' in calc.prefix:
        calc.parameters.do_orbdep = True
        calc.parameters.do_bare_eigs = True
    else:
        calc.parameters.do_orbdep = False
    if calc.prefix in ['dft_init', 'pz_init', 'pz_innerloop_init', 'dft_dummy',
                       'ki', 'kipz', 'pz_print', 'kipz_print', 'dft_n+1_dummy',
                       'ki_final', 'kipz_final', 'pkipz_final']:
        calc.parameters.fixed_state = False
    else:
        calc.parameters.fixed_state = True
        if '-1' in calc.prefix:
            calc.parameters.f_cutoff = 1e-5
        else:
            calc.parameters.f_cutoff = 1.0
    if 'n+1' in calc.prefix:
        calc.parameters.nelec += 1
        if add_to_spin_up:
            calc.parameters.nelup += 1
        else:
            calc.parameters.neldw += 1
        if 'dummy' not in calc.prefix:
            calc.parameters.restart_from_wannier_pwscf = True

    # electrons
    # For all calculations calculating alpha, remove the empty states and
    # increase the energy thresholds
    if not any([s in calc.prefix for s in ['init', 'print', 'final']]) and \
       calc.prefix not in ['ki', 'kipz', 'dft_dummy']:
        calc.parameters.nbnd = None
        calc.parameters.conv_thr *= 100
        calc.parameters.esic_conv_thr *= 100

    # For the dft_dummy calculation, we don't need empty states because these will be overwritten by the w90
    # wavefunctions
    if calc.prefix == 'dft_dummy':
        calc.parameters.nbnd = None

    if all(workflow.atoms.pbc) and not any([s == calc.prefix for s in ['dft_init', 'dft_n-1', 'dft_n+1',
                                                                       'kipz', 'kipz_n-1', 'kipz_n+1']]):
        calc.parameters.do_outerloop = False
        calc.parameters.do_innerloop = False
    elif any([s in calc.prefix for s in ['frozen', 'dummy', 'print', 'innerloop']]) or calc.prefix == 'pkipz_final':
        calc.parameters.do_outerloop = False
        if calc.has_empty_states():
            calc.parameters.do_outerloop_empty = False
    elif calc.prefix in ['ki', 'ki_final']:
        calc.parameters.do_outerloop = False
        if calc.has_empty_states():
            if workflow.parameters.init_empty_orbitals == 'pz':
                calc.parameters.do_outerloop_empty = True
            else:
                calc.parameters.do_outerloop_empty = False
    else:
        calc.parameters.do_outerloop = True
        if calc.has_empty_states():
            calc.parameters.do_outerloop_empty = True

    if calc.parameters.maxiter is None and calc.parameters.do_outerloop:
        calc.parameters.maxiter = 300
    if calc.parameters.empty_states_maxstep is None and calc.parameters.do_outerloop_empty:
        calc.parameters.empty_states_maxstep = 300

    # No empty states minimization in the solids workflow for the moment
    if all(workflow.atoms.pbc) and calc.has_empty_states():
        calc.parameters.do_outerloop_empty = False
        calc.parameters.do_innerloop_empty = False

    # nksic
    if calc.parameters.do_orbdep:
        calc.parameters.odd_nkscalfact = True
        calc.parameters.odd_nkscalfact_empty = True
    calc.parameters.do_innerloop_cg = True
    if calc.prefix[:2] == 'pz' and 'print' not in calc.prefix:
        calc.parameters.do_innerloop = True
    else:
        calc.parameters.do_innerloop = False
    if calc.has_empty_states():
        calc.parameters.do_innerloop_empty = False
    if 'kipz' in calc.prefix:
        calc.parameters.which_orbdep = 'nkipz'
    elif 'pz' in calc.prefix:
        calc.parameters.which_orbdep = 'pz'
    elif 'ki' in calc.prefix:
        calc.parameters.which_orbdep = 'nki'
    if 'print' in calc.prefix:
        calc.parameters.print_wfc_anion = True

    if workflow.parameters.mt_correction:
        calc.parameters.which_compensation = 'tcc'
    else:
        calc.parameters.which_compensation = 'none'

    # If we are using frozen orbitals, we override the above logic and freeze the variational orbitals
    # post-initialization
    if workflow.parameters.frozen_orbitals and 'init' not in calc.prefix and not any([s == calc.prefix for s in
                                                                                      ['dft_n-1', 'dft_n+1', 'kipz_n-1',
                                                                                       'kipz_n+1']]):
        calc.parameters.do_outerloop = False
        calc.parameters.do_innerloop = False
        if calc.has_empty_states():
            calc.parameters.do_outerloop_empty = False
            calc.parameters.do_innerloop_empty = False

    # Handle any keywords provided by kwargs
    # Note that since this is performed after the above logic this can (deliberately
    # or accidentally) overwrite the above settings
    calc.parameters.update(**kwargs)

    # Sanity checking
    if calc.parameters.print_wfc_anion and calc.parameters.index_empty_to_save is None:
        raise ValueError('`print_wfc_anion` is set to `True` but you have not selected '
                         'an `index_empty_to_save`. Provide this as an argument to `new_cp_calculator`')

    # don't print QC in some cases
    if 'dummy' in calc.prefix or calc.prefix[-2:] == '+1':
        calc.skip_qc = True

    return calc


class InitializationWorkflow(Workflow[KoopmansDSCFOutputs]):

    output_model = KoopmansDSCFOutputs

    def _run(self) -> None:
        wannier_hamiltonian_files: Dict[BlockID, File] | None = None

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs'] or \
                (all(self.atoms.pbc) and self.parameters.init_orbitals == 'kohn-sham'):
            # Wannier functions using pw.x, wannier90.x and pw2wannier90.x (pw.x only for Kohn-Sham states)
            wannier_workflow = WannierizeWorkflow.fromparent(self)
            if wannier_workflow.parameters.calculate_bands:
                wannier_workflow.parameters.calculate_bands = \
                    not self.calculator_parameters['ui'].do_smooth_interpolation

            # Perform the wannierization workflow
            wannier_workflow.run()
            if wannier_workflow.status != Status.COMPLETED:
                return

            # Store the Hamitonian files
            hr_file_keys: List[BlockID]
            if self.parameters.spin_polarized:
                hr_file_ids = [BlockID(filled=True, spin='up'),
                               BlockID(filled=False, spin='up'),
                               BlockID(filled=True, spin='down'),
                               BlockID(filled=False, spin='down')]
            else:
                hr_file_ids = [BlockID(filled=True), BlockID(filled=False)]

            wannier_hamiltonian_files = {}
            for b_id in hr_file_ids:
                hr_file = wannier_workflow.outputs.hr_files[b_id]
                assert hr_file is not None
                wannier_hamiltonian_files[b_id] = hr_file

            # Convert the files over from w90 format to kcp format
            nscf_calc = wannier_workflow.outputs.nscf_calculation
            nscf_outdir = File(nscf_calc, nscf_calc.parameters.outdir)
            fold_workflow = FoldToSupercellWorkflow.fromparent(self, nscf_outdir=nscf_outdir,
                                                               hr_files=wannier_workflow.outputs.hr_files,
                                                               wannier90_calculations=wannier_workflow.outputs.wannier90_calculations,
                                                               wannier90_pp_calculations=wannier_workflow.outputs.preprocessing_calculations)
            fold_workflow.run()
            if fold_workflow.status != Status.COMPLETED:
                return

            # Convert self.atoms to the supercell
            self.primitive_to_supercell()

            # We need a dummy calc before the real dft_init in order
            # to copy the previously calculated Wannier functions
            dummy_calc = internal_new_kcp_calculator(self, 'dft_dummy')
            status = self.run_steps(dummy_calc)
            if status != Status.COMPLETED:
                return

            # DFT restarting from Wannier functions (after copying the Wannier functions)
            calc = internal_new_kcp_calculator(self, 'dft_init', restart_mode='restart',
                                               restart_from_wannier_pwscf=True, do_outerloop=True)
            calc.link(File(dummy_calc, dummy_calc.parameters.outdir), calc.parameters.outdir)
            for filename, f in fold_workflow.outputs.kcp_files.items():
                calc.link(f, calc.read_directory / 'K00001' / filename, symlink=True)
            status = self.run_steps(calc)
            if status != Status.COMPLETED:
                return

            # Check the consistency between the PW and CP band gaps
            pw_calc = [c for c in self.calculations if isinstance(
                c, calculators.PWCalculator) and c.parameters.calculation == 'nscf'][-1]
            pw_gap = pw_calc.results['lumo_energy'] - pw_calc.results['homo_energy']
            cp_gap = calc.results['lumo_energy'] - calc.results['homo_energy']
            if abs(pw_gap - cp_gap) > 2e-2 * pw_gap:
                raise ValueError(f'PW and CP band gaps are not consistent: {pw_gap} {cp_gap}')

            # The CP restarting from Wannier functions must be already converged
            Eini = calc.results['convergence']['filled'][0]['Etot']
            Efin = calc.results['energy']
            if abs(Efin - Eini) > 1e-6 * abs(Efin):
                raise ValueError(f'Too much difference between the initial and final CP energies: {Eini} {Efin}')

        elif self.parameters.functional in ['ki', 'pkipz']:
            calc_dft = internal_new_kcp_calculator(self, 'dft_init')
            if self.parameters.fix_spin_contamination:
                status = spin_symmetrize(self, calc_dft)
                if status != Status.COMPLETED:
                    return
            else:
                status = self.run_steps(calc_dft)
                if status != Status.COMPLETED:
                    return

            if self.parameters.init_orbitals == 'kohn-sham':
                calc = calc_dft
            elif self.parameters.init_orbitals == 'pz':
                if self.calculations[-1].parameters.nelec == 2:
                    # If we only have two electrons, then the filled manifold is trivially invariant under unitary
                    # transformations. Furthermore, the PZ functional is invariant w.r.t. unitary rotations of the
                    # empty states. Thus in this instance we can skip the initialization of the manifold entirely
                    self.print('Skipping the optimisation of the variational orbitals since they are invariant under '
                               'unitary transformations')
                else:
                    assert self.bands is not None
                    calc = internal_new_kcp_calculator(self, 'pz_innerloop_init', alphas=self.bands.alphas)

                    # Use the KS eigenfunctions as initial guesses for the variational orbitals
                    calc.link(File(calc_dft, calc_dft.parameters.outdir),
                              calc.parameters.outdir, recursive_symlink=True)
                    ndw_dir = calc_dft.write_directory / 'K00001'
                    ndr_dir = calc.read_directory / 'K00001'
                    calc.link(ndw_dir / 'evc1.dat', ndr_dir / 'evc01.dat', symlink=True, overwrite=True)
                    calc.link(ndw_dir / 'evc2.dat', ndr_dir / 'evc02.dat', symlink=True, overwrite=True)
                    status = self.run_steps(calc)
                    if status != Status.COMPLETED:
                        return
            else:
                raise ValueError('Should not arrive here')

        elif self.parameters.functional == 'kipz':
            # DFT from scratch
            calc_dft = internal_new_kcp_calculator(self, 'dft_init')
            if self.parameters.fix_spin_contamination:
                status = spin_symmetrize(self, calc_dft)
                if status != Status.COMPLETED:
                    return
            else:
                status = self.run_steps(calc_dft)
                if status != Status.COMPLETED:
                    return

            if self.parameters.init_orbitals == 'kohn-sham':
                # Initialize the density with DFT and use the KS eigenfunctions as guesses for the variational orbitals
                raise NotImplementedError()
                # self._overwrite_canonical_with_variational_orbitals(calc)
                # self._copy_most_recent_calc_to_ndw(ndw_final)
            elif self.parameters.init_orbitals == 'pz':
                # PZ from DFT (generating PZ density and PZ orbitals)
                calc = internal_new_kcp_calculator(self, 'pz_init')
                calc.link(calc_dft.write_directory, calc.read_directory, recursive_symlink=True)
                status = self.run_steps(calc)
                if status != Status.COMPLETED:
                    return
            else:
                raise ValueError('Should not arrive here')

        else:
            raise ValueError("Should not arrive here; there must be an inconsistency between the above code and \
                             `workflow.valid_settings`")

        # Compiling a dictionary of the variational orbital files
        variational_orbitals = {}
        if self.parameters.init_orbitals == 'kohn-sham':
            prefix = 'evc'
        else:
            prefix = 'evc0'
        ndw_dir = calc.write_directory / 'K00001'
        ndr_dir = calc.read_directory / 'K00001'
        for ispin in range(1, 3):
            if self.parameters.init_orbitals in ['mlwfs', 'projwfs'] or \
                    (all(self.atoms.pbc) and self.parameters.init_orbitals == 'kohn-sham'):
                variational_orbitals[f'evc_occupied{ispin}.dat'] = ndr_dir / f'evc_occupied{ispin}.dat'
            else:
                variational_orbitals[f'evc0{ispin}.dat'] = ndw_dir / f'{prefix}{ispin}.dat'
            if calc.has_empty_states(ispin - 1):
                variational_orbitals[f'evc0_empty{ispin}.dat'] = ndw_dir / f'{prefix}_empty{ispin}.dat'

        # If fixing spin contamination, copy the spin-up variational orbitals to the spin-down channel
        if self.parameters.fix_spin_contamination:
            for spin_up_file in ['evc01.dat', 'evc0_empty1.dat', 'evc_occupied1.dat']:
                if spin_up_file in variational_orbitals:
                    variational_orbitals[spin_up_file.replace('1', '2')] = variational_orbitals[spin_up_file]

        self.outputs = self.output_model(variational_orbital_files=variational_orbitals, final_calc=calc,
                                         wannier_hamiltonian_files=wannier_hamiltonian_files)

        self.status = Status.COMPLETED
        return


def print_alpha_history(wf: Workflow):
    # Printing out a progress summary
    assert wf.bands is not None
    if not wf.ml.predict:
        wf.print(f'\n**α**')
        if wf.parameters.spin_polarized:
            wf.print('\n**spin up**')
            wf.print(wf.bands.alpha_history(spin=0).to_markdown(), wrap=False)
            wf.print('\n**spin down**')
            wf.print(wf.bands.alpha_history(spin=1).to_markdown(), wrap=False)
        else:
            wf.print(wf.bands.alpha_history().to_markdown(), wrap=False)

    if None not in [b.predicted_alpha for b in wf.bands]:
        wf.print(f'\n**predicted α**')
        if wf.parameters.spin_polarized:
            wf.print('\n**spin up**')
            wf.print(wf.bands.predicted_alpha_history(spin=0).to_markdown(), wrap=False)
            wf.print('\n**spin down**')
            wf.print(wf.bands.predicted_alpha_history(spin=1).to_markdown(), wrap=False)
        else:
            wf.print(wf.bands.predicted_alpha_history().to_markdown(), wrap=False)

    if not wf.bands.error_history().empty:
        wf.print(f'\n**ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**')
        if wf.parameters.spin_polarized:
            wf.print('\n**spin up**')
            wf.print(wf.bands.error_history(spin=0).to_markdown(), wrap=False)
            wf.print('\n**spin down**')
            wf.print(wf.bands.error_history(spin=1).to_markdown(), wrap=False)
        else:
            wf.print(wf.bands.error_history().to_markdown(), wrap=False)
    wf.print('')
