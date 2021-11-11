"""

Workflow module for koopmans, containing the workflow for performing KI and KIPZ calculations

Written by Edward Linscott Jan 2020
Split off from workflow.py Oct 2020

"""

import os
import copy
import numpy as np
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import pandas as pd
from ase import Atoms
from ase.dft.dos import DOS
from ase.spectrum.band_structure import BandStructure
from koopmans import utils
from koopmans.settings import KoopmansCPSettingsDict
from koopmans.bands import Band, Bands
from koopmans import calculators
from ._generic import Workflow


class KoopmansDSCFWorkflow(Workflow):

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

        if 'kcp' not in self.master_calc_params:
            raise ValueError(
                'Performing a KC calculation requires a "kcp" block in the .json input file')

        # If periodic, convert the kcp calculation into a Γ-only supercell calculation
        kcp_params = self.master_calc_params['kcp']
        if self.parameters.periodic:
            if self.parameters.spin_polarised:
                raise NotImplementedError('Yet to implement spin-polarised calculations for periodic systems')

            # Update the KCP settings to correspond to a supercell (leaving self.atoms unchanged for the moment)
            self.convert_kcp_to_supercell()
            nocc = self.master_calc_params['w90_occ'].num_wann
            nemp = self.master_calc_params['w90_emp'].num_wann
            # Note that self.parameters.orbital_groups does not have a spin index, as opposed to self.bands.groups
            if self.parameters.orbital_groups is None:
                self.parameters.orbital_groups = list(range(0, nocc + nemp))
            self.parameters.orbital_groups = [i for _ in range(np.prod(self.kgrid))
                                              for i in self.parameters.orbital_groups[:nocc]] \
                + [i for _ in range(np.prod(self.kgrid))
                   for i in self.parameters.orbital_groups[nocc:]]

            # Check the number of empty states has been correctly configured
            w90_emp_params = self.master_calc_params['w90_emp']
            expected_empty_states_nbnd = w90_emp_params.num_wann * np.prod(self.kgrid)
            if kcp_params.empty_states_nbnd == 0:
                # 0 is the default value
                kcp_params.empty_states_nbnd = expected_empty_states_nbnd
            elif kcp_params.empty_states_nbnd != expected_empty_states_nbnd:
                raise ValueError('kcp empty_states_nbnd and wannier90 num_wann (emp) are inconsistent')

        # Initialise the bands object
        if self.parameters.spin_polarised:
            filling = [[True for _ in range(kcp_params.nelup)] + [False for _ in range(kcp_params.empty_states_nbnd)],
                       [True for _ in range(kcp_params.neldw)] + [False for _ in range(kcp_params.empty_states_nbnd)]]
            groups = self.parameters.orbital_groups
            if len(groups) != 2 or not isinstance(groups[0], list):
                raise ValueError('If spin_polarised = True, orbital_groups should be a list containing two sublists '
                                 '(one per spin channel)')
        else:
            filling = [[True for _ in range(kcp_params.nelec // 2)]
                       + [False for _ in range(kcp_params.empty_states_nbnd)] for _ in range(2)]
            # self.parameters.orbital_groups does not have a spin index
            if self.parameters.orbital_groups is None:
                groups = None
            else:
                groups = [self.parameters.orbital_groups for _ in range(2)]

        # Sanity checking
        if groups is not None:
            for g, f in zip(groups, filling):
                assert len(g) == len(f), 'orbital_groups is the wrong dimension; it should have dimensions (2, num_bands)'

        self.bands = Bands(n_bands=len(filling[0]), n_spin=2, spin_polarised=self.parameters.spin_polarised,
                           filling=filling, groups=groups,
                           self_hartree_tol=self.parameters.orbital_groups_self_hartree_tol)

        if self.parameters.alpha_from_file:
            # Reading alpha values from file
            self.bands.alphas = self.read_alphas_from_file()
        else:
            # Initialising alpha with a guess
            self.bands.alphas = self.parameters.alpha_guess

        # Raise errors if any UI keywords are provided but will be overwritten by the workflow
        for ui_keyword in ['kc_ham_file', 'w90_seedname', 'kpts', 'dft_ham_file', 'dft_smooth_ham_file']:
            for ui_kind in ['occ', 'emp']:
                value = getattr(self.master_calc_params[f'ui_{ui_kind}'], ui_keyword)
                [default_value] = [s.default for s in self.master_calc_params['ui'].settings if s.name == ui_keyword]
                if value != default_value:
                    raise ValueError(f'UI keyword {ui_keyword} has been set in the input file, but this will be '
                                     'automatically set by the Koopmans workflow. Remove this keyword from the input '
                                     'file')

        # Initialise self.init_empty_orbitals if it has not been set
        if self.parameters.init_empty_orbitals == 'same':
            self.parameters.init_empty_orbitals = self.parameters.init_orbitals
        if self.parameters.init_empty_orbitals != self.parameters.init_orbitals:
            raise NotImplementedError(f'The combination init_orbitals = {self.parameters.init_orbitals} '
                                      f'and init_empty_orbitals = {self.parameters.init_empty_orbitals} '
                                      'has not yet been implemented')

        # By default, we will always run the higher-res DFT calculations. Changing this flag allows for calculations
        # where this is unnecessary (such as pKIPZ) to skip this step.
        self.redo_preexisting_smooth_dft_calcs = True

    def convert_kcp_to_supercell(self):
        # Multiply all extensive KCP settings by the appropriate prefactor
        prefactor = np.prod(self.kgrid)
        for attr in ['nelec', 'nelup', 'neldw', 'empty_states_nbnd', 'conv_thr', 'esic_conv_thr']:
            value = getattr(self.master_calc_params['kcp'], attr, None)
            if value is not None:
                setattr(self.master_calc_params['kcp'], attr, prefactor * value)

    def read_alphas_from_file(self, directory: Path = Path()):
        '''
        This routine reads in the contents of file_alpharef.txt and file_alpharef_empty.txt and
        stores the result in self.bands.alphas

        Since utils.read_alpha_file provides a flattened list of alphas so we must convert this
        to a nested list using convert_flat_alphas_for_kcp()
        '''

        flat_alphas = utils.read_alpha_file(directory)
        params = self.master_calc_params['kcp']
        assert isinstance(params, KoopmansCPSettingsDict)
        alphas = calculators.convert_flat_alphas_for_kcp(flat_alphas, params)

        if self.parameters.spin_polarised:
            raise NotImplementedError('Need to check implementation')

        return alphas

    def run(self) -> None:
        '''
        This function runs a KI/pKIPZ/KIPZ workflow from start to finish

        Running this function will generate several directories:
            init/                 -- the density and manifold initialisation calculations
            calc_alpha/orbital_#/ -- calculations where we have fixed a particular orbital
                                     in order to calculate alpha
            final/                -- the final KI/KIPZ calculation
        '''

        # Removing old directories
        if self.parameters.from_scratch:
            if self.parameters.init_orbitals != 'from old ki':
                # if self.init_orbitals == "from old ki" we don't want to delete the directory containing
                # the KI calculation we're reading the manifold from, or the TMP files
                utils.system_call('rm -r init 2>/dev/null', False)
                utils.system_call(f'rm -r {self.master_calc_params["kcp"].outdir} 2>/dev/null', False)
            utils.system_call('rm -r calc_alpha 2>/dev/null', False)
            utils.system_call('rm -r final 2>/dev/null', False)
            if getattr(self, 'redo_preexisting_smooth_dft_calcs', True):
                utils.system_call('rm -r postproc 2>/dev/null', False)

        self.print('Initialisation of density and variational orbitals', style='heading')
        self.perform_initialisation()

        if self.parameters.from_scratch and self.parameters.init_orbitals != 'from old ki' \
                and self.parameters.fix_spin_contamination \
                and self.parameters.init_orbitals not in ['mlwfs', 'projwfs']:
            self.print('Copying the spin-up variational orbitals over to the spin-down channel')
            calc = self.calculations[-1]
            savedir = f'{calc.parameters.outdir}/{calc.parameters.prefix}_{calc.parameters.ndw}.save/K00001'
            utils.system_call(f'cp {savedir}/evc01.dat {savedir}/evc02.dat')
            if calc.parameters.empty_states_nbnd is not None and calc.parameters.empty_states_nbnd > 0:
                utils.system_call(
                    f'cp {savedir}/evc0_empty1.dat {savedir}/evc0_empty2.dat')

        self.print('Calculating screening parameters', style='heading')
        if self.parameters.calculate_alpha:
            self.perform_alpha_calculations()
        else:
            self.print('Skipping calculation of screening parameters', end='')
            if len(self.bands.alpha_history) == 0:
                self.print('; reading values from file')
                self.bands.alphas = self.read_alphas_from_file()
            else:
                self.print()
            self.bands.print_history(indent=self.print_indent + 1)

        # Final calculation
        self.print(f'Final {self.parameters.functional.upper().replace("PK","pK")} calculation', style='heading')
        self.perform_final_calculations()

        # Postprocessing
        if self.parameters.periodic and self.kpath is not None:
            self.print(f'\nPostprocessing', style='heading')
            self.perform_postprocessing()

    def perform_initialisation(self) -> None:
        # Import these here so that if these have been monkey-patched, we get the monkey-patched version
        from koopmans.workflows import WannierizeWorkflow, FoldToSupercellWorkflow

        # The final calculation during the initialisation, regardless of the workflow settings, should write to ndw = 51
        ndw_final = 51

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs']:
            # Wannier functions using pw.x, wannier90.x and pw2wannier90.x
            wannier_workflow = WannierizeWorkflow(**self.wf_kwargs)

            # Perform the wannierisation workflow within the init directory
            self.run_subworkflow(wannier_workflow, subdirectory='init')

            # Now, convert the files over from w90 format to (k)cp format
            fold_workflow = FoldToSupercellWorkflow(**self.wf_kwargs)

            # Do this in the same directory as the wannierisation
            self.run_subworkflow(fold_workflow, subdirectory='init/wannier')

            # Convert self.atoms to the supercell
            self.primitive_to_supercell()

            # We need a dummy calc before the real dft_init in order
            # to copy the previously calculated Wannier functions
            calc = self.new_kcp_calculator('dft_dummy')
            calc.directory = Path('init')
            self.run_calculator(calc, enforce_ss=False)

            # DFT restarting from Wannier functions (after copying the Wannier functions)
            calc = self.new_kcp_calculator('dft_init', restart_mode='restart',
                                           restart_from_wannier_pwscf=True, do_outerloop=True, ndw=ndw_final)
            calc.directory = Path('init')
            restart_dir = Path(f'{calc.parameters.outdir}/{calc.parameters.prefix}_{calc.parameters.ndr}.save/K00001')
            for typ in ['occ', 'emp']:
                if typ == 'occ':
                    evcw_file = Path('init/wannier/occ/evcw.dat')
                    dest_file = restart_dir / 'evc_occupied.dat'
                    if evcw_file.is_file():
                        shutil.copy(evcw_file, dest_file)
                    else:
                        raise OSError(f'Could not find {evcw_file}')
                if typ == 'emp':
                    for i_spin in ['1', '2']:
                        evcw_file = Path('init/wannier/emp') / f'evcw{i_spin}.dat'
                        dest_file = restart_dir / f'evc0_empty{i_spin}.dat'
                        if evcw_file.is_file():
                            shutil.copy(evcw_file, dest_file)
                        else:
                            raise OSError(f'Could not find {evcw_file}')

            self.run_calculator(calc, enforce_ss=False)

            # Check the consistency between the PW and CP band gaps
            pw_calc = [c for c in self.calculations if isinstance(
                c, calculators.PWCalculator) and c.parameters.calculation == 'nscf'][-1]
            pw_gap = pw_calc.results['lumo_ene'] - pw_calc.results['homo_ene']
            cp_gap = calc.results['lumo_energy'] - calc.results['homo_energy']
            if abs(pw_gap - cp_gap) > 2e-2 * pw_gap:
                raise ValueError(f'PW and CP band gaps are not consistent: {pw_gap} {cp_gap}')

            # The CP restarting from Wannier functions must be already converged
            Eini = calc.results['convergence']['filled'][0]['Etot']
            Efin = calc.results['energy']
            if abs(Efin - Eini) > 1e-6 * abs(Efin):
                raise ValueError(f'Too much difference between the initial and final CP energies: {Eini} {Efin}')

        elif self.parameters.init_orbitals == 'from old ki':
            self.print('Copying the density and orbitals from a pre-existing KI calculation')

            # Read the .cpi file to work out the value for ndw
            calc = calculators.KoopmansCPCalculator.fromfile('init/ki_init')

            # Move the old save directory to correspond to ndw_final, using pz_innerloop_init to work out where the
            # code will expect the tmp files to be
            old_savedir = Path(calc.parameters.outdir) / f'{calc.parameters.prefix}_{calc.parameters.ndw}.save'
            savedir = Path(f'{self.new_kcp_calculator("pz_innerloop_init").parameters.outdir}') / \
                f'{calc.parameters.prefix}_{ndw_final}.save'
            if not old_savedir.is_dir():
                raise ValueError(f'{old_savedir} does not exist; a previous '
                                 'and complete KI calculation is required '
                                 'if init_orbitals="from old ki"')
            if savedir.is_dir():
                savedir.rmdir()

            # Use the chdir construct in order to create the directory savedir if it does not already exist and is
            # nested
            with utils.chdir(savedir):
                utils.system_call(f'rsync -a {old_savedir}/ .')

            # Check that the files defining the variational orbitals exist
            savedir /= 'K00001'
            files_to_check = [Path('init/ki_init.cpo'), savedir / 'evc01.dat', savedir / 'evc02.dat']

            if calc.parameters.empty_states_nbnd is not None and calc.parameters.empty_states_nbnd > 0:
                files_to_check += [savedir / f'evc0_empty{i}.dat' for i in [1, 2]]

            for fname in files_to_check:
                if not fname.is_file():
                    raise ValueError(f'Could not find {fname}')

            if not calc.is_complete():
                raise ValueError('init/ki_init.cpo is incomplete so cannot be used '
                                 'to initialise the density and orbitals')

            self.calculations.append(calc)

        elif self.parameters.functional in ['ki', 'pkipz']:
            calc = self.new_kcp_calculator('dft_init')
            calc.directory = Path('init')
            self.run_calculator(calc, enforce_ss=self.parameters.fix_spin_contamination)

            # Use the KS eigenfunctions as better guesses for the variational orbitals
            self._overwrite_canonical_with_variational_orbitals(calc)

            if self.parameters.init_orbitals == 'kohn-sham':
                self._copy_most_recent_calc_to_ndw(ndw_final)
            elif self.parameters.init_orbitals == 'pz':
                calc = self.new_kcp_calculator('pz_innerloop_init', alphas=self.bands.alphas, ndw=ndw_final)
                calc.directory = Path('init')
                if self.calculations[-1].parameters.nelec == 2:
                    # If we only have two electrons, then the filled manifold is trivially invariant under unitary
                    # transformations. Furthermore, the PZ functional is invariant w.r.t. unitary rotations of the
                    # empty states. Thus in this instance we can skip the initialisation of the manifold entirely
                    self.print('Skipping the optimisation of the variational orbitals since they are invariant under '
                               'unitary transformations')
                    self._copy_most_recent_calc_to_ndw(ndw_final)
                else:
                    self.run_calculator(calc)
            else:
                raise ValueError('Should not arrive here')

        elif self.parameters.functional == 'kipz':
            # DFT from scratch
            calc = self.new_kcp_calculator('dft_init')
            calc.directory = Path('init')
            self.run_calculator(calc, enforce_ss=self.parameters.fix_spin_contamination)

            if self.parameters.init_orbitals == 'kohn-sham':
                # Initialise the density with DFT and use the KS eigenfunctions as guesses for the variational orbitals
                self._overwrite_canonical_with_variational_orbitals(calc)
                self._copy_most_recent_calc_to_ndw(ndw_final)
            elif self.parameters.init_orbitals == 'pz':
                # PZ from DFT (generating PZ density and PZ orbitals)
                calc = self.new_kcp_calculator('pz_init', ndw=ndw_final)
                calc.directory = Path('init')
                self.run_calculator(calc)
            else:
                raise ValueError('Should not arrive here')

        else:
            raise ValueError("Should not arrive here; there must be an inconsistency between the above code and \
                             workflow.valid_settings")

        return

    def _copy_most_recent_calc_to_ndw(self, ndw):
        calc = self.calculations[-1]
        if calc.parameters.ndw != ndw:
            assert calc.is_complete(), 'Cannot copy results of a previous calculation that is not itself complete'
            save_prefix = f'{calc.parameters.outdir}/{calc.parameters.prefix}'
            utils.system_call(f'cp -r {save_prefix}_{calc.parameters.ndw}.save {save_prefix}_{ndw}.save')

    def _overwrite_canonical_with_variational_orbitals(self, calc: calculators.CalculatorExt) -> None:
        self.print('Overwriting the variational orbitals with Kohn-Sham orbitals')
        savedir = f'{calc.parameters.outdir}/{calc.parameters.prefix}_{calc.parameters.ndw}.save/K00001'
        utils.system_call(f'cp {savedir}/evc1.dat {savedir}/evc01.dat')
        utils.system_call(f'cp {savedir}/evc2.dat {savedir}/evc02.dat')
        if calc.parameters.empty_states_nbnd is not None and calc.parameters.empty_states_nbnd > 0:
            utils.system_call(f'cp {savedir}/evc_empty1.dat {savedir}/evc0_empty1.dat')
            utils.system_call(f'cp {savedir}/evc_empty2.dat {savedir}/evc0_empty2.dat')

    def perform_alpha_calculations(self) -> None:
        # Set up directories
        Path('calc_alpha').mkdir(exist_ok=True)

        converged = False
        i_sc = 0

        alpha_indep_calcs = []

        while not converged and i_sc < self.parameters.n_max_sc_steps:
            i_sc += 1

            # Setting up directories
            iteration_directory = Path('calc_alpha')
            outdir = self.master_calc_params['kcp'].outdir.name
            outdir = Path.cwd() / iteration_directory / outdir

            if not outdir.is_dir():
                outdir.mkdir()

            if self.parameters.n_max_sc_steps > 1:
                self.print('SC iteration {}'.format(i_sc), style='subheading')
                iteration_directory /= f'iteration_{i_sc}'
                if not iteration_directory.is_dir():
                    iteration_directory.mkdir()

            # Do a KI/KIPZ calculation with the updated alpha values
            trial_calc = self.new_kcp_calculator(calc_presets=self.parameters.functional.replace('pkipz', 'ki'),
                                                 alphas=self.bands.alphas)
            trial_calc.directory = iteration_directory

            if i_sc == 1:
                if self.parameters.functional == 'kipz' and not self.parameters.periodic:
                    # For the first KIPZ trial calculation, do the innerloop
                    trial_calc.parameters.do_innerloop = True
            else:
                # For later SC loops, read in the matching calculation from the
                # previous loop rather than the initialisation calculations
                trial_calc.parameters.ndr = trial_calc.parameters.ndw

            # Run the calculation and store the result. Note that we only need to continue
            # enforcing the spin symmetry if the density will change
            self.run_calculator(trial_calc, enforce_ss=self.parameters.fix_spin_contamination and i_sc > 1)
            alpha_dep_calcs = [trial_calc]

            # Update the bands' self-Hartree and energies (assuming spin-symmetry)
            self.bands.self_hartrees = trial_calc.results['orbital_data']['self-Hartree']

            # Group the bands
            self.bands.assign_groups(allow_reassignment=True)

            skipped_orbitals = []
            first_band_of_each_channel = [self.bands.get(spin=spin)[0] for spin in range(2)]
            # Loop over removing/adding an electron from/to each orbital
            for band in self.bands:
                if self.parameters.spin_polarised and band in first_band_of_each_channel:
                    self.print(f'Spin {band.spin + 1}', style='subheading')

                # Working out what to print for the orbital heading (grouping skipped bands together)
                if band in self.bands.to_solve or band == self.bands.get(spin=band.spin)[-1]:
                    if band not in self.bands.to_solve and (self.parameters.spin_polarised or band.spin == 0):
                        skipped_orbitals.append(band.index)
                    if len(skipped_orbitals) > 0:
                        if len(skipped_orbitals) == 1:
                            self.print(f'Orbital {skipped_orbitals[0]}', style='subheading')
                        else:
                            orb_range = f'{skipped_orbitals[0]}-{skipped_orbitals[-1]}'
                            self.print(f'Orbitals {orb_range}', style='subheading')
                        self.print(f'Skipping; will use the screening parameter of an equivalent orbital')
                        skipped_orbitals = []
                    if band not in self.bands.to_solve:
                        continue
                elif not self.parameters.spin_polarised and band.spin == 1:
                    # In this case, skip over the bands entirely and don't include it in the printout about which
                    # bands we've skipped
                    continue
                else:
                    # Skip the bands which can copy the screening parameter from another
                    # calculation in the same orbital group
                    skipped_orbitals.append(band.index)
                    continue

                self.print(f'Orbital {band.index}', style='subheading')

                # Set up directories
                if self.parameters.spin_polarised:
                    directory = Path(f'{iteration_directory}/spin_{band.spin + 1}/orbital_{band.index}')
                    outdir_band = outdir / f'spin_{band.spin + 1}/orbital_{band.index}'
                else:
                    directory = Path(f'{iteration_directory}/orbital_{band.index}')
                    outdir_band = outdir / f'orbital_{band.index}'
                if not directory.is_dir():
                    directory.mkdir(parents=True)

                # Link tmp files from band-independent calculations
                if not outdir_band.is_dir():
                    outdir_band.mkdir(parents=True)

                    utils.symlink(f'{trial_calc.parameters.outdir}/*.save', outdir_band)

                # Don't repeat if this particular alpha_i was converged
                if i_sc > 1 and abs(band.error) < self.parameters.alpha_conv_thr:
                    self.print(f'Skipping band {band.index} since this alpha is already converged')
                    # if self.parameters.from_scratch:
                    for b in self.bands:
                        if b == band or (band.group is not None and b.group == band.group):
                            b.alpha = band.alpha
                            b.error = band.error
                    continue

                # When we write/update the alpharef files in the work directory
                # make sure to include the fixed band alpha in file_alpharef.txt
                # rather than file_alpharef_empty.txt
                if band.filled:
                    index_empty_to_save = None
                else:
                    if self.parameters.spin_polarised:
                        raise NotImplementedError()
                    index_empty_to_save = band.index - self.bands.num(filled=True, spin=band.spin)

                # Perform the fixed-band-dependent calculations
                if self.parameters.functional in ['ki', 'pkipz']:
                    if band.filled:
                        calc_types = ['dft_n-1']
                    else:
                        calc_types = ['pz_print', 'dft_n+1_dummy', 'dft_n+1']
                else:
                    if band.filled:
                        calc_types = ['kipz_n-1']
                    else:
                        calc_types = ['kipz_print', 'dft_n+1_dummy', 'kipz_n+1']

                for calc_type in calc_types:
                    if self.parameters.functional in ['ki', 'pkipz']:
                        # The calculations whose results change with alpha are...
                        #  - the KI calculations
                        #  - DFT calculations on empty variational orbitals
                        # We don't need to redo any of the others
                        if trial_calc.parameters.empty_states_nbnd == 0 or band.filled:
                            if i_sc > 1 and 'ki' not in calc_type:
                                continue
                    else:
                        # No need to repeat the dummy calculation; all other
                        # calculations are dependent on the screening parameters so
                        # will need updating at each step
                        if i_sc > 1 and calc_type == 'dft_n+1_dummy':
                            continue

                    if 'print' in calc_type:
                        # Note that the 'print' calculations for empty bands do not
                        # in fact involve the fixing of that band (and thus for the
                        # 'fixed' band the corresponding alpha should be in
                        # file_alpharef_empty.txt)
                        alphas = self.bands.alphas
                        filling = self.bands.filling
                    elif not band.filled:
                        # In the case of empty orbitals, we gain an extra orbital in
                        # the spin-up channel, so we explicitly construct both spin
                        # channels for "alphas" and "filling"
                        alphas = self.bands.alphas
                        alphas[band.spin].append(alphas[band.spin][-1])
                        filling = self.bands.filling
                        filling[band.spin][band.index - 1] = True
                        filling[band.spin].append(False)
                    else:
                        alphas = self.bands.alphas
                        filling = self.bands.filling

                    if self.parameters.spin_polarised and (not band.filled or 'print' in calc_type):
                        raise NotImplementedError()

                    # Work out the index of the band that is fixed (noting that we will be throwing away all empty
                    # bands)
                    fixed_band = min(band.index, self.bands.num(filled=True, spin=band.spin) + 1)
                    if self.parameters.spin_polarised and band.spin == 1:
                        fixed_band += self.bands.num(filled=True, spin=0)

                    # Set up calculator
                    calc = self.new_kcp_calculator(calc_type, alphas=alphas, filling=filling, fixed_band=fixed_band,
                                                   index_empty_to_save=index_empty_to_save, outdir=outdir_band)
                    calc.directory = directory

                    # Run kcp.x
                    self.run_calculator(calc)

                    # Reset the value of 'fixed_band' so we can keep track of which calculation
                    # is which. This is important for empty orbital calculations, where fixed_band
                    # is always set to the LUMO but in reality we're fixing the band corresponding
                    # to index_empty_to_save from an earlier calculation
                    calc.parameters.fixed_band = band

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
                            if trial_calc.parameters.empty_states_nbnd > 0 and not band.filled:
                                alpha_dep_calcs.append(calc)
                            else:
                                alpha_indep_calcs.append(calc)
                        else:
                            # For KIPZ, the DFT calculations are dependent on alpha via
                            # the definition of the variational orbitals. We only want to
                            # store the calculations that used the most recent value of alpha

                            alpha_dep_calcs.append(calc)

                    # Copying of evcfixed_empty.dat to evc_occupied.dat
                    if calc_type in ['pz_print', 'kipz_print']:
                        evcempty_dir = f'{outdir_band}/{calc.parameters.prefix}_{calc.parameters.ndw}.save/K00001/'
                    elif calc_type == 'dft_n+1_dummy':
                        evcocc_dir = f'{outdir_band}/{calc.parameters.prefix}_{calc.parameters.ndr}.save/K00001/'
                        if os.path.isfile(f'{evcempty_dir}/evcfixed_empty.dat'):
                            utils.system_call(f'cp {evcempty_dir}/evcfixed_empty.dat {evcocc_dir}/evc_occupied.dat')
                        else:
                            raise OSError(f'Could not find {evcempty_dir}/evcfixed_empty.dat')

                # Calculate an updated alpha and a measure of the error
                # E(N) - E_i(N - 1) - lambda^alpha_ii(1)     (filled)
                # E_i(N + 1) - E(N) - lambda^alpha_ii(0)     (empty)
                #
                # Note that we can do this even from calculations that have been skipped because
                # we read in all the requisite information from the output files and .pkl files
                # that do not get overwritten

                calcs = [c for calc_set in [alpha_dep_calcs, alpha_indep_calcs]
                         for c in calc_set if c.parameters.fixed_band == band]

                alpha, error = self.calculate_alpha_from_list_of_calcs(
                    calcs, trial_calc, band, filled=band.filled)

                for b in self.bands:
                    if b == band or (b.group is not None and b.group == band.group):
                        b.alpha = alpha
                        b.error = error

            self.bands.print_history(indent=self.print_indent + 1)

            converged = all([abs(b.error) < 1e-3 for b in self.bands])

        if converged:
            self.print('Screening parameters have been converged')
        else:
            self.print('Screening parameters have been determined but are not necessarily converged')

    def perform_final_calculations(self) -> None:

        directory = Path('final')
        if not directory.is_dir():
            directory.mkdir()

        if self.parameters.functional == 'pkipz':
            final_calc_types = ['ki', 'pkipz']
        else:
            final_calc_types = [self.parameters.functional]

        for final_calc_type in final_calc_types:

            final_calc_type += '_final'

            # For pKIPZ, the appropriate ndr can change but it is always ndw of the previous
            # KI calculation
            if final_calc_type == 'pkipz_final':
                ndr = [c.parameters.ndw for c in self.calculations if c.prefix in [
                    'ki', 'ki_final'] and hasattr(c.parameters, 'ndw')][-1]
                calc = self.new_kcp_calculator(final_calc_type, ndr=ndr, write_hr=True)
            else:
                calc = self.new_kcp_calculator(final_calc_type, write_hr=True)

            calc.directory = directory

            self.run_calculator(calc)

    def perform_postprocessing(self) -> None:
        # Import these here so that if these have been monkey-patched, we get the monkey-patched version
        from koopmans.workflows import WannierizeWorkflow

        # Transform self.atoms back to the primitive cell
        self.supercell_to_primitive()

        calc: calculators.UnfoldAndInterpolateCalculator

        if self.master_calc_params['ui'].do_smooth_interpolation:
            wf_kwargs = self.wf_kwargs
            wf_kwargs['kgrid'] = [x * y for x,
                                  y in zip(wf_kwargs['kgrid'], self.master_calc_params['ui'].smooth_int_factor)]
            wannier_workflow = WannierizeWorkflow(**wf_kwargs)

            # Here, we allow for skipping of the smooth dft calcs (assuming they have been already run)
            # This is achieved via the optional argument of from_scratch in run_subworkflow(), which
            # overrides the value of wannier_workflow.from_scratch, as well as preventing the inheritance of
            # self.from_scratch to wannier_workflow.from_scratch and back again after the subworkflow finishes
            from_scratch = getattr(self, 'redo_preexisting_smooth_dft_calcs', None)
            self.run_subworkflow(wannier_workflow, subdirectory='postproc', from_scratch=from_scratch)

        for calc_presets in ['occ', 'emp']:
            calc = self.new_ui_calculator(calc_presets)
            self.run_calculator(calc, enforce_ss=False)

        # Merge the two calculations to print out the DOS and bands
        calc = self.new_ui_calculator('merge')

        # Merge the bands
        energies = [c.results['band structure'].energies for c in self.calculations[-2:]]
        reference = np.max(energies[0])
        calc.results['band structure'] = BandStructure(self.kpath, np.concatenate(energies, axis=2) - reference)

        if calc.parameters.do_dos:
            # Generate the DOS
            calc.calc_dos()

        # Print out the merged bands and DOS
        if self.parameters.from_scratch:
            with utils.chdir('postproc'):
                calc.write_results()

        # Store the calculator in the workflow's list of all the calculators
        self.calculations.append(calc)

    def new_kcp_calculator(self, calc_presets: str = 'dft_init',
                           alphas: Optional[List[List[float]]] = None,
                           filling: Optional[List[List[bool]]] = None,
                           **kwargs) -> calculators.KoopmansCPCalculator:
        """

        Generates a new KCP calculator based on the self.master_calc_params["kcp"]
        parameters, modifying the appropriate settings to match the
        chosen calc_presets, and altering any Quantum Espresso keywords
        specified as kwargs

        Arguments:

            calc_presets
                The set of preset values to use; must be one of the following strings:

                Initialisation
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
            alphas = self.bands.alphas

        # Generate a new kcp calculator copied from the master calculator
        calc: calculators.KoopmansCPCalculator = self.new_calculator('kcp', alphas=alphas, filling=filling)

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
            if self.parameters.calculate_alpha:
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
            calc.parameters.nelup += 1
            if 'dummy' not in calc.prefix:
                calc.parameters.restart_from_wannier_pwscf = True

        # electrons
        # For all calculations calculating alpha, remove the empty states and
        # increase the energy thresholds
        if not any([s in calc.prefix for s in ['init', 'print', 'final']]) and \
           calc.prefix not in ['ki', 'kipz', 'dft_dummy']:
            calc.parameters.empty_states_nbnd = 0
            calc.parameters.conv_thr *= 100
            calc.parameters.esic_conv_thr *= 100

        if self.parameters.periodic and not any([s == calc.prefix for s in ['dft_init', 'dft_n-1', 'dft_n+1',
                                                                            'kipz', 'kipz_n-1', 'kipz_n+1']]):
            calc.parameters.do_outerloop = False
            calc.parameters.do_innerloop = False
        elif any([s in calc.prefix for s in ['frozen', 'dummy', 'print', 'innerloop']]) or calc.prefix == 'pkipz_final':
            calc.parameters.do_outerloop = False
            if calc.parameters.empty_states_nbnd > 0:
                calc.parameters.do_outerloop_empty = False
        elif calc.prefix in ['ki', 'ki_final']:
            calc.parameters.do_outerloop = False
            if calc.parameters.empty_states_nbnd > 0:
                if self.parameters.init_empty_orbitals == 'pz':
                    calc.parameters.do_outerloop_empty = True
                else:
                    calc.parameters.do_outerloop_empty = False
        else:
            calc.parameters.do_outerloop = True
            if calc.parameters.empty_states_nbnd > 0:
                calc.parameters.do_outerloop_empty = True

        if calc.parameters.maxiter is None and calc.parameters.do_outerloop:
            calc.parameters.maxiter = 300
        if calc.parameters.empty_states_maxstep is None and calc.parameters.do_outerloop_empty:
            calc.parameters.empty_states_maxstep = 300

        # No empty states minimization in the solids workflow for the moment
        if self.parameters.periodic and calc.parameters.empty_states_nbnd > 0:
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
        if calc.parameters.empty_states_nbnd > 0:
            calc.parameters.do_innerloop_empty = False
        if 'kipz' in calc.prefix:
            calc.parameters.which_orbdep = 'nkipz'
        elif 'pz' in calc.prefix:
            calc.parameters.which_orbdep = 'pz'
        elif 'ki' in calc.prefix:
            calc.parameters.which_orbdep = 'nki'
        if 'print' in calc.prefix:
            calc.parameters.print_wfc_anion = True

        if self.parameters.mt_correction:
            calc.parameters.which_compensation = 'tcc'
        else:
            calc.parameters.which_compensation = 'none'

        # If we are using frozen orbitals, we override the above logic and freeze the variational orbitals
        # post-initialisation
        if self.parameters.frozen_orbitals and 'init' not in calc.prefix and not any([s == calc.prefix for s in
                                                                                      ['dft_n-1', 'dft_n+1', 'kipz_n-1',
                                                                                       'kipz_n+1']]):
            calc.parameters.do_outerloop = False
            calc.parameters.do_innerloop = False
            if calc.parameters.empty_states_nbnd > 0:
                calc.parameters.do_outerloop_empty = False
                calc.parameters.do_innerloop_empty = False

        # Handle any keywords provided by kwargs
        # Note that since this is performed after the above logic this can (deliberately
        # or accidentally) overwrite the above settings
        calc.parameters.update(**kwargs)

        # Sanity checking
        if calc.parameters.print_wfc_anion and calc.parameters.index_empty_to_save is None:
            raise ValueError('Error: print_wfc_anion is set to true but you have not selected '
                             'an index_empty_to_save. Provide this as an argument to new_cp_calculator')

        # avoid innerloops for one-orbital-manifolds
        if calc.parameters.nelup in [0, 1] and calc.parameters.neldw in [0, 1]:
            calc.parameters.do_innerloop = False
        if calc.parameters.empty_states_nbnd == 1:
            calc.parameters.do_innerloop_empty = False

        # don't print QC in some cases
        if 'dummy' in calc.prefix:
            calc.skip_qc = True
        elif calc.prefix[-2:] == '+1':
            # Don't check N+1 energies because they're known to be unreliable
            if 'energy' in calc.results_for_qc:
                calc.results_for_qc.remove('energy')

        return calc

    def new_ui_calculator(self, calc_presets: str, **kwargs) -> calculators.UnfoldAndInterpolateCalculator:

        valid_calc_presets = ['occ', 'emp', 'merge']
        assert calc_presets in valid_calc_presets, 'In KoopmansDSCFWorkflow.new_ui_calculator() calc_presets must be ' \
            '/'.join([f'"{s}"' for s in valid_calc_presets]) + ', but you have tried to set it equal to {calc_presets}'

        if calc_presets == 'merge':
            # Dummy calculator for merging bands and dos
            kwargs['directory'] = Path('postproc')
            pass
        else:
            # Automatically generating UI calculator settings
            kwargs['directory'] = Path(f'postproc/{calc_presets}')
            kwargs['kc_ham_file'] = Path(f'final/ham_{calc_presets}_1.dat').resolve()
            kwargs['w90_seedname'] = Path(f'init/wannier/{calc_presets}/wann').resolve()
            if self.master_calc_params['ui'].do_smooth_interpolation:
                kwargs['dft_smooth_ham_file'] = Path(f'postproc/wannier/{calc_presets}/wann_hr.dat').resolve()
                kwargs['dft_ham_file'] = Path(f'init/wannier/{calc_presets}/wann_hr.dat').resolve()

        calc: calculators.UnfoldAndInterpolateCalculator = self.new_calculator('ui', **kwargs)
        calc.prefix = self.parameters.functional

        return calc

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
                utils.warn('Could not find 2nd order Makov-Payne energy; applying first order only')
                mp_energy = mp1
            else:
                mp_energy = mp1 + mp2

            dE -= np.sign(charge) * mp_energy / self.parameters.eps_inf

        alpha = alpha_guess * (dE - lambda_0) / (lambda_a - lambda_0)

        # The error is lambda^alpha(1) - lambda^alpha_i(1)
        error = dE - lambda_a

        return alpha, error
