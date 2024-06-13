"""

Workflow module for koopmans, containing the workflow for generating maximally localized
Wannier functions (MLWFs), using the Wannier90 code

Written by Riccardo De Gennaro Nov 2020

"""

import copy
import math
import shutil
from pathlib import Path
from typing import List, Optional, TypeVar

# isort: off
import koopmans.mpl_config
import matplotlib.pyplot as plt
# isort: on

import numpy as np
from aiida_koopmans.helpers import get_wannier90bandsworkchain_builder_from_ase

from koopmans import calculators, projections, utils
from koopmans.pseudopotentials import nelec_from_pseudos, read_pseudo_file

from ._workflow import Workflow

CalcExtType = TypeVar('CalcExtType', bound='calculators.CalculatorExt')


class WannierizeWorkflow(Workflow):

    def __init__(self, *args, force_nspin2=False, scf_kgrid=None, **kwargs):
        super().__init__(*args, **kwargs)

        pw_params = self.calculator_parameters['pw']

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs'] \
                and self.parameters.init_empty_orbitals in ['mlwfs', 'projwfs']:

            if self.parameters.spin_polarized:
                spins = ['up', 'down']
            else:
                spins = [None]

            for spin in spins:
                # Work out where we have breaks between blocks of projections, and check that this is commensurate
                # with the number of electrons in this spin channel. Note that this code can be replaced when we have
                # an algorithm for occ-emp separation within Wannier90
                num_bands_occ = self.number_of_electrons(spin)
                if not spin:
                    num_bands_occ /= 2
                divs = self.projections.divisions(spin)
                cumulative_divs = [sum(divs[:i+1]) for i in range(len(divs))]
                if num_bands_occ not in cumulative_divs:
                    message = 'The provided Wannier90 projections are not commensurate with the number of ' \
                              'electrons; divide your list of projections into sublists that represent blocks ' \
                              'of bands to Wannierize separately'
                    raise ValueError(message)

                # Compare the number of bands from PW to Wannier90
                num_bands_w90 = self.projections.num_bands(spin=spin)
                if num_bands_w90 > pw_params.nbnd:
                    raise ValueError(f'You have provided more bands to the Wannier90 calculator ({num_bands_w90}) '
                                     f'than the preceeding PW calculation ({pw_params.nbnd})')
                elif num_bands_w90 == pw_params.nbnd:
                    pass
                else:
                    # Update the projections_blocks to account for additional empty bands
                    self.projections.num_extra_bands[spin] = pw_params.nbnd - num_bands_w90

                # Sanity checking
                assert pw_params.nbnd == self.projections.num_bands(spin=spin)

        elif self.parameters.init_orbitals == 'kohn-sham' and self.parameters.init_empty_orbitals == 'kohn-sham':
            pass

        else:
            raise NotImplementedError('WannierizeWorkflow only supports setting init_orbitals and init_empty_orbitals '
                                      'to "mlwfs"/"projwfs" or "kohn-sham"')

        # Spin-polarization
        self._force_nspin2 = force_nspin2
        if force_nspin2:
            pw_params.nspin = 2
        elif self.parameters.spin_polarized:
            pw_params.nspin = 2
        else:
            pw_params.nspin = 1

        # When running a smooth PW calculation the k-grid for the scf calculation
        # must match the original k-grid
        self._scf_kgrid = scf_kgrid

        # Calculate by default the band structure
        if self.parameters.calculate_bands is None:
            if len(self.kpoints.path.kpts) > 1:
                self.parameters.calculate_bands = True
            else:
                self.parameters.calculate_bands = False

        # This workflow only makes sense for DFT, not an ODD
        self.parameters.functional = 'dft'

    def _run(self):
        '''

        Wrapper for the calculation of (maximally localized) Wannier functions
        using PW and Wannier90

        '''
        # MB mod
        if self.parameters.mode == "ase":
            if self.parameters.init_orbitals in ['mlwfs', 'projwfs']:
                self.print('Wannierization', style='heading')
            else:
                self.print('Kohn-Sham orbitals', style='heading')

            if self.parameters.from_scratch:
                utils.system_call("rm -rf wannier", False)
        else: 
            self.dft_wchains = {}
            self.w90_wchains = {"occ":[],"emp":[]}
            self.wannier90_files = {"occ":{},"emp":{}}

        # Run PW scf and nscf calculations
        # PWscf needs only the valence bands
        calc_pw = self.new_calculator('pw')
        calc_pw.parameters.pop('nbnd', None)
        calc_pw.directory = 'wannier'
        calc_pw.prefix = 'scf'
        if self._scf_kgrid:
            calc_pw.parameters.kpts = self._scf_kgrid
        self.run_calculator(calc_pw)
        
        # MB mod
        if not self.parameters.mode == "ase":
            self.dft_wchains["scf"] = calc_pw.wchain

        calc_pw = self.new_calculator('pw', calculation='nscf', nosym=True, noinv=True)
        calc_pw.directory = 'wannier'
        calc_pw.prefix = 'nscf'
        self.run_calculator(calc_pw)
        # MB mod
        if not self.parameters.mode == "ase":
            self.dft_wchains["nscf"] = calc_pw.wchain

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs'] \
                and self.parameters.init_empty_orbitals in ['mlwfs', 'projwfs']:
            # Loop over the various subblocks that we must wannierize separately
            for block in self.projections:
                n_occ_bands = self.number_of_electrons(block.spin)
                if not block.spin:
                    n_occ_bands /= 2

                if max(block.include_bands) <= n_occ_bands:
                    # Block consists purely of occupied bands
                    init_orbs = self.parameters.init_orbitals
                    # MB mod
                    label_aiida = "occ"
                elif min(block.include_bands) > n_occ_bands:
                    # Block consists purely of empty bands
                    init_orbs = self.parameters.init_empty_orbitals
                    # MB mod
                    label_aiida = "emp"
                else:
                    # Block contains both occupied and empty bands
                    raise ValueError(f'{block} contains both occupied and empty bands. This should not happen.')
                # Store the number of electrons in the ProjectionBlocks object so that it can work out which blocks to
                # merge with one another
                self.projections.num_occ_bands[block.spin] = n_occ_bands

                calc_type = 'w90'
                if block.spin:
                    calc_type += f'_{block.spin}'

                # MB mod: start of construction of W90 workflow and submission
                # Construct the subdirectory label
                w90_dir = Path('wannier') / block.directory

                # MB mod: after the initialization of this calculator, run the Wannier90BandsWorkChain
                # 1) pre-processing Wannier90 calculation
                calc_w90 = self.new_calculator(calc_type, init_orbitals=init_orbs, directory=w90_dir,
                                               **block.w90_kwargs)
                
                # MB mod
                if not self.parameters.mode == "ase":
                    builder = get_wannier90bandsworkchain_builder_from_ase(
                        wannierize_workflow=self, w90_calculator=calc_w90
                    )
                    calc_w90.builder_aiida = builder
                    self.run_calculator(calc_w90)
                    self.w90_wchains[label_aiida].append(calc_w90.wchain)
                else:
                    calc_w90.prefix = 'wann_preproc'
                    calc_w90.command.flags = '-pp'
                    self.run_calculator(calc_w90)
                    utils.system_call(f'rsync -a {calc_w90.directory}/wann_preproc.nnkp {calc_w90.directory}/wann.nnkp')

                    # 2) standard pw2wannier90 calculation
                    calc_p2w = self.new_calculator('pw2wannier', directory=w90_dir,
                                                spin_component=block.spin,
                                                outdir=calc_pw.parameters.outdir)
                    calc_p2w.prefix = 'pw2wan'
                    self.run_calculator(calc_p2w)

                    # 3) Wannier90 calculation
                    calc_w90 = self.new_calculator(calc_type, directory=w90_dir,
                                                init_orbitals=init_orbs,
                                                bands_plot=self.parameters.calculate_bands,
                                                **block.w90_kwargs)
                    calc_w90.prefix = 'wann'
                    self.run_calculator(calc_w90)
                    # MB mod: end W90 workflow and submission
                    
                if hasattr(self, 'bands'):
                    # Add centers and spreads info to self.bands
                    if block.spin is None:
                        remaining_bands = [b for b in self.bands if b.center is None and b.spin == 0]
                    else:
                        if block.spin == 'up':
                            i_spin = 0
                        else:
                            i_spin = 1
                        remaining_bands = [b for b in self.bands if b.center is None and b.spin == i_spin]

                    centers = calc_w90.results['centers']
                    spreads = calc_w90.results['spreads']
                    for band, center, spread in zip(remaining_bands, centers, spreads):
                        band.center = center
                        band.spread = spread

                        if block.spin is None and len(self.bands.get(spin=1)) > 0:
                            # Copy over spin-up results to spin-down
                            [match] = [b for b in self.bands if b.index == band.index and b.spin == 1]
                            match.center = center
                            match.spread = spread

            # Merging Hamiltonian files, U matrix files, centers files if necessary
            if self.parent is not None:
                for merge_directory, block in self.projections.to_merge.items():
                    if len(block) == 1:
                        # MB mod
                        if self.parameters.mode == "ase":
                            with utils.chdir('wannier'):
                                utils.symlink(block[0].directory, merge_directory,
                                            exist_ok=not self.parameters.from_scratch)          
                        else:
                            from aiida_koopmans.data.utils import (
                                produce_wannier90_files,
                            )
                            self.wannier90_files[merge_directory.name] = produce_wannier90_files(self,merge_directory.name)                   
                    else:
                        # MB mod
                        if not self.parameters.mode == "ase":
                            calc_w90.prefix = merge_directory.name 
                            
                        self.merge_wannier_files(block, merge_directory, prefix=calc_w90.prefix)

                        # Extending the U_dis matrix file, if necessary
                        num_wann = sum([b.w90_kwargs['num_wann'] for b in block])
                        num_bands = sum([b.w90_kwargs['num_bands'] for b in block])
                        if num_bands > num_wann and self.parameters.method == 'dfpt':
                            self.extend_wannier_u_dis_file(block, merge_directory, prefix=calc_w90.prefix)
                     
        # MB mod: we can skip this for now.
        if self.parameters.mode == "ase":
            if self.parameters.calculate_bands:
                # Run a "bands" calculation, making sure we don't overwrite
                # the scf/nscf tmp files by setting a different prefix
                calc_pw_bands = self.new_calculator('pw', calculation='bands', kpts=self.kpoints.path)
                calc_pw_bands.directory = 'wannier'
                calc_pw_bands.prefix = 'bands'
                calc_pw_bands.parameters.prefix += '_bands'

                # Link the save directory so that the bands calculation can use the old density
                if self.parameters.from_scratch:
                    [src, dest] = [(c.parameters.outdir / c.parameters.prefix).with_suffix('.save')
                                for c in [calc_pw, calc_pw_bands]]

                    if dest.exists():
                        shutil.rmtree(str(dest))
                    shutil.copytree(src, dest)
                self.run_calculator(calc_pw_bands)

                # Calculate a projected DOS
                pseudos = [read_pseudo_file(calc_pw_bands.parameters.pseudo_dir / p) for p in
                        self.pseudopotentials.values()]
                if all([p['header']['number_of_wfc'] > 0 for p in pseudos]):
                    calc_dos = self.new_calculator('projwfc', filpdos=self.name)
                    calc_dos.directory = 'pdos'
                    calc_dos.pseudopotentials = self.pseudopotentials
                    calc_dos.spin_polarized = self.parameters.spin_polarized
                    calc_dos.pseudo_dir = calc_pw_bands.parameters.pseudo_dir
                    calc_dos.parameters.prefix = calc_pw_bands.parameters.prefix
                    self.run_calculator(calc_dos)

                    # Prepare the DOS for plotting
                    dos = copy.deepcopy(calc_dos.results['dos'])
                else:
                    # Skip if the pseudos don't have the requisite PP_PSWFC blocks
                    utils.warn('Some of the pseudopotentials do not have PP_PSWFC blocks, which means a projected DOS '
                            'calculation is not possible. Skipping...')
                    dos = None

                # Select those calculations that generated a band structure (and are part of this wannierize workflow)
                i_scf = [i for i, c in enumerate(self.calculations) if isinstance(c, calculators.PWCalculator)
                        and c.parameters.calculation == 'scf'][-1]
                selected_calcs = [c for c in self.calculations[i_scf:-1]
                                if 'band structure' in c.results and c != calc_pw_bands]

                # Store the pw BandStructure (for working out the vertical shift to set the valence band edge to zero)
                pw_bands = calc_pw_bands.results['band structure']

                # Prepare the band structures for plotting
                ax = None
                labels = ['explicit'] \
                    + [f'interpolation ({c.directory.name.replace("block_", "block ").replace("spin_", "spin ").replace("_",", ")})'
                    for c in selected_calcs]
                color_cycle = plt.rcParams['axes.prop_cycle']()
                bs_list = []
                bsplot_kwargs_list = []
                colors = {}
                for calc, label in zip([calc_pw_bands] + selected_calcs, labels):
                    if 'band structure' in calc.results:
                        # Load the bandstructure, shifted by the valence band maximum of the pw bands calculation
                        bs = calc.results['band structure'].subtract_reference(pw_bands.reference)

                        # Tweaking the plot aesthetics
                        kwargs = {'label': label}
                        up_label = label.replace(', down', ', up')
                        if ', down' in label:
                            kwargs['ls'] = '--'
                        if ', down' in label and up_label in colors:
                            colors[label] = colors[up_label]
                        else:
                            colors[label] = [next(color_cycle)['color'] for _ in range(bs.energies.shape[0])]
                        if 'explicit' in label:
                            kwargs['ls'] = 'none'
                            kwargs['marker'] = 'x'
                        kwargs['colors'] = colors[label]

                        # Store
                        bs_list.append(bs)
                        bsplot_kwargs_list.append(kwargs)

                # Shift the DOS, too
                if dos is not None:
                    dos._energies -= pw_bands.reference

                # Plot
                self.plot_bandstructure(bs_list, dos, bsplot_kwargs=bsplot_kwargs_list)

        return

    def new_calculator(self, calc_type, *args, **kwargs) -> CalcExtType:  # type: ignore[type-var, misc]
        init_orbs = kwargs.pop('init_orbitals', None)
        calc: CalcExtType = super().new_calculator(calc_type, *args, **kwargs)

        # Extra tweaks for Wannier90 calculations
        if calc_type.startswith('w90'):
            if calc.parameters.num_bands != calc.parameters.num_wann:
                if calc.parameters.dis_num_iter is None:
                    calc.parameters.dis_num_iter = 5000
            else:
                calc.parameters.dis_win_min = None
                calc.parameters.dis_win_max = None
                calc.parameters.dis_froz_min = None
                calc.parameters.dis_froz_max = None

            # mlwfs/projwfs
            if init_orbs == 'projwfs':
                calc.parameters.num_iter = 0
            elif init_orbs == 'mlwfs':
                pass
            else:
                raise ValueError(f'Unrecognized orbital type {init_orbs} (must be "mlwfs" or "projwfs")')

            if calc.parameters.gamma_only != self.kpoints.gamma_only:
                # forcing W90 to follow the same logic of PW for the gamma_trick
                calc.parameters.gamma_only = self.kpoints.gamma_only
        if calc_type == 'pw2wannier':
            if self._force_nspin2 and not self.parameters.spin_polarized:
                calc.parameters.spin_component = 'up'

        return calc

    def merge_wannier_files(self, block: List[projections.ProjectionBlock], merge_directory: Path, prefix: str = 'wann'):
        """
        Merges the hr (Hamiltonian), u (rotation matrix), and wannier centers files of a collection of blocks that
        share the same filling and spin
        """
        # MB mod
        if self.parameters.mode == "ase":
            # Working out the directories where to read in files and where to write out files to
            dirs_in: List[Path] = []
            for b in block:
                assert b.directory is not None, 'The block which you are trying to merge is missing a directory; this ' \
                    'should not happen'
                dirs_in.append(Path('wannier') / b.directory)
            assert merge_directory is not None, 'The block which you are trying to merge is missing a ' \
                'merge_directory; this should not happen'
            dir_out = Path('wannier') / merge_directory

            # Merging the hr (Hamiltonian) files
            self.merge_wannier_hr_files(dirs_in, dir_out, prefix)

            if self.parameters.method == 'dfpt':
                # Merging the U (rotation matrix) files
                self.merge_wannier_u_files(dirs_in, dir_out, prefix)

                # Merging the wannier centers files
                self.merge_wannier_centers_files(dirs_in, dir_out, prefix)
        else:
            # Merging the hr (Hamiltonian) files
            self.wannier90_files[merge_directory.name]["hr_dat"] = self.merge_wannier_hr_files(dirs_in, dir_out, prefix=merge_directory.name, aiida_calculations = self.w90_wchains[merge_directory.name])
            
            if self.parameters.method == 'dfpt':
                # Merging the U (rotation matrix) files
                self.wannier90_files[merge_directory.name]["u_mat"] = self.merge_wannier_u_files(dirs_in, dir_out, prefix=merge_directory.name, aiida_calculations = self.w90_wchains[merge_directory.name])

                # Merging the wannier centers files
                self.wannier90_files[merge_directory.name]["centres_xyz"] = self.merge_wannier_centers_files(dirs_in, dir_out, prefix=merge_directory.name, aiida_calculations = self.w90_wchains[merge_directory.name])
                
            

    @staticmethod
    def merge_wannier_hr_files(dirs_in: List[Path], dir_out: Path, prefix: str, aiida_calculations: list = None):
        # Reading in each hr file in turn
        hr_list = []
        weights_out = None
        rvect_out = None
        
        # MB mod - large changes, I add a long block almost duplicated. We should put this elsewhere.
        if aiida_calculations:
            # READING THE AIIDA FILES
            list_file_in = []
            for block in aiida_calculations:
                hr = block.outputs.wannier90.retrieved.get_object_content('aiida' + '_hr.dat')
                hr = hr.split("\n")
                for line in range(len(hr)):
                    hr[line] += "\n"
                list_file_in.append(hr[:-1])
            
            for file in list_file_in:
                # Reading the hr file
                hr, rvect, weights, nrpts = utils.read_wannier_hr_file(file=file)

                # Sanity checking
                if weights_out is None:
                    weights_out = weights
                elif weights != weights_out:
                    raise ValueError(f'{file} contains weights that differ from the other blocks. This should not '
                                        'happen.')
                if rvect_out is None:
                    rvect_out = rvect
                elif np.all(rvect != rvect_out):
                    raise ValueError(f'{file} contains a set of R-vectors that differ from the other blocks. This '
                                        'should not happen.')

                # Reshaping this block of the Hamiltonian in preparation for constructing the block matrix, and storing it
                num_wann2 = hr.size // nrpts
                num_wann = int(math.sqrt(num_wann2))
                hr_list.append(hr.reshape(nrpts, num_wann, num_wann))

            # Constructing the block matrix hr_out which is dimensions (nrpts, num_wann_tot, num_wann_tot)
            num_wann_tot = sum([hr.shape[-1] for hr in hr_list])
            hr_out = np.zeros((nrpts, num_wann_tot, num_wann_tot), dtype=complex)
            start = 0
            for hr in hr_list:
                end = start + hr.shape[1]
                for irpt in range(nrpts):
                    hr_out[irpt, start:end, start:end] = hr[irpt, :, :]
                start = end

            assert rvect_out is not None
            assert weights_out is not None

            hr = utils.write_wannier_hr_file(dir_out / (prefix + '_hr.dat'), hr_out, rvect_out.tolist(), weights_out, "aiida_hr.dat")
            
            return hr # SingleFileData.
            
        for dir_in in dirs_in:
            # Reading the hr file
            fname_in = dir_in / (prefix + '_hr.dat')
            hr, rvect, weights, nrpts = utils.read_wannier_hr_file(fname_in)

            # Sanity checking
            if weights_out is None:
                weights_out = weights
            elif weights != weights_out:
                raise ValueError(f'{fname_in} contains weights that differ from the other blocks. This should not '
                                 'happen.')
            if rvect_out is None:
                rvect_out = rvect
            elif np.all(rvect != rvect_out):
                raise ValueError(f'{fname_in} contains a set of R-vectors that differ from the other blocks. This '
                                 'should not happen.')

            # Reshaping this block of the Hamiltonian in preparation for constructing the block matrix, and storing it
            num_wann2 = hr.size // nrpts
            num_wann = int(math.sqrt(num_wann2))
            hr_list.append(hr.reshape(nrpts, num_wann, num_wann))

        # Constructing the block matrix hr_out which is dimensions (nrpts, num_wann_tot, num_wann_tot)
        num_wann_tot = sum([hr.shape[-1] for hr in hr_list])
        hr_out = np.zeros((nrpts, num_wann_tot, num_wann_tot), dtype=complex)
        start = 0
        for hr in hr_list:
            end = start + hr.shape[1]
            for irpt in range(nrpts):
                hr_out[irpt, start:end, start:end] = hr[irpt, :, :]
            start = end

        assert rvect_out is not None
        assert weights_out is not None

        utils.write_wannier_hr_file(dir_out / (prefix + '_hr.dat'), hr_out, rvect_out.tolist(), weights_out)

    @staticmethod
    def merge_wannier_u_files(dirs_in: List[Path], dir_out: Path, prefix: str, aiida_calculations: list = None):
        u_list = []
        kpts_master = None
        
        # MB mod - large changes, I add a long block almost duplicated. We should put this elsewhere.
        if aiida_calculations:
            # READING THE AIIDA FILES
            list_file_in = []
            for block in aiida_calculations:
                hr = block.outputs.wannier90.retrieved.get_object_content('aiida' + '_u.mat')
                hr = hr.split("\n")
                for line in range(len(hr)):
                    hr[line] += "\n"
                list_file_in.append(hr[:-1])
            
            for file in list_file_in:
                # Reading the U file
                fname_in = dir_in / (prefix + '_u.mat')

                umat, kpts, nkpts = utils.read_wannier_u_file(file=file)

                if kpts_master is None:
                    kpts_master = kpts
                elif nkpts == len(kpts_master) and np.allclose(kpts, kpts_master):
                    pass
                else:
                    raise ValueError(f'{fname_in} has an inconsistent set of k-points with the other files you are merging')

                u_list.append(umat)

            shape_u_merged = [nkpts] + np.sum([u.shape for u in u_list], axis=0)[1:].tolist()
            u_merged = np.zeros(shape_u_merged, dtype=complex)

            # Constructing a large block-diagonal U matrix from all the individual matrices in u_list
            i_start = 0
            j_start = 0
            for u in u_list:
                i_end = i_start + u.shape[1]
                j_end = j_start + u.shape[2]
                u_merged[:, i_start:i_end, j_start:j_end] = u
                i_start = i_end
                j_start = j_end

            # Writing out the large U file
            u = utils.write_wannier_u_file(dir_out / (prefix + '_u.mat'), u_merged, kpts)
            
            return u 
        
        for dir_in in dirs_in:
            # Reading the U file
            fname_in = dir_in / (prefix + '_u.mat')

            umat, kpts, nkpts = utils.read_wannier_u_file(fname_in)

            if kpts_master is None:
                kpts_master = kpts
            elif nkpts == len(kpts_master) and np.allclose(kpts, kpts_master):
                pass
            else:
                raise ValueError(f'{fname_in} has an inconsistent set of k-points with the other files you are merging')

            u_list.append(umat)

        shape_u_merged = [nkpts] + np.sum([u.shape for u in u_list], axis=0)[1:].tolist()
        u_merged = np.zeros(shape_u_merged, dtype=complex)

        # Constructing a large block-diagonal U matrix from all the individual matrices in u_list
        i_start = 0
        j_start = 0
        for u in u_list:
            i_end = i_start + u.shape[1]
            j_end = j_start + u.shape[2]
            u_merged[:, i_start:i_end, j_start:j_end] = u
            i_start = i_end
            j_start = j_end

        # Writing out the large U file
        utils.write_wannier_u_file(dir_out / (prefix + '_u.mat'), u_merged, kpts)

    def merge_wannier_centers_files(self, dirs_in: List[Path], dir_out: Path, prefix: str, aiida_calculations: list = None):
        centers_list = []
        
        # MB mod - large changes, I add a long block almost duplicated. We should put this elsewhere.
        if aiida_calculations:
            # READING THE AIIDA FILES
            list_file_in = []
            for block in aiida_calculations:
                hr = block.outputs.wannier90.retrieved.get_object_content('aiida' + '_centers.xyz')
                hr = hr.split("\n")
                for line in range(len(hr)):
                    hr[line] += "\n"
                list_file_in.append(hr[:-1])
            
            for file in list_file_in:
            # Reading the centers file
                fname_in = dir_in / (prefix + '_centres.xyz')

                centers, _ = utils.read_wannier_centers_file(fname_in)

                centers_list += centers

            # Writing the centers file
            centers = utils.write_wannier_centers_file(dir_out / (prefix + '_centres.xyz'), centers_list, self.atoms)
            return centers
        
        for dir_in in dirs_in:
            # Reading the centers file
            fname_in = dir_in / (prefix + '_centres.xyz')

            centers, _ = utils.read_wannier_centers_file(fname_in)

            centers_list += centers

        # Writing the centers file
        utils.write_wannier_centers_file(dir_out / (prefix + '_centres.xyz'), centers_list, self.atoms)

    def extend_wannier_u_dis_file(self, block: List[projections.ProjectionBlock], merge_directory: Path, prefix: str = 'wann'):
        # Read in
        assert block[-1].directory is not None
        
        # MB mod
        if self.parameters.mode == "ase":
            fname_in = Path('wannier') / block[-1].directory / (prefix + '_u_dis.mat')
            udis_mat, kpts, _ = utils.read_wannier_u_file(fname_in)
        else:
            u_dis = self.w90_wchains["emp"][-1].outputs.wannier90.retrieved.get_object_content('aiida' + '_u_dis.mat')
            u_dis = u_dis.split("\n")
            for line in range(len(u_dis)):
                u_dis[line] += "\n"
            u_dis = u_dis[:-1]
            udis_mat, kpts, _ = utils.read_wannier_u_file(fname_in, file=u_dis)

        # Calculate how many empty bands we have
        spin = block[0].spin
        if spin:
            nbnd_occ = self.number_of_electrons(spin)
        else:
            nbnd_occ = self.number_of_electrons() // 2
        nbnd_tot = self.calculator_parameters['pw'].nbnd - nbnd_occ

        # Calculate how many empty wannier functions we have
        nwann_tot = sum([len(p) for p in block])

        # Build up the larger U_dis matrix, which is a nkpts x nwann_emp x nbnd_emp matrix...
        udis_mat_large = np.zeros((len(kpts), nwann_tot, nbnd_tot), dtype=complex)
        # ... with the diagonal entries equal to 1...
        udis_mat_large[:, :nwann_tot, :nwann_tot] = np.identity(nwann_tot)
        # ... except for the last block, where we insert the contents of the corresponding u_dis file
        udis_mat_large[:, -udis_mat.shape[1]:, -udis_mat.shape[2]:] = udis_mat

        # Write out
        fname_out = Path('wannier') / merge_directory / (prefix + '_u_dis.mat')
        # MB mod
        if self.parameters.mode == "ase":
            utils.write_wannier_u_file(fname_out, udis_mat_large, kpts)
        else:
            self.wannier90_files[merge_directory.name]["u_dis_mat"] = utils.write_wannier_u_file(fname_out, udis_mat_large, kpts, "aiida_u_dis.mat")