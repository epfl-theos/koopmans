"""

Workflow module for koopmans, containing the workflow for generating maximally localized
Wannier functions (MLWFs), using the Wannier90 code

Written by Riccardo De Gennaro Nov 2020

"""

import copy
import math
import shutil
from functools import partial
from pathlib import Path
from typing import Callable, List, Optional, Tuple, TypeVar

from ase import Atoms
from pydantic import BaseModel

# isort: off
import koopmans.mpl_config
import matplotlib.pyplot as plt
# isort: on

import numpy as np

from koopmans import calculators, projections, utils
from koopmans.process import Process
from koopmans.pseudopotentials import nelec_from_pseudos, read_pseudo_file

from ._merge_wannier import (MergeProcess, merge_wannier_centers_file_contents,
                             merge_wannier_hr_file_contents,
                             merge_wannier_u_file_contents)
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
        if self.parameters.init_orbitals in ['mlwfs', 'projwfs']:
            self.print('Wannierization', style='heading')
        else:
            self.print('Kohn-Sham orbitals', style='heading')

        if self.parameters.from_scratch:
            utils.system_call("rm -rf wannier", False)

        # Run PW scf and nscf calculations
        # PWscf needs only the valence bands
        calc_scf = self.new_calculator('pw')
        calc_scf.parameters.pop('nbnd', None)
        calc_scf.directory = 'wannier'
        calc_scf.prefix = 'scf'
        if self._scf_kgrid:
            calc_scf.parameters.kpts = self._scf_kgrid
        self.run_calculator(calc_scf)

        calc_nscf = self.new_calculator('pw', calculation='nscf', nosym=True, noinv=True)
        calc_nscf.directory = 'wannier'
        calc_nscf.prefix = 'nscf'
        self.link(calc_scf, calc_scf.parameters.outdir, calc_nscf, calc_nscf.parameters.outdir)
        self.run_calculator(calc_nscf)

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
                elif min(block.include_bands) > n_occ_bands:
                    # Block consists purely of empty bands
                    init_orbs = self.parameters.init_empty_orbitals
                else:
                    # Block contains both occupied and empty bands
                    raise ValueError(f'{block} contains both occupied and empty bands. This should not happen.')
                # Store the number of electrons in the ProjectionBlocks object so that it can work out which blocks to
                # merge with one another
                self.projections.num_occ_bands[block.spin] = n_occ_bands

                calc_type = 'w90'
                if block.spin:
                    calc_type += f'_{block.spin}'

                # Construct the subdirectory label
                w90_dir = Path('wannier') / block.directory

                # 1) pre-processing Wannier90 calculation
                calc_w90 = self.new_calculator(calc_type, init_orbitals=init_orbs, directory=w90_dir,
                                               **block.w90_kwargs)
                calc_w90.prefix = 'wann_preproc'
                calc_w90.command.flags = '-pp'
                self.run_calculator(calc_w90)
                utils.system_call(f'rsync -a {calc_w90.directory}/wann_preproc.nnkp {calc_w90.directory}/wann.nnkp')

                # 2) standard pw2wannier90 calculation
                calc_p2w = self.new_calculator('pw2wannier', directory=w90_dir,
                                               spin_component=block.spin)
                self.link(calc_nscf, calc_nscf.parameters.outdir, calc_p2w, calc_p2w.parameters.outdir)
                calc_p2w.prefix = 'pw2wan'
                self.run_calculator(calc_p2w)

                # 3) Wannier90 calculation
                calc_w90 = self.new_calculator(calc_type, directory=w90_dir,
                                               init_orbitals=init_orbs,
                                               bands_plot=self.parameters.calculate_bands,
                                               **block.w90_kwargs)
                calc_w90.prefix = 'wann'
                self.run_calculator(calc_w90)

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
            if True:  # self.parent is not None:
                utils.warn('manually forcing the merging of the wannier files')
                for merge_directory, block in self.projections.to_merge.items():
                    if len(block) == 1:
                        with utils.chdir('wannier'):
                            utils.symlink(block[0].directory, merge_directory,
                                          exist_ok=not self.parameters.from_scratch)
                    else:
                        self.merge_wannier_files(block, merge_directory, prefix=calc_w90.prefix)

                        # Extending the U_dis matrix file, if necessary
                        num_wann = sum([b.w90_kwargs['num_wann'] for b in block])
                        num_bands = sum([b.w90_kwargs['num_bands'] for b in block])
                        if num_bands > num_wann and self.parameters.method == 'dfpt':
                            raise NotImplementedError()
                            # extend_proc = Process(merge_function=extend_wannier_u_dis_file)
                            # self.extend_wannier_u_dis_file(block, merge_directory, prefix=calc_w90.prefix)

        if self.parameters.calculate_bands:
            # Run a "bands" calculation, making sure we don't overwrite
            # the scf/nscf tmp files by setting a different prefix
            calc_pw_bands = self.new_calculator('pw', calculation='bands', kpts=self.kpoints.path)
            calc_pw_bands.directory = 'wannier'
            calc_pw_bands.prefix = 'bands'
            calc_pw_bands.parameters.prefix += '_bands'

            # Link the save directory so that the bands calculation can use the old density
            self.link(calc_nscf, (calc_nscf.parameters.outdir / calc_nscf.parameters.prefix).with_suffix('.save'),
                      calc_pw_bands, (calc_pw_bands.parameters.outdir / calc_pw_bands.parameters.prefix).with_suffix('.save'))
            self.run_calculator(calc_pw_bands)

            # Calculate a projected DOS
            pseudos = [read_pseudo_file(calc_pw_bands.directory / calc_pw_bands.parameters.pseudo_dir / p) for p in
                       self.pseudopotentials.values()]
            if all([p['header']['number_of_wfc'] > 0 for p in pseudos]) and False:
                calc_dos = self.new_calculator('projwfc', filpdos=self.name)
                calc_dos.directory = 'pdos'
                calc_dos.pseudopotentials = self.pseudopotentials
                calc_dos.spin_polarized = self.parameters.spin_polarized
                calc_dos.pseudo_dir = calc_pw_bands.parameters.pseudo_dir
                calc_dos.parameters.prefix = calc_pw_bands.parameters.prefix
                self.link(calc_pw_bands, calc_pw_bands.parameters.outdir, calc_dos, calc_dos.parameters.outdir)
                self.run_calculator(calc_dos)

                # Prepare the DOS for plotting
                dos = copy.deepcopy(calc_dos.results['dos'])
            else:
                # Skip if the pseudos don't have the requisite PP_PSWFC blocks
                utils.warn('Manually disabling pDOS')
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

        # Working out the calculations from which to fetch the files, where to write out files to
        src_calcs: List[calculators.Wannier90Calculator] = []
        for b in block:
            assert b.directory is not None, 'The block which you are trying to merge is missing a directory; this ' \
                'should not happen'
            matching_calc = [c for c in self.calculations[::-1]
                             if isinstance(c, calculators.Wannier90Calculator) and c.directory.name == b.directory.name][0]
            src_calcs.append(matching_calc)

        assert merge_directory is not None, 'The block which you are trying to merge is missing a ' \
            'merge_directory; this should not happen'

        # Merging the wannier_hr (Hamiltonian) files
        merge_hr_proc = MergeProcess(merge_function=merge_wannier_hr_file_contents,
                                     src_files=[(calc, Path(prefix + '_hr.dat')) for calc in src_calcs],
                                     dst_file=Path('wannier') / merge_directory / (prefix + '_hr.dat'))
        self.run_process(merge_hr_proc)

        if self.parameters.method == 'dfpt' and self.parent is not None:
            # Merging the U (rotation matrix) files
            merge_u_proc = MergeProcess(merge_function=merge_wannier_u_file_contents,
                                        src_files=[(calc, Path(prefix + '_u.mat')) for calc in src_calcs],
                                        dst_file=Path('wannier') / merge_directory / (prefix + '_u.mat'))
            self.run_process(merge_u_proc)

            # Merging the wannier centers files
            merge_centers_proc = MergeProcess(merge_function=partial(merge_wannier_centers_file_contents, atoms=self.atoms),
                                              src_files=[(calc, Path(prefix + '_centres.xyz')) for calc in src_calcs],
                                              dst_file=Path('wannier') / merge_directory / (prefix + '_centres.xyz'))
            self.run_process(merge_centers_proc)
