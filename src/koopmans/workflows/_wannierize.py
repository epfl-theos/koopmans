"""Wannierize workflow module for koopmans.

Workflow module for koopmans, containing the workflow for generating maximally localized
Wannier functions (MLWFs), using the Wannier90 code

Written by Riccardo De Gennaro Nov 2020

"""

import copy
from enum import Enum
from functools import partial
from typing import (Annotated, Any, Dict, List, Optional, Sequence, Type,
                    TypeVar)

import numpy as np
import numpy.typing as npt
from ase_koopmans.dft.kpoints import BandPath
from ase_koopmans.spectrum.band_structure import BandStructure
from ase_koopmans.spectrum.doscollection import GridDOSCollection

# isort: off
import koopmans.mpl_config  # noqa: F401
import matplotlib.pyplot as plt
# isort: on

from koopmans import calculators, utils
from koopmans.calculators import Wannier90Calculator
from koopmans.files import File
from koopmans.process_io import IOModel
from koopmans.processes.wannier import (ExtendProcess, MergeProcess,
                                        extend_wannier_u_dis_file_content,
                                        merge_wannier_centers_file_contents,
                                        merge_wannier_hr_file_contents,
                                        merge_wannier_u_file_contents)
from koopmans.processes.wjl import (WannierJLCheckNeighborsProcess,
                                    WannierJLGenerateNeighborsProcess,
                                    WannierJLSplitProcess)
from koopmans.projections import (BlockID, ExplicitProjections,
                                  ImplicitProjectionsBlock, Projections,
                                  ProjectionsBlock)
from koopmans.status import Status
from koopmans.utils import Spin
from koopmans.utils.warnings import IncommensurateProjectionsWarning, warn

from ._workflow import Workflow

CalcExtType = TypeVar('CalcExtType', bound='calculators.CalculatorExt')


class WannierizeOutput(IOModel):
    """Output model for the WannierizeWorkflow."""

    band_structures: List[BandStructure]
    dos: Optional[GridDOSCollection] = None
    u_files: Dict[BlockID, File | None]
    hr_files: Dict[BlockID, File]
    centers_files: Dict[BlockID, File | None]
    u_dis_files: Dict[BlockID, File | None]
    nnkp_files: Dict[BlockID, File | None]
    nscf_calculation: calculators.PWCalculator
    wannier90_calculations: List[calculators.Wannier90Calculator]
    projections: Projections


class MergeableFile(Enum):
    """Wannier90 output files that can be merged or extended when Wannierizing blocks of bands separately."""

    HR = '_hr.mat'
    U = '_u.mat'
    CENTERS = '_centres.xyz'
    U_DIS = '_u_dis.mat'


class WannierizeWorkflow(Workflow[WannierizeOutput]):
    """Workflow for Wannierizing an entire system using Quantum ESPRESSO and Wannier90.

    The bands are possibly split into blocks of bands (separated in energy) that are Wannierized separately

    """

    output_model = WannierizeOutput

    def __init__(self, *args, force_nspin2: bool = False, scf_kgrid: List[float] | None = None,
                 files_to_merge: List[MergeableFile] = [], merge_occ_and_empty: bool = False, **kwargs):
        super().__init__(*args, **kwargs)

        if self.projections is None:
            raise ValueError(f'{self.__class__.__name__} requires projections but none were provided.')

        pw_params = self.calculator_parameters['pw']

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs'] \
                and self.parameters.init_empty_orbitals in ['mlwfs', 'projwfs']:
            pass
        elif self.parameters.init_orbitals == 'kohn-sham' and self.parameters.init_empty_orbitals == 'kohn-sham':
            pass
        else:
            raise NotImplementedError('`WannierizeWorkflow` only supports setting `init_orbitals` and '
                                      '`init_empty_orbitals` to `mlwfs`/`projwfs`/`kohn-sham`')

        # Sanity checks
        if self.parameters.init_orbitals in ['mlwfs', 'projwfs'] \
                and self.parameters.init_empty_orbitals in ['mlwfs', 'projwfs']:

            spins: List[Spin] = [Spin.UP, Spin.DOWN] if self.parameters.spin_polarized else [Spin.NONE]

            for spin in spins:
                # Work out where we have breaks between blocks of projections, and check that this is commensurate
                # with the number of electrons in this spin channel. Note that this code can be replaced when we have
                # an algorithm for occ-emp separation within Wannier90
                num_bands_occ = self.number_of_electrons(spin)
                if spin == Spin.NONE:
                    num_bands_occ //= 2
                divs = self.projections.divisions(spin)
                cumulative_divs = [sum(divs[:i + 1]) for i in range(len(divs))]
                if num_bands_occ not in cumulative_divs and isinstance(self.projections, ExplicitProjections):
                    message = 'The provided Wannier90 projections are not commensurate with the number of ' \
                              'electrons; divide your list of projections into sublists ' \
                              'that represent blocks of bands to Wannierize separately'
                    warn(message, IncommensurateProjectionsWarning)

                # Compare the number of bands from PW to Wannier90
                num_bands_w90 = self.projections.num_bands(spin=spin)
                if num_bands_w90 > pw_params.nbnd:
                    warn(f'You have provided more bands to the `Wannier90 calculator` ({num_bands_w90}) '
                         f'than the preceeding PW calculation ({pw_params.nbnd})')
                    pw_params.nbnd = num_bands_w90
                elif num_bands_w90 == pw_params.nbnd:
                    pass
                else:
                    # Update the projections_blocks to account for additional empty bands
                    self.projections.num_extra_bands[spin] = pw_params.nbnd - num_bands_w90

                # Sanity checking
                assert pw_params.nbnd == self.projections.num_bands(spin=spin)

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
            assert isinstance(self.kpoints.path, BandPath)
            if len(self.kpoints.path.kpts) > 1:
                self.parameters.calculate_bands = True
            else:
                self.parameters.calculate_bands = False

        # This workflow only makes sense for DFT, not an ODD
        self.parameters.functional = 'dft'

        # Save the attributes that define the behavior of the workflow when merging/extending files
        self._files_to_merge = files_to_merge
        self._merge_occ_and_empty = merge_occ_and_empty

    def _run(self) -> None:
        """Run the workflow."""
        if self.projections is None:
            raise ValueError(f'{self.__class__.__name__} requires projections but none were provided.')

        # Run PW scf and nscf calculations
        # PWscf needs only the valence bands
        calc_scf = self.new_calculator('pw')
        calc_scf.parameters.pop('nbnd', None)
        calc_scf.prefix = 'scf'
        if self._scf_kgrid:
            calc_scf.parameters.kpts = self._scf_kgrid
        status = self.run_steps(calc_scf)
        if status != Status.COMPLETED:
            return

        calc_nscf = self.new_calculator('pw', calculation='nscf', nosym=True, noinv=True)
        calc_nscf.prefix = 'nscf'
        calc_nscf.link(File(calc_scf, calc_scf.parameters.outdir))
        status = self.run_steps(calc_nscf)
        if status != Status.COMPLETED:
            return

        if self.parameters.calculate_bands or self.parameters.block_wannierization_threshold:
            # Run a "bands" calculation
            calc_pw_bands = self.new_calculator('pw', calculation='bands', kpts=self.kpoints.path)
            calc_pw_bands.prefix = 'bands'
            calc_pw_bands.parameters.prefix += '_bands'

            # Link the save directory so that the bands calculation can use the old density
            calc_pw_bands.link(
                File(calc_nscf, (calc_nscf.parameters.outdir / calc_nscf.parameters.prefix).with_suffix('.save')),
                (calc_pw_bands.parameters.outdir / calc_pw_bands.parameters.prefix).with_suffix('.save')
            )
            status = self.run_steps(calc_pw_bands)
            if status != Status.COMPLETED:
                return
        else:
            calc_pw_bands = None

        u_files = {}
        hr_files: Dict[BlockID, File] = {}
        centers_files = {}
        u_dis_files = {}
        nnkp_files = {}
        wannier90_calculations: List[Wannier90Calculator] = []
        new_projections = None

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs'] \
                and self.parameters.init_empty_orbitals in ['mlwfs', 'projwfs']:

            # Store the number of electrons for each spin channel in self.projections. It needs to know this
            # to work out which bands to merge with one another.
            for spin in Spin:
                nelec = self.number_of_electrons(spin)
                if spin == Spin.NONE:
                    nelec //= 2
                self.projections.num_occ_bands[spin] = nelec

            block_subworkflows: list[Workflow] = []
            groups: list[list[int] | None] = []
            for block in self.projections:
                # Check if the block should be split. It will be split if...
                # - the block includes both occupied and empty bands
                # - the are gaps between bands that are wider than `block_wannierization_threshold`
                n_occ_bands = self.number_of_electrons(block.spin)
                if block.spin == Spin.NONE:
                    n_occ_bands //= 2

                if calc_pw_bands is not None:
                    bs = calc_pw_bands.results['band structure']
                    ispin = 1 if block.spin == Spin.DOWN else 0
                    all_groups = detect_band_blocks(bs._energies[ispin, :, :self.projections.num_wann(block.spin)],
                                                    tol=self.parameters.block_wannierization_threshold,
                                                    num_occ_bands=n_occ_bands)

                    # Restrict the groups to only those bands that are in this particular block
                    groups = []
                    for group in all_groups:
                        # Calculate the overlap with the bands of this block
                        overlap = [g for g in group if g in block.include_bands]
                        if overlap:
                            groups.append(overlap)
                else:
                    groups = [None]

                # Construct the subworkflow
                kwargs: dict[str, Any] = {}
                subworkflow_class: Type[Workflow]
                if len(groups) > 1:
                    # Need to split
                    subworkflow_class = WannierizeAndSplitBlockWorkflow
                    kwargs['groups'] = groups
                    kwargs['files_to_merge'] = self._files_to_merge + [MergeableFile.HR]
                    minimize = self.parameters.init_orbitals == 'mlwfs'
                else:
                    subworkflow_class = WannierizeBlockWorkflow
                    if max(block.include_bands) <= n_occ_bands:
                        # Block consists purely of occupied bands
                        minimize = self.parameters.init_orbitals == 'mlwfs'
                    else:
                        # Block consists purely of empty bands
                        minimize = self.parameters.init_empty_orbitals == 'mlwfs'

                subworkflow = subworkflow_class.fromparent(self,
                                                           minimize=minimize,
                                                           force_nspin2=self._force_nspin2,
                                                           block=block,
                                                           pw_outdir=File(calc_nscf, calc_nscf.parameters.outdir),
                                                           **kwargs)
                if len(self.projections) > 1:
                    subworkflow.name = f'Wannierize {block.id.label.replace("_", " ").replace("block", "Block")}'
                else:
                    subworkflow.name = 'Wannierize'
                block_subworkflows.append(subworkflow)

            for wf in block_subworkflows:
                wf.proceed()
            if any([wf.status != Status.COMPLETED for wf in block_subworkflows]):
                return

            # Assemble the outputs of each finalized block in a flattened list
            finalized_blocks_and_outputs = []
            for subwf in block_subworkflows:
                if isinstance(subwf, WannierizeBlockWorkflow):
                    finalized_blocks_and_outputs.append((subwf.block, subwf.outputs))
                elif isinstance(subwf, WannierizeAndSplitBlockWorkflow):
                    for subblock, outputs in zip(subwf.outputs.blocks, subwf.outputs.block_outputs):
                        finalized_blocks_and_outputs.append((subblock, outputs))

                    # Add the unique U and U_dis files
                    u_dis_files[subwf.block.id] = subwf.outputs.u_dis_file
                    u_files[subwf.block.id] = subwf.outputs.u_file
                    centers_files[subwf.block.id] = subwf.outputs.centers_file
                    assert subwf.outputs.hr_file is not None
                    hr_files[subwf.block.id] = subwf.outputs.hr_file

            # Construct a new Projections object with the finalized blocks
            if len(finalized_blocks_and_outputs) != len(self.projections):
                new_projections = self.projections.__class__(
                    blocks=[b[0] for b in finalized_blocks_and_outputs],
                    atoms=self.atoms,
                    num_occ_bands=self.projections.num_occ_bands,
                    exclude_bands=self.projections.exclude_bands,
                )
            else:
                new_projections = self.projections

            # Store the outputs of the finalized blocks
            for block, outputs in finalized_blocks_and_outputs:
                assert outputs.hr_file is not None
                hr_files[block.id] = outputs.hr_file
                centers_files[block.id] = outputs.centers_file
                u_files[block.id] = outputs.u_file
                u_dis_files[block.id] = outputs.u_dis_file
                nnkp_files[block.id] = outputs.nnkp_file
                wannier90_calculations.append(outputs.wannier90_calculation)

            # Merging Hamiltonian files, U matrix files, centers files if necessary
            for merged_block_id, blocks in new_projections.to_merge(
                    merge_occ_and_empty=self._merge_occ_and_empty).items():
                if len(blocks) == 1:
                    # If the set of blocks to merge only contains one block, we don't need to merge anything. In this
                    # case, the "merged" outputs are just the outputs of that single block
                    [block] = blocks
                    hr_files[merged_block_id] = hr_files[block.id]
                    u_files[merged_block_id] = u_files[block.id]
                    u_dis_files[merged_block_id] = u_dis_files[block.id]
                    centers_files[merged_block_id] = centers_files[block.id]

                else:
                    # Fetching the list of calculations for this block
                    src_calcs: List[calculators.Wannier90Calculator] = [
                        b.w90_calc for b in blocks if b.w90_calc is not None]
                    prefix = src_calcs[-1].prefix

                    if MergeableFile.HR in self._files_to_merge:
                        # Merging the wannier_hr (Hamiltonian) files
                        merge_hr_proc = MergeProcess(merge_function=merge_wannier_hr_file_contents,
                                                     src_files=[File(calc, calc.prefix + '_hr.dat')
                                                                for calc in src_calcs],
                                                     dst_file=f'{prefix}_hr.dat')
                        merge_hr_proc.name = f'merge_{merged_block_id.label}_wannier_hamiltonian'
                        status = self.run_steps(merge_hr_proc)
                        if status != Status.COMPLETED:
                            return
                        hr_files[merged_block_id] = merge_hr_proc.outputs.dst_file

                    if MergeableFile.U in self._files_to_merge:
                        # Merging the U (rotation matrix) files
                        merge_u_proc = MergeProcess(merge_function=merge_wannier_u_file_contents,
                                                    src_files=[File(calc, calc.prefix + '_u.mat')
                                                               for calc in src_calcs],
                                                    dst_file=f'{prefix}_u.mat')
                        merge_u_proc.name = f'merge_{merged_block_id.label}_wannier_u'
                        status = self.run_steps(merge_u_proc)
                        if status != Status.COMPLETED:
                            return
                        u_files[merged_block_id] = merge_u_proc.outputs.dst_file

                    if MergeableFile.CENTERS in self._files_to_merge:
                        # Merging the wannier centers files
                        merge_centers_proc = MergeProcess(
                            merge_function=partial(merge_wannier_centers_file_contents, atoms=self.atoms),
                            src_files=[File(calc, calc.prefix + '_centres.xyz') for calc in src_calcs],
                            dst_file=f'{prefix}_centres.xyz')
                        merge_centers_proc.name = f'merge_{merged_block_id.label}_wannier_centers'
                        status = self.run_steps(merge_centers_proc)
                        if status != Status.COMPLETED:
                            return
                        centers_files[merged_block_id] = merge_centers_proc.outputs.dst_file

            if MergeableFile.U_DIS in self._files_to_merge:
                # For the last block (per spin channel), extend the U_dis matrix file if necessary
                spins: List[Spin] = [Spin.UP, Spin.DOWN] if self.parameters.spin_polarized else [Spin.NONE]
                final_label_blocks = [
                    [(block_id, block) for block_id, block
                     in new_projections.to_merge(merge_occ_and_empty=self._merge_occ_and_empty).items()
                     if block_id.spin == s][-1] for s in spins
                ]

                if self._merge_occ_and_empty:
                    raise NotImplementedError("Need to work out what to do here")

                for merged_block_id, block in final_label_blocks:
                    u_dis_file = None
                    num_wann = sum([b.w90_kwargs['num_wann'] for b in block])
                    num_bands = sum([b.w90_kwargs['num_bands'] for b in block])
                    if num_bands > num_wann:
                        calc_with_u_dis = block[-1].w90_calc
                        assert calc_with_u_dis is not None
                        if len(block) == 1:
                            u_dis_file = File(calc_with_u_dis, calc_with_u_dis.prefix + '_u_dis.mat')
                        else:
                            # First, calculate how many empty bands we have
                            spin = merged_block_id.spin
                            if spin != Spin.NONE:
                                nbnd_occ = self.number_of_electrons(spin)
                            else:
                                nbnd_occ = self.number_of_electrons() // 2
                            nbnd_tot = self.calculator_parameters['pw'].nbnd - nbnd_occ

                            # Second, calculate how many empty wannier functions we have
                            nwann_tot = 0
                            for p in block:
                                assert p.num_wann is not None
                                nwann_tot += p.num_wann

                            # Finally, construct and run a Process to perform the file manipulation
                            filling_label = '' if merged_block_id.filled else '_emp'
                            extend_function = partial(extend_wannier_u_dis_file_content, nbnd=nbnd_tot, nwann=nwann_tot)
                            extend_proc = ExtendProcess(extend_function=extend_function,
                                                        src_file=File(calc_with_u_dis,
                                                                      calc_with_u_dis.prefix + '_u_dis.mat'),
                                                        dst_file=calc_with_u_dis.prefix + f'{filling_label}_u_dis.mat')
                            extend_proc.name = f'extend_{merged_block_id.label}_wannier_u_dis'
                            status = self.run_steps(extend_proc)
                            if status != Status.COMPLETED:
                                return
                            u_dis_file = extend_proc.outputs.dst_file
                    u_dis_files[merged_block_id] = u_dis_file

        dos = None
        bs_list = []
        if self.parameters.calculate_bands:
            assert calc_pw_bands is not None

            # Calculate a projected DOS
            if all([p['header']['number_of_wfc'] > 0 for p in self.pseudopotentials.values()]):
                calc_dos = self.new_calculator('projwfc', filpdos=self.name)
                calc_dos.pseudopotentials = self.pseudopotentials
                calc_dos.spin_polarized = self.parameters.spin_polarized
                calc_dos.pseudo_dir = calc_pw_bands.parameters.pseudo_dir
                calc_dos.parameters.prefix = calc_pw_bands.parameters.prefix
                calc_dos.link(File(calc_pw_bands, calc_pw_bands.parameters.outdir), calc_dos.parameters.outdir)
                status = self.run_steps(calc_dos)
                if status != Status.COMPLETED:
                    return

                # Prepare the DOS for plotting
                dos = copy.deepcopy(calc_dos.results['dos'])
            else:
                # Skip if the pseudos don't have the requisite PP_PSWFC blocks
                warn('Some of the pseudopotentials do not have `PP_PSWFC` blocks, which means a projected DOS '
                     'calculation is not possible. Skipping...')

            # Select those calculations that generated a band structure (and are part of this wannierize workflow)
            i_scf = [i for i, c in enumerate(self.calculations) if isinstance(
                c, calculators.PWCalculator) and c.parameters.calculation == 'scf'][-1]
            selected_calcs = [c for c in self.calculations[i_scf:-1]
                              if 'band structure' in c.results and c != calc_pw_bands]

            # Store the pw BandStructure (for working out the vertical shift to set the valence band edge to zero)
            pw_bands = calc_pw_bands.results['band structure']

            # Prepare the band structures for plotting
            labels = ['explicit']
            for c in selected_calcs:
                assert c.parent_process is not None
                assert hasattr(c.parent_process, 'block')
                labels.append(f'interpolation ({c.parent_process.block.id.label.replace("_", " ")})')
            color_cycle = plt.rcParams['axes.prop_cycle']()
            bsplot_kwargs_list = []
            colors: dict[str, str | list[str]] = {}
            for calc, label in zip([calc_pw_bands] + selected_calcs, labels):
                if 'band structure' in calc.results:
                    # Load the bandstructure, shifted by the valence band maximum of the pw bands calculation
                    bs = calc.results['band structure'].subtract_reference(pw_bands.reference)

                    # Tweaking the plot aesthetics
                    plot_kwargs: dict[str, Any] = {'label': label}
                    up_label = label.replace('down', 'up')
                    if 'down' in label:
                        plot_kwargs['ls'] = '--'
                    if 'down' in label and up_label in colors:
                        colors[label] = colors[up_label]
                    else:
                        colors[label] = [next(color_cycle)['color'] for _ in range(bs.energies.shape[0])]
                    if 'explicit' in label:
                        plot_kwargs['ls'] = 'none'
                        plot_kwargs['marker'] = 'x'
                    plot_kwargs['colors'] = colors[label]

                    # Store
                    bs_list.append(bs)
                    bsplot_kwargs_list.append(plot_kwargs)

            # Shift the DOS, too
            if dos is not None:
                dos._energies -= pw_bands.reference

            # Plot
            self.plot_bandstructure(bs_list, dos, bsplot_kwargs=bsplot_kwargs_list)

        # Store the results
        self.outputs = WannierizeOutput(
            band_structures=bs_list, dos=dos, u_files=u_files,
            hr_files=hr_files, centers_files=centers_files, u_dis_files=u_dis_files,
            nnkp_files=nnkp_files,
            nscf_calculation=calc_nscf,
            wannier90_calculations=wannier90_calculations,
            projections=new_projections if new_projections is not None else self.projections)

        self.status = Status.COMPLETED

        return

    def merge_wannier_files(self, block: List[ProjectionsBlock], filling_label: str, prefix: str = 'wann'):
        """Merge various Wannier files for a block of bands.

        This function merges hr (Hamiltonian), u (rotation matrix), and wannier centers files of a collection of blocks
        that share the same filling and spin.
        """
        raise NotImplementedError()


class WannierizeBlockOutput(IOModel):
    """Output model for the WannierizeBlockWorkflow."""

    hr_file: File | None = None
    centers_file: File | None = None
    u_file: File | None = None
    u_dis_file: File | None = None
    amn_file: File | None = None
    eig_file: File | None = None
    mmn_file: File | None = None
    nnkp_file: File | None = None
    chk_file: File
    wannier90_input_file: File
    wannier90_calculation: calculators.Wannier90Calculator


class WannierizeBlockWorkflow(Workflow[WannierizeBlockOutput]):
    """Workflow that Wannierizes a block of bands using Wannier90."""

    output_model = WannierizeBlockOutput

    def __init__(self, *args, block: ProjectionsBlock, force_nspin2=False, minimize=True, pw_outdir: File | None = None,
                 amn_file: File | None = None, eig_file: File | None = None, mmn_file: File | None = None,
                 nnkp_file: File | None = None, **kwargs):
        self.block = block
        self._force_nspin2 = force_nspin2
        self.minimize = minimize
        self.pw_outdir = pw_outdir
        self.amn_file = amn_file
        self.eig_file = eig_file
        self.mmn_file = mmn_file
        self.nnkp_file = nnkp_file
        super().__init__(*args, **kwargs)

    def _run(self) -> None:

        calc_type = 'w90'
        if self.block.spin != Spin.NONE:
            calc_type += f'_{self.block.spin}'

        calc_w90_pp: calculators.Wannier90Calculator | None = None
        if self.amn_file is None or self.eig_file is None or self.mmn_file is None:
            if self.pw_outdir is None:
                raise ValueError(
                    'Wannierization requires either a `pw_outdir` or `.amn`, `.eig`, `.mmn`, and `.nnkp`, files')

            # pre-processing Wannier90 calculation
            calc_w90_pp = self.new_calculator(calc_type, **self.block.w90_kwargs)
            assert isinstance(calc_w90_pp, calculators.Wannier90Calculator)
            calc_w90_pp.prefix = 'wannier90_preproc'
            status = self.run_steps(calc_w90_pp, additional_flags=['-pp'])
            if status != Status.COMPLETED:
                return

            # pw2wannier90 calculation
            calc_p2w: calculators.PW2WannierCalculator = self.new_calculator(
                'pw2wannier', spin_component=self.block.spin)
            calc_p2w.prefix = 'pw2wannier90'
            if calc_w90_pp.parameters.wannier_plot:
                # If we want to plot the Wannier functions, we need to write the UNK files
                calc_p2w.parameters.write_unk = True
            calc_p2w.link(self.pw_outdir, calc_p2w.parameters.outdir, symlink=True)
            calc_p2w.link(File(calc_w90_pp, calc_w90_pp.prefix + '.nnkp'), calc_p2w.parameters.seedname + '.nnkp')
            status = self.run_steps(calc_p2w)
            if status != Status.COMPLETED:
                return

            self.amn_file = File(calc_p2w, calc_p2w.parameters.seedname + '.amn')
            self.eig_file = File(calc_p2w, calc_p2w.parameters.seedname + '.eig')
            self.mmn_file = File(calc_p2w, calc_p2w.parameters.seedname + '.mmn')
            self.nnkp_file = File(calc_p2w, calc_p2w.parameters.seedname + '.nnkp')

        # Wannier90 calculation
        calc_w90: calculators.Wannier90Calculator = self.new_calculator(
            calc_type, bands_plot=self.parameters.calculate_bands, **self.block.w90_kwargs)
        calc_w90.prefix = 'wannier90'
        for f in [self.amn_file, self.eig_file, self.mmn_file, self.nnkp_file]:
            if f is not None:
                calc_w90.link(f, calc_w90.prefix + f.suffix, symlink=True)
        if calc_w90.parameters.wannier_plot:
            for src in File(calc_p2w, '.').glob('UNK*'):
                calc_w90.link(src, src.name, symlink=True)
        status = self.run_steps(calc_w90)
        if status != Status.COMPLETED:
            return
        self.block.w90_calc = calc_w90

        if self.variational_orbitals is not None:
            # Add centers and spreads info to self.variational_orbitals
            exclude_bands = calc_w90.parameters.exclude_bands or []
            orbs = [b for b in self.variational_orbitals if b.spin == self.block.spin
                    and b.index not in exclude_bands]

            centers = calc_w90.results['centers']
            spreads = calc_w90.results['spreads']
            for orb, center, spread in zip(orbs, centers, spreads):
                orb.center = center
                orb.spread = spread

        hr_file = File(calc_w90, calc_w90.prefix + '_hr.dat') if calc_w90.parameters.write_hr else None
        u_file = File(calc_w90, calc_w90.prefix + '_u.mat') if calc_w90.parameters.write_u_matrices else None

        centers_file = File(calc_w90, calc_w90.prefix + '_centres.xyz') if calc_w90.parameters.write_xyz else None
        u_dis_file: File | None = File(calc_w90, calc_w90.prefix + '_u_dis.mat')
        chk_file = File(calc_w90, calc_w90.prefix + '.chk')
        wannier90_input_file = File(calc_w90, calc_w90.prefix + '.win')

        if not u_dis_file.exists():  # type: ignore
            u_dis_file = None
        self.outputs = WannierizeBlockOutput(hr_file=hr_file, u_file=u_file, centers_file=centers_file,
                                             amn_file=self.amn_file, eig_file=self.eig_file, mmn_file=self.mmn_file,
                                             nnkp_file=self.nnkp_file, wannier90_calculation=calc_w90,
                                             u_dis_file=u_dis_file, chk_file=chk_file,
                                             wannier90_input_file=wannier90_input_file)

        self.status = Status.COMPLETED

        return

    def new_calculator(self, *args, **kwargs):
        """Create a new calculator for this workflow."""
        return _internal_new_calculator(self, *args, **kwargs)


class WannierizeAndSplitBlockOutput(IOModel):
    """Output model for the WannierizeAndSplitBlockWorkflow."""

    block_outputs: List[WannierizeBlockOutput]
    blocks: Sequence[ProjectionsBlock]
    u_file: File | None = None
    u_dis_file: File | None = None
    centers_file: File | None = None
    hr_file: File | None = None


class WannierizeAndSplitBlockWorkflow(Workflow[WannierizeAndSplitBlockOutput]):
    """Workflow that Wannierizes a block of bands, splits it using WannierJL, and then Wannierizes the split blocks."""

    output_model = WannierizeAndSplitBlockOutput

    def __init__(self, *args, pw_outdir: File, block: ProjectionsBlock, groups: List[List[int]],
                 force_nspin2: bool = False, minimize: bool = True, files_to_merge: List[MergeableFile] = [], **kwargs):
        self.pw_outdir = pw_outdir
        self.block = block
        self.groups = groups
        self._force_nspin2 = force_nspin2
        self.minimize = minimize
        self._files_to_merge = files_to_merge
        super().__init__(*args, **kwargs)

    def _run(self) -> None:
        # Perform a Block Wannierization
        wannierize_wf = WannierizeBlockWorkflow.fromparent(
            self, pw_outdir=self.pw_outdir, block=self.block, force_nspin2=self._force_nspin2, minimize=self.minimize,
            write_u_matrices=True)
        wannierize_wf.name = f'wannierize-{self.block.id}'
        wannierize_wf.proceed()
        if wannierize_wf.status != Status.COMPLETED:
            return

        # Check that the .nnkp file has the bvectors needed by the parallel transport algorithm
        # For some oblique cells this is not the case by default and we need to construct a cubic .nnkp file
        check_nnkp_process = WannierJLCheckNeighborsProcess(
            wannier90_input_file=wannierize_wf.outputs.wannier90_input_file,
            chk_file=wannierize_wf.outputs.chk_file,
            mmn_file=wannierize_wf.outputs.mmn_file,
            amn_file=wannierize_wf.outputs.amn_file,
            eig_file=wannierize_wf.outputs.eig_file,
        )
        check_nnkp_process.name = 'check_wjl_compatibility'
        status = self.run_steps(check_nnkp_process)
        if status != Status.COMPLETED:
            return

        # Create a cubic nnkp
        if not check_nnkp_process.outputs.has_cubic_neighbors:
            # The .nnkp file is missing some b-vectors that wjl will require to perform the parallel transport, so
            # we must regenerate it
            generate_nnkp_process = WannierJLGenerateNeighborsProcess(
                wannier90_input_file=wannierize_wf.outputs.wannier90_input_file)
            status = self.run_steps(generate_nnkp_process)
            if status != Status.COMPLETED:
                return

            calc_p2w: calculators.PW2WannierCalculator = self.new_calculator(
                'pw2wannier', spin_component=self.block.spin, atom_proj=False)
            calc_p2w.prefix = 'pw2wannier90-cubic'
            calc_p2w.link(self.pw_outdir, calc_p2w.parameters.outdir, symlink=True)
            calc_p2w.link(generate_nnkp_process.outputs.nnkp_file, calc_p2w.parameters.seedname + '.nnkp')
            status = self.run_steps(calc_p2w)
            if status != Status.COMPLETED:
                return
            cubic_nnkp_file = File(calc_p2w, calc_p2w.parameters.seedname + '.nnkp')
            cubic_mmn_file = File(calc_p2w, calc_p2w.parameters.seedname + '.mmn')
        else:
            # The files from the original Wannierization should be sufficient for wjl to perform the parallel transport
            cubic_nnkp_file = None
            cubic_mmn_file = None

        # Determine how to split the blocks
        assert isinstance(self.block, ImplicitProjectionsBlock)
        new_blocks = self.block.split(self.groups)

        # Split the blocks using `wjl`
        split_process = WannierJLSplitProcess(
            indices=self.groups,
            outdirs=[new_block.label for new_block in new_blocks],
            wannier90_input_file=wannierize_wf.outputs.wannier90_input_file,
            chk_file=wannierize_wf.outputs.chk_file,
            mmn_file=wannierize_wf.outputs.mmn_file,
            amn_file=wannierize_wf.outputs.amn_file,
            eig_file=wannierize_wf.outputs.eig_file,
            cubic_nnkp_file=cubic_nnkp_file,
            cubic_mmn_file=cubic_mmn_file
        )

        status = self.run_steps(split_process)
        if status != Status.COMPLETED:
            return

        # Merge the u_dis files
        u_dis_file = None
        if MergeableFile.U_DIS in self._files_to_merge:

            assert wannierize_wf.outputs.u_file is not None
            _, kpoint_list, _ = utils.parse_wannier_u_file_contents(wannierize_wf.outputs.u_file.read_text())

            def combine_u_dis_files(list_of_file_contents: list[str]) -> str:
                # Parse the contents of the U dis files (dimension nkpts x n_bands x n_wann)
                u_dis_matrices = [utils.parse_wannier_amn_file_contents(
                    f, check_square=False) for f in list_of_file_contents]

                # Concatenate the U dis matrices along the n_wann axis
                merged_u_dis_matrix = np.concatenate(u_dis_matrices, axis=2)

                # U_dis files are expected to be of shape (nkpts x n_wann x n_bands)
                merged_u_dis_matrix = np.transpose(merged_u_dis_matrix, (0, 2, 1))

                # Return the merged U dis matrix as a string
                return utils.generate_wannier_u_file_contents(merged_u_dis_matrix, kpoint_list)

            merge_u_dis_process = MergeProcess(
                merge_function=combine_u_dis_files,
                src_files=[b.u_file for b in split_process.outputs.blocks],
                dst_file="wannier90_u_dis.mat"
            )

            merge_u_dis_process.name = 'merge_u_dis'

            status = self.run_steps(merge_u_dis_process)
            if status != Status.COMPLETED:
                return

            u_dis_file = merge_u_dis_process.outputs.dst_file

        # Construct workflows to Wannierize each subblock
        subwfs: List[WannierizeBlockWorkflow] = []
        for new_block, block_files in zip(new_blocks, split_process.outputs.blocks):
            # Wannierize without preprocessing and without disentanglement
            wf = WannierizeBlockWorkflow.fromparent(self,
                                                    block=new_block,
                                                    amn_file=block_files.amn_file,
                                                    eig_file=block_files.eig_file,
                                                    mmn_file=block_files.mmn_file,
                                                    write_u_matrices=True,
                                                    )
            wf.name = f'wannierize-{new_block.id}'
            subwfs.append(wf)

        # Run each subblock
        for wf in subwfs:
            wf.proceed()
        if any([wf.status != Status.COMPLETED for wf in subwfs]):
            return

        # Conbine the U matrices
        u_file = None
        if MergeableFile.U in self._files_to_merge:
            def combine_u_matrices(file_contents: List[str]) -> str:
                u_contents = [utils.parse_wannier_u_file_contents(f) for f in file_contents]

                # Get the k-point information
                _, kpts, n_kpts = u_contents[0]

                # Work out the total number of Wannier functions
                n_wann: int = sum([u[0].shape[1] for u in u_contents])

                # Construct the block-diagonal U matrix
                u_merged = np.zeros((n_kpts, n_wann, n_wann), dtype=np.complex128)
                i_start = 0
                for umat, _, _ in u_contents:
                    u_merged[:, i_start: i_start + umat.shape[1], i_start:i_start + umat.shape[1]] = umat
                    i_start += umat.shape[1]

                return utils.generate_wannier_u_file_contents(u_merged, kpts)

            merge_u_process = MergeProcess(merge_function=combine_u_matrices,
                                           src_files=[subwf.outputs.u_file for subwf in subwfs],
                                           dst_file='wannier90_u.mat')
            merge_u_process.name = 'merge_u'

            status = self.run_steps(merge_u_process)
            if status != Status.COMPLETED:
                return

            u_file = merge_u_process.outputs.dst_file

        # Merge the wannier centers files
        centers_file = None
        if MergeableFile.CENTERS in self._files_to_merge:
            centers_files = [subwf.outputs.centers_file for subwf in subwfs]
            assert isinstance(centers_files[0], File)
            dst_file = centers_files[0].name
            merge_centers_proc = MergeProcess(
                merge_function=partial(merge_wannier_centers_file_contents, atoms=self.atoms),
                src_files=centers_files,
                dst_file=dst_file)
            merge_centers_proc.name = 'merge_centers'
            status = self.run_steps(merge_centers_proc)
            if status != Status.COMPLETED:
                return
            centers_file = merge_centers_proc.outputs.dst_file

        hr_file = None
        if MergeableFile.HR in self._files_to_merge:
            # Merging the wannier_hr (Hamiltonian) files
            hr_files = [subwf.outputs.hr_file for subwf in subwfs]
            assert isinstance(hr_files[0], File)
            dst_file = hr_files[0].name
            merge_hr_proc = MergeProcess(merge_function=merge_wannier_hr_file_contents,
                                         src_files=hr_files,
                                         dst_file=dst_file)
            merge_hr_proc.name = 'merge_hamiltonian'
            status = self.run_steps(merge_hr_proc)
            if status != Status.COMPLETED:
                return
            hr_file = merge_hr_proc.outputs.dst_file

        self.outputs = WannierizeAndSplitBlockOutput(block_outputs=[subwf.outputs for subwf in subwfs],
                                                     blocks=new_blocks,
                                                     u_file=u_file,
                                                     u_dis_file=u_dis_file,
                                                     centers_file=centers_file,
                                                     hr_file=hr_file)

        self.status = Status.COMPLETED

    def new_calculator(self, *args, **kwargs):
        """Create a new calculator for this workflow."""
        return _internal_new_calculator(self, *args, **kwargs)


def _internal_new_calculator(wf, calc_type, *args, **kwargs) -> CalcExtType:  # type: ignore[type-var, misc]
    calc: CalcExtType = super(wf.__class__, wf).new_calculator(calc_type, *args, **kwargs)

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

        if not wf.minimize:
            calc.parameters.num_iter = 0

        if calc.parameters.gamma_only != wf.kpoints.gamma_only:
            # forcing W90 to follow the same logic of PW for the gamma_trick
            calc.parameters.gamma_only = wf.kpoints.gamma_only
    if calc_type == 'pw2wannier':
        if wf._force_nspin2 and not wf.parameters.spin_polarized:
            calc.parameters.spin_component = 'up'

    return calc


def detect_band_blocks(energies: Annotated[npt.NDArray[np.float64], (None, None)],
                       tol: Optional[float] = None,
                       num_occ_bands: Optional[int] = None) -> List[List[int]]:
    """Determine the block of bands in a bandstructure that are separated from one another by at least "tol" eV.

    If "num_occ_bands" is provided, it will also split the bands between occupied and empty bands.
    """
    # If num_occ_bands is not provided, the following lines make sure the subsequent code works
    if num_occ_bands is None:
        num_occ_bands = -1

    # For each band, starting from the second...
    i_bands: List[List[int]] = [[1]]
    for i in range(1, energies.shape[1]):
        if i == num_occ_bands:
            # ... if we have crossed over to the unoccupied bands, create a new group
            i_bands.append([i + 1])
        elif tol and energies[:, i].min() - energies[:, i - 1].max() > tol:
            # ... if it is well-separated, create a new group
            i_bands.append([i + 1])
        else:
            # ... if it is not, add it to the current group
            i_bands[-1].append(i + 1)

    return i_bands
