"""Wannierize workflow module for koopmans.

Workflow module for koopmans, containing the workflow for generating maximally localized
Wannier functions (MLWFs), using the Wannier90 code

Written by Riccardo De Gennaro Nov 2020

"""

import copy
from functools import partial
from typing import Any, Dict, List, Optional, TypeVar

from ase_koopmans.dft.kpoints import BandPath
from ase_koopmans.spectrum.band_structure import BandStructure
from ase_koopmans.spectrum.doscollection import GridDOSCollection

# isort: off
import koopmans.mpl_config  # noqa: F401
import matplotlib.pyplot as plt
# isort: on

from koopmans import calculators, projections, utils
from koopmans.files import File
from koopmans.process_io import IOModel
from koopmans.processes.wannier import (ExtendProcess, MergeProcess,
                                        extend_wannier_u_dis_file_content,
                                        merge_wannier_centers_file_contents,
                                        merge_wannier_hr_file_contents,
                                        merge_wannier_u_file_contents)
from koopmans.projections import BlockID
from koopmans.status import Status
from koopmans.utils import SpinType

from ._workflow import Workflow

CalcExtType = TypeVar('CalcExtType', bound='calculators.CalculatorExt')


class WannierizeOutput(IOModel):
    """Output model for the WannierizeWorkflow."""

    band_structures: List[BandStructure]
    dos: Optional[GridDOSCollection] = None
    u_matrices_files: Dict[BlockID, File | None]
    hr_files: Dict[BlockID, File]
    centers_files: Dict[BlockID, File | None]
    u_dis_files: Dict[BlockID, File | None]
    preprocessing_calculations: List[calculators.Wannier90Calculator]
    nscf_calculation: calculators.PWCalculator
    wannier90_calculations: List[calculators.Wannier90Calculator]


class WannierizeWorkflow(Workflow[WannierizeOutput]):
    """Workflow for Wannierizing an entire system using Quantum ESPRESSO and Wannier90.

    The bands are possibly split into blocks of bands (separated in energy) that are Wannierized separately

    """

    output_model = WannierizeOutput

    def __init__(self, *args, force_nspin2: bool = False, scf_kgrid: List[float] | None = None, **kwargs):
        super().__init__(*args, **kwargs)

        pw_params = self.calculator_parameters['pw']

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs'] \
                and self.parameters.init_empty_orbitals in ['mlwfs', 'projwfs']:

            spins: List[SpinType] = ["up", "down"] if self.parameters.spin_polarized else [None]

            for spin in spins:
                # Work out where we have breaks between blocks of projections, and check that this is commensurate
                # with the number of electrons in this spin channel. Note that this code can be replaced when we have
                # an algorithm for occ-emp separation within Wannier90
                num_bands_occ = self.number_of_electrons(spin)
                if not spin:
                    num_bands_occ //= 2
                divs = self.projections.divisions(spin)
                cumulative_divs = [sum(divs[:i + 1]) for i in range(len(divs))]
                if num_bands_occ not in cumulative_divs:
                    message = 'The provided `Wannier90` projections are not commensurate with the number of ' \
                              'electrons; divide your list of projections into sublists that represent blocks ' \
                              'of bands to Wannierize separately'
                    raise ValueError(message)

                # Compare the number of bands from PW to Wannier90
                num_bands_w90 = self.projections.num_bands(spin=spin)
                if num_bands_w90 > pw_params.nbnd:
                    raise ValueError(f'You have provided more bands to the `Wannier90 calculator` ({num_bands_w90}) '
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
            raise NotImplementedError('`WannierizeWorkflow` only supports setting `init_orbitals` and '
                                      '`init_empty_orbitals` to `mlwfs`/`projwfs`/`kohn-sham`')

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

    def _run(self) -> None:
        """Run thw workflow."""
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

        u_matrices_files = {}
        hr_files: Dict[BlockID, File] = {}
        centers_files = {}
        u_dis_files = {}
        wannier90_calculations = []
        preprocessing_calculations = []

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs'] \
                and self.parameters.init_empty_orbitals in ['mlwfs', 'projwfs']:

            block_subworkflows: list[Workflow] = []
            for block in self.projections:
                wannierize_block_subworkflow = WannierizeBlockWorkflow.fromparent(
                    self, force_nspin2=self._force_nspin2, block=block,
                    pw_outdir=File(calc_nscf, calc_nscf.parameters.outdir))
                wannierize_block_subworkflow.name = \
                    f'Wannierize {block.id.label.replace("_", " ").replace("block", "Block")}'
                block_subworkflows.append(wannierize_block_subworkflow)

            for wf in block_subworkflows:
                wf.run()
            if any([wf.status != Status.COMPLETED for wf in block_subworkflows]):
                return

            for block, subwf in zip(self.projections, block_subworkflows):
                assert isinstance(subwf, WannierizeBlockWorkflow)

                # Store the results
                assert subwf.outputs.hr_file is not None
                hr_files[block.id] = subwf.outputs.hr_file
                centers_files[block.id] = subwf.outputs.centers_file
                u_matrices_files[block.id] = subwf.outputs.u_matrices_file
                wannier90_calculations.append(subwf.outputs.wannier90_calculation)
                preprocessing_calculations.append(subwf.outputs.preprocessing_calculation)

            # Merging Hamiltonian files, U matrix files, centers files if necessary
            if self.parent_process is not None:

                for block_id, block in self.projections.to_merge.items():
                    if len(block) == 1:
                        # If there is only one block, we don't need to merge anything
                        calc = block[0].w90_calc
                        assert calc is not None
                        if calc.parameters.write_hr:
                            hr_files[block_id] = File(calc, calc.prefix + '_hr.dat')
                        if calc.parameters.write_u_matrices:
                            u_matrices_files[block_id] = File(calc, calc.prefix + '_u.mat')
                        if calc.parameters.write_xyz:
                            centers_files[block_id] = File(calc, calc.prefix + '_centres.xyz')
                    else:
                        # Fetching the list of calculations for this block
                        src_calcs: List[calculators.Wannier90Calculator] = [
                            b.w90_calc for b in block if b.w90_calc is not None]
                        prefix = src_calcs[-1].prefix

                        # Merging the wannier_hr (Hamiltonian) files
                        merge_hr_proc = MergeProcess(merge_function=merge_wannier_hr_file_contents,
                                                     src_files=[File(calc, calc.prefix + '_hr.dat')
                                                                for calc in src_calcs],
                                                     dst_file=f'{prefix}_hr.dat')
                        merge_hr_proc.name = f'merge_{block_id.label}_wannier_hamiltonian'
                        status = self.run_steps(merge_hr_proc)
                        if status != Status.COMPLETED:
                            return
                        hr_files[block_id] = merge_hr_proc.outputs.dst_file

                        if self.parameters.method == 'dfpt' and self.parent_process is not None:
                            # Merging the U (rotation matrix) files
                            merge_u_proc = MergeProcess(merge_function=merge_wannier_u_file_contents,
                                                        src_files=[File(calc, calc.prefix + '_u.mat')
                                                                   for calc in src_calcs],
                                                        dst_file=f'{prefix}_u.mat')
                            merge_u_proc.name = f'merge_{block_id.label}_wannier_u'
                            status = self.run_steps(merge_u_proc)
                            if status != Status.COMPLETED:
                                return
                            u_matrices_files[block_id] = merge_u_proc.outputs.dst_file

                            # Merging the wannier centers files
                            merge_centers_proc = MergeProcess(
                                merge_function=partial(merge_wannier_centers_file_contents, atoms=self.atoms),
                                src_files=[File(calc, calc.prefix + '_centres.xyz') for calc in src_calcs],
                                dst_file=f'{prefix}_centres.xyz')
                            merge_centers_proc.name = f'merge_{block_id.label}_wannier_centers'
                            status = self.run_steps(merge_centers_proc)
                            if status != Status.COMPLETED:
                                return
                            centers_files[block_id] = merge_centers_proc.outputs.dst_file

                # For the last block (per spin channel), extend the U_dis matrix file if necessary
                spins: List[SpinType] = ['up', 'down'] if self.parameters.spin_polarized else [None]
                final_label_blocks = [
                    [(block_id, block) for block_id, block in self.projections.to_merge.items()
                     if block_id.spin == s][-1] for s in spins]

                for block_id, block in final_label_blocks:
                    u_dis_file = None
                    num_wann = sum([b.w90_kwargs['num_wann'] for b in block])
                    num_bands = sum([b.w90_kwargs['num_bands'] for b in block])
                    if num_bands > num_wann and self.parameters.method == 'dfpt':
                        calc_with_u_dis = block[-1].w90_calc
                        assert calc_with_u_dis is not None
                        if len(block) == 1:
                            u_dis_file = File(calc_with_u_dis, calc_with_u_dis.prefix + '_u_dis.mat')
                        else:
                            # First, calculate how many empty bands we have
                            spin = block_id.spin
                            if spin:
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
                            filling_label = '' if block_id.filled else '_emp'
                            extend_function = partial(extend_wannier_u_dis_file_content, nbnd=nbnd_tot, nwann=nwann_tot)
                            extend_proc = ExtendProcess(extend_function=extend_function,
                                                        src_file=(calc_with_u_dis,
                                                                  calc_with_u_dis.prefix + '_u_dis.mat'),
                                                        dst_file=calc_with_u_dis.prefix + f'{filling_label}_u_dis.mat')
                            extend_proc.name = f'extend_{block_id.label}_wannier_u_dis'
                            status = self.run_steps(extend_proc)
                            if status != Status.COMPLETED:
                                return
                            u_dis_file = extend_proc.outputs.dst_file
                    u_dis_files[block_id] = u_dis_file

        dos = None
        bs_list = []
        if self.parameters.calculate_bands:
            # Run a "bands" calculation, making sure we don't overwrite
            # the scf/nscf tmp files by setting a different prefix
            calc_pw_bands = self.new_calculator('pw', calculation='bands', kpts=self.kpoints.path)
            calc_pw_bands.prefix = 'bands'
            calc_pw_bands.parameters.prefix += '_bands'

            # Link the save directory so that the bands calculation can use the old density
            src = File(calc_nscf, (calc_nscf.parameters.outdir / calc_nscf.parameters.prefix).with_suffix('.save'))
            dst = (calc_pw_bands.parameters.outdir / calc_pw_bands.parameters.prefix).with_suffix('.save')
            calc_pw_bands.link(src, dst)

            # Run the bands calculation
            status = self.run_steps(calc_pw_bands)
            if status != Status.COMPLETED:
                return

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
                utils.warn('Some of the pseudopotentials do not have `PP_PSWFC` blocks, which means a projected DOS '
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
                    kwargs: dict[str, Any] = {'label': label}
                    up_label = label.replace('down', 'up')
                    if 'down' in label:
                        kwargs['ls'] = '--'
                    if 'down' in label and up_label in colors:
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

        # Store the results
        self.outputs = self.output_model(band_structures=bs_list, dos=dos, u_matrices_files=u_matrices_files,
                                         hr_files=hr_files, centers_files=centers_files, u_dis_files=u_dis_files,
                                         preprocessing_calculations=preprocessing_calculations,
                                         nscf_calculation=calc_nscf,
                                         wannier90_calculations=wannier90_calculations)

        self.status = Status.COMPLETED

        return

    def merge_wannier_files(self, block: List[projections.ProjectionBlock], filling_label: str, prefix: str = 'wann'):
        """Merge various Wannier files for a block of bands.

        This function merges hr (Hamiltonian), u (rotation matrix), and wannier centers files of a collection of blocks
        that share the same filling and spin.
        """
        raise NotImplementedError()


class WannierizeBlockOutput(IOModel):
    """Output model for the WannierizeBlockWorkflow."""

    hr_file: File | None = None
    centers_file: File | None = None
    u_matrices_file: File | None = None
    wannier90_calculation: calculators.Wannier90Calculator
    preprocessing_calculation: calculators.Wannier90Calculator


class WannierizeBlockWorkflow(Workflow[WannierizeBlockOutput]):
    """Workflow that Wannierizes a block of bands using Wannier90."""

    output_model = WannierizeBlockOutput

    def __init__(self, *args, pw_outdir: File, block: projections.ProjectionBlock, force_nspin2=False, **kwargs):
        self.pw_outdir = pw_outdir
        self.block = block
        self._force_nspin2 = force_nspin2
        super().__init__(*args, **kwargs)

    def _run(self) -> None:
        n_occ_bands = self.number_of_electrons(self.block.spin)
        if not self.block.spin:
            n_occ_bands //= 2

        assert self.block.include_bands is not None
        if max(self.block.include_bands) <= n_occ_bands:
            # Block consists purely of occupied bands
            init_orbs = self.parameters.init_orbitals
        else:
            init_orbs = self.parameters.init_empty_orbitals
        # Store the number of electrons in the ProjectionBlocks object so that it can work out which blocks to
        # merge with one another
        self.projections.num_occ_bands[self.block.spin] = n_occ_bands

        calc_type = 'w90'
        if self.block.spin:
            calc_type += f'_{self.block.spin}'

        # 1) pre-processing Wannier90 calculation
        calc_w90_pp: calculators.Wannier90Calculator = self.new_calculator(
            calc_type, init_orbitals=init_orbs, **self.block.w90_kwargs)
        calc_w90_pp.prefix = 'wannier90_preproc'
        status = self.run_steps(calc_w90_pp, additional_flags=['-pp'])
        if status != Status.COMPLETED:
            return

        # 2) standard pw2wannier90 calculation
        calc_p2w: calculators.PW2WannierCalculator = self.new_calculator('pw2wannier', spin_component=self.block.spin)
        calc_p2w.prefix = 'pw2wannier90'
        if calc_w90_pp.parameters.wannier_plot:
            # If we want to plot the Wannier functions, we need to write the UNK files
            calc_p2w.parameters.write_unk = True
        calc_p2w.link(self.pw_outdir, calc_p2w.parameters.outdir, symlink=True)
        calc_p2w.link(File(calc_w90_pp, calc_w90_pp.prefix + '.nnkp'), calc_p2w.parameters.seedname + '.nnkp')
        status = self.run_steps(calc_p2w)
        if status != Status.COMPLETED:
            return

        # 3) Wannier90 calculation
        calc_w90: calculators.Wannier90Calculator = \
            self.new_calculator(calc_type, init_orbitals=init_orbs,
                                bands_plot=self.parameters.calculate_bands, **self.block.w90_kwargs)
        calc_w90.prefix = 'wannier90'

        # Link the pw2wannier90 output files
        for ext in ['.eig', '.amn', '.mmn']:
            calc_w90.link(File(calc_p2w, calc_p2w.parameters.seedname + ext), calc_w90.prefix + ext, symlink=True)
        if calc_w90.parameters.wannier_plot:
            for src in File(calc_p2w, '.').glob('UNK*'):
                calc_w90.link(src, src.name, symlink=True)

        status = self.run_steps(calc_w90)
        if status != Status.COMPLETED:
            return
        self.block.w90_calc = calc_w90

        if self.bands is not None:
            # Add centers and spreads info to self.bands
            if self.block.spin is None:
                remaining_bands = [b for b in self.bands if b.center is None and b.spin == 0]
            else:
                if self.block.spin == 'up':
                    i_spin = 0
                else:
                    i_spin = 1
                remaining_bands = [b for b in self.bands if b.center is None and b.spin == i_spin]

            centers = calc_w90.results['centers']
            spreads = calc_w90.results['spreads']
            for band, center, spread in zip(remaining_bands, centers, spreads):
                band.center = center
                band.spread = spread

                if self.block.spin is None and len(self.bands.get(spin=1)) > 0:
                    # Copy over spin-up results to spin-down
                    [match] = [b for b in self.bands if b.index == band.index and b.spin == 1]
                    match.center = center
                    match.spread = spread

        hr_file = File(calc_w90, calc_w90.prefix + '_hr.dat') if calc_w90.parameters.write_hr else None
        u_file = File(calc_w90, calc_w90.prefix + '_u.mat') if calc_w90.parameters.write_u_matrices else None
        centers_file = File(calc_w90, calc_w90.prefix + '_centres.xyz') if calc_w90.parameters.write_xyz else None
        self.outputs = self.output_model(hr_file=hr_file, u_matrices_file=u_file, centers_file=centers_file,
                                         preprocessing_calculation=calc_w90_pp, wannier90_calculation=calc_w90)

        self.status = Status.COMPLETED

        return

    def new_calculator(self, calc_type, *args, **kwargs) -> CalcExtType:  # type: ignore[type-var, misc]
        """Create a new calculator, with some extra tweaks for Wannier90 calculations."""
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
                raise ValueError(f'Unrecognized orbital type `{init_orbs}` (must be `mlwfs`/`projwfs`)')

            if calc.parameters.gamma_only != self.kpoints.gamma_only:
                # forcing W90 to follow the same logic of PW for the gamma_trick
                calc.parameters.gamma_only = self.kpoints.gamma_only
        if calc_type == 'pw2wannier':
            if self._force_nspin2 and not self.parameters.spin_polarized:
                calc.parameters.spin_component = 'up'

        return calc
