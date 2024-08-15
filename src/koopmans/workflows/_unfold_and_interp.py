"""

Workflow module for koopmans, containing the workflow for unfolding and interpolating
results from the supercell to generate a bandstructure

For the moment the code works only with cubic, tetragonal and orthorhombic systems.

Originally written by Riccardo De Gennaro as the standalone 'unfolding and interpolate' code
Integrated within koopmans by Edward Linscott Jan 2021
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from ase.dft.dos import DOS
from ase.spectrum.band_structure import BandStructure

from koopmans import calculators, outputs, processes, utils
from koopmans.files import FilePointer

from ._workflow import Workflow


class UnfoldAndInterpolateOutput(outputs.OutputModel):
    band_structure: BandStructure
    dos: DOS | None

    class Config:
        arbitrary_types_allowed = True


class UnfoldAndInterpolateWorkflow(Workflow):

    output_model = UnfoldAndInterpolateOutput  # type: ignore

    def __init__(self, *args, koopmans_ham_files: Dict[Tuple[str, str | None], FilePointer],
                 dft_ham_files: Dict[Tuple[str, str | None], FilePointer],
                 redo_smooth_dft: Optional[bool] = None, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._dft_ham_files = dft_ham_files
        self._koopmans_ham_files = koopmans_ham_files
        self._redo_smooth_dft = redo_smooth_dft

    def _run(self) -> None:
        '''

        Wrapper for the whole unfolding and interpolation workflow, consisting of:
        - a smooth WannierizeWorkflow (if required)
        - an UnfoldAndInterpolateProcess for occ states
        - an UnfoldAndInterpolateProcess for emp states
        - an UnfoldAndInterpolateProcess to merge occ and emp results

        '''
        # Import these here so that if these have been monkey-patched, we get the monkey-patched version
        from koopmans.workflows import WannierizeWorkflow

        # Transform self.atoms back to the primitive cell
        if self.parameters.method == 'dscf':
            self.supercell_to_primitive()

        # Store the original w90 calculations
        w90_calcs = [c for c in self.calculations if isinstance(c, calculators.Wannier90Calculator)
                     and c.command.flags == ''][-len(self.projections):]

        if self.calculator_parameters['ui'].do_smooth_interpolation:
            wannier_workflow = WannierizeWorkflow.fromparent(self, scf_kgrid=self.kpoints.grid)
            assert self.kpoints.grid is not None
            wannier_workflow.kpoints.grid = [x * y for x, y in zip(self.kpoints.grid,
                                             self.calculator_parameters['ui'].smooth_int_factor)]

            # Here, we allow for skipping of the smooth dft calcs (assuming they have been already run)
            # This is achieved via the optional argument of from_scratch in run(), which
            # overrides the value of wannier_workflow.from_scratch, as well as preventing the inheritance of
            # self.from_scratch to wannier_workflow.from_scratch and back again after the subworkflow finishes
            wannier_workflow.run(from_scratch=self._redo_smooth_dft)

        process: processes.UnfoldAndInterpolateProcess
        spins: List[Optional[str]]
        if self.parameters.spin_polarized:
            spins = ['up', 'down']
        else:
            spins = [None]

        for spin, band_filling in zip(spins, self.bands.filling):
            # Extract the centers and spreads corresponding to this particular spin
            centers = np.array([center for c, p in zip(w90_calcs, self.projections)
                               for center in c.results['centers'] if p.spin == spin])
            spreads = np.array([spread for c, p in zip(w90_calcs, self.projections)
                               for spread in c.results['spreads'] if p.spin == spin])

            for filled, filling in zip([True, False], ['occ', 'emp']):
                presets = filling
                if spin:
                    presets += '_' + spin

                # Extract the centers and spreads that have this particular filling
                if self.parameters.method == 'dscf':
                    # For dscf, self.bands correspond to the supercell so band_filling involves many copies of each
                    # band
                    assert self.kpoints.grid is not None
                    ngrid = np.prod(self.kpoints.grid, dtype=int)
                else:
                    # For dfpt, self.bands correspond to the primitive cell so band_filling is already the correct
                    # dimensions
                    ngrid = 1
                mask = np.array(band_filling[::ngrid]) == filled

                # Add the smooth DFT Hamiltonian file if relevant
                if self.calculator_parameters['ui'].do_smooth_interpolation:
                    dft_smooth_ham_file = wannier_workflow.outputs.hr_files[presets]
                else:
                    dft_smooth_ham_file = None

                process = self.new_ui_process(presets, centers=centers[mask], spreads=spreads[mask].tolist(),
                                              dft_smooth_ham_file=dft_smooth_ham_file)

                # Run the process
                self.run_process(process)

        # Merge the bands
        if self.parameters.spin_polarized:
            energies = [[p.outputs.band_structure.energies for p in subset]
                        for subset in [self.processes[-4:-2], self.processes[-2:]]]
            reference = np.max([e[0] for e in energies])
            energies_np = np.concatenate([np.concatenate(e, axis=2) for e in energies], axis=0)
        else:
            energies = [p.outputs.band_structure.energies for p in self.processes[-2:]]
            reference = np.max(energies[0])
            energies_np = np.concatenate(energies, axis=2)
        merged_bs = BandStructure(self.kpoints.path, energies_np, reference=reference)

        if process.inputs.parameters.do_dos:
            merged_dos = processes.generate_dos(merged_bs, self.plotting)

        # Plot the band structure and DOS
        if process.inputs.parameters.do_dos:
            # Add the DOS only if the k-path is sufficiently sampled to mean the individual Gaussians are not visible
            # (by comparing the median jump between adjacent eigenvalues to the smearing width)
            median_eval_gap = max([np.median(e[1:] - e[:-1])
                                  for e in [np.sort(ekn.flatten()) for ekn in merged_dos.e_skn]])
            if merged_dos.width < 5 * median_eval_gap:
                merged_dos = None
                utils.warn('The DOS will not be plotted, because the Brillouin zone is too poorly sampled for the '
                           'specified value of smearing. In order to generate a DOS, increase the k-point density '
                           '(`kpath_density` in the `setup` `k_points` subblock) and/or the smearing (`degauss` '
                           'in the `plot` block)')
        else:
            merged_dos = None

        # Shift the DOS to align with the band structure
        if merged_dos is not None:
            merged_dos.e_skn -= merged_bs.reference

        self.plot_bandstructure(merged_bs.subtract_reference(), merged_dos)

        # Shift the DOS back
        if merged_dos is not None:
            merged_dos.e_skn += merged_bs.reference

        # Store the results
        self.outputs = self.output_model(band_structure=merged_bs, dos=merged_dos)

    def new_ui_process(self, presets: str, **kwargs) -> processes.UnfoldAndInterpolateProcess:
        valid_presets = ['occ', 'occ_up', 'occ_down', 'emp', 'emp_up', 'emp_down']
        assert presets in valid_presets, \
            'In UnfoldAndInterpolateWorkflow.new_ui_process() presets must be ' \
            + '/'.join([f'"{s}"' for s in valid_presets]) + \
            f', but you have tried to set it equal to {presets}'

        preset_tuple: Tuple[str, str | None]
        if '_' in presets:
            preset_tuple = tuple(presets.split('_'))  # type: ignore
        else:
            preset_tuple = (presets, None)

        kwargs['kc_ham_file'] = self._koopmans_ham_files[preset_tuple]
        if self.calculator_parameters['ui'].do_smooth_interpolation:
            kwargs['dft_ham_file'] = self._dft_ham_files[preset_tuple]

        # # Automatically generating UI calculator settings
        # if self.parameters.method == 'dscf':
        #     # DSCF case
        #     import ipdb; ipdb.set_trace()
        #     if '_' in presets:
        #         ham_prefix = presets.replace('up', '1').replace('down', '2')
        #     else:
        #         ham_prefix = presets + '_1'
        # else:
        #     # DFPT case
        #     if self.parameters.spin_polarized:
        #         raise NotImplementedError()

        #     # Add the Koopmans Hamiltonian
        #     kc_ham_calc = [c for c in self.calculations if isinstance(c, calculators.KoopmansHamCalculator)][-1]
        #     ham_filename = f'{kc_ham_calc.parameters.prefix}.kcw_hr_{presets}.dat'
        #     kwargs['kc_ham_file'] = (kc_ham_calc, ham_filename)

        #     # Add the DFT Hamiltonian files if performing smooth interpolation
        #     if self.calculator_parameters['ui'].do_smooth_interpolation:
        #         suffix = '' if presets == 'occ' else '_emp'
        #         ham_filename = f'wannier90{suffix}_hr.dat'
        #         [wann_calc, smooth_wann_calc] = [p for p in self.processes if str(
        #             getattr(p.inputs, 'dst_file', '')) == ham_filename][-2:]
        #         kwargs['dft_ham_file'] = (wann_calc, ham_filename)
        #         kwargs['dft_smooth_ham_file'] = (smooth_wann_calc, ham_filename)

        parameters = self.calculator_parameters['ui']
        parameters.kgrid = self.kpoints.grid
        parameters.kpath = self.kpoints.path

        process = processes.UnfoldAndInterpolateProcess(atoms=self.atoms,
                                                        parameters=parameters,
                                                        plotting_parameters=self.plotting,
                                                        **kwargs)

        process.name += f'_{presets}'

        return process
