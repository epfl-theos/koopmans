"""

Workflow module for koopmans, containing the workflow for unfolding and interpolating
results from the supercell to generate a bandstructure

For the moment the code works only with cubic, tetragonal and orthorhombic systems.

Originally written by Riccardo De Gennaro as the standalone 'unfolding and interpolate' code
Integrated within koopmans by Edward Linscott Jan 2021
"""

from typing import Dict, Generator, List, Literal, Optional, Tuple
from pydantic import ConfigDict

import numpy as np
from ase_koopmans.dft.dos import DOS
from ase_koopmans.spectrum.band_structure import BandStructure

from koopmans import calculators, utils
from koopmans.process_io import IOModel
from koopmans.files import File
from koopmans.processes.ui import UnfoldAndInterpolateProcess, generate_dos
from koopmans.projections import BlockID
from koopmans.status import Status

from ._wannierize import WannierizeWorkflow
from ._workflow import Workflow


class UnfoldAndInterpolateOutput(IOModel):
    band_structure: BandStructure
    dos: Optional[DOS]
    smooth_dft_ham_files: Optional[Dict[BlockID, File]]
    model_config = ConfigDict(arbitrary_types_allowed=True)


class UnfoldAndInterpolateWorkflow(Workflow):

    output_model = UnfoldAndInterpolateOutput  # type: ignore

    def __init__(self, *args, koopmans_ham_files: Dict[BlockID, File],
                 dft_ham_files: Dict[BlockID, File],
                 smooth_dft_ham_files: Dict[BlockID, File] | None = None, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._dft_ham_files = dft_ham_files
        self._koopmans_ham_files = koopmans_ham_files
        self._smooth_dft_ham_files = smooth_dft_ham_files

    def _run(self) -> None:
        '''

        Wrapper for the whole unfolding and interpolation workflow, consisting of:
        - a smooth WannierizeWorkflow (if required)
        - an UnfoldAndInterpolateProcess for occ states
        - an UnfoldAndInterpolateProcess for emp states
        - an UnfoldAndInterpolateProcess to merge occ and emp results

        '''

        # Transform self.atoms back to the primitive cell
        if self.parameters.method == 'dscf':
            self.supercell_to_primitive()

        # Store the original w90 calculations
        w90_calcs = [c for c in self.calculations if isinstance(c, calculators.Wannier90Calculator)
                     and c.command.flags == ''][-len(self.projections):]

        if self.calculator_parameters['ui'].do_smooth_interpolation and self._smooth_dft_ham_files is None:

            assert self.kpoints.grid is not None
            wannier_workflow = WannierizeWorkflow.fromparent(self, scf_kgrid=self.kpoints.grid)
            wannier_workflow.kpoints.grid = [x * y for x, y in zip(self.kpoints.grid,
                                             self.calculator_parameters['ui'].smooth_int_factor)]

            wannier_workflow.run()
            if wannier_workflow.status != Status.COMPLETED:
                return

            # Save the smooth DFT Hamiltonian files
            self._smooth_dft_ham_files = wannier_workflow.outputs.hr_files

        process: UnfoldAndInterpolateProcess
        spins: List[Literal[None, "up", "down", "spinor"]]
        if self.parameters.spin_polarized:
            spins = ['up', 'down']
        else:
            spins = [None]

        assert self.bands is not None
        processes = []
        for spin, band_filling in zip(spins, self.bands.filling):
            # Extract the centers and spreads corresponding to this particular spin
            centers = np.array([center for c, p in zip(w90_calcs, self.projections)
                               for center in c.results['centers'] if p.spin == spin])
            spreads = np.array([spread for c, p in zip(w90_calcs, self.projections)
                               for spread in c.results['spreads'] if p.spin == spin])

            for filled in [True, False]:
                block_id = BlockID(filled=filled, spin=spin)

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
                    assert self._smooth_dft_ham_files is not None
                    dft_smooth_ham_file = self._smooth_dft_ham_files[block_id]
                else:
                    dft_smooth_ham_file = None

                process = self.new_ui_process(block_id, centers=centers[mask], spreads=spreads[mask].tolist(),
                                              dft_smooth_ham_file=dft_smooth_ham_file)

                processes.append(process)

        # Run the processes
        status = self.run_steps(processes)
        if status != Status.COMPLETED:
            return

        # Merge the bands
        if self.parameters.spin_polarized:
            energies = [[p.outputs.band_structure.energies for p in subset]
                        for subset in [self.processes[-4:-2], self.processes[-2:]]]
            energies_np = np.concatenate([np.concatenate(e, axis=2) for e in energies], axis=0)
            reference = np.max([np.max(e[0]) for e in energies])
        else:
            energies = [p.outputs.band_structure.energies for p in self.processes[-2:]]
            reference = np.max(energies[0])
            energies_np = np.concatenate(energies, axis=2)
        merged_bs = BandStructure(self.kpoints.path, energies_np, reference=reference)

        if process.inputs.parameters.do_dos:
            merged_dos = generate_dos(merged_bs, self.plotting, self.parameters.spin_polarized)

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
        self.outputs = self.output_model(band_structure=merged_bs, dos=merged_dos,
                                         smooth_dft_ham_files=self._smooth_dft_ham_files)

        self.status = Status.COMPLETED

        return

    def new_ui_process(self, block_id: BlockID, **kwargs) -> UnfoldAndInterpolateProcess:
        kwargs['kc_ham_file'] = self._koopmans_ham_files[block_id]
        kwargs['dft_ham_file'] = self._dft_ham_files[block_id]

        parameters = self.calculator_parameters['ui']
        parameters.kgrid = self.kpoints.grid
        parameters.kpath = self.kpoints.path

        process = UnfoldAndInterpolateProcess(atoms=self.atoms,
                                              parameters=parameters,
                                              plotting_parameters=self.plotting,
                                              **kwargs)

        process.name += f'_{block_id.label}'

        return process
