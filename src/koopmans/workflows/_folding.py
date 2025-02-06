"""

Workflow module for koopmans, containing the workflow for converting W90 or PW files to
kcp.x friendly-format using wann2kcp.x

Written by Edward Linscott Feb 2021

"""

import os
from pathlib import Path
from typing import Dict, Generator, List, Tuple

import numpy as np

from koopmans import calculators, utils
from koopmans.files import FilePointer
from koopmans.outputs import OutputModel
from koopmans.processes.merge_evc import MergeEVCProcess
from koopmans.projections import BlockID, ProjectionBlock
from koopmans.status import Status
from koopmans.step import Step

from ._workflow import Workflow


class FoldToSupercellOutputs(OutputModel):
    kcp_files: Dict[str, FilePointer]

    class Config:
        arbitrary_types_allowed = True


class FoldToSupercellWorkflow(Workflow):

    output_model = FoldToSupercellOutputs  # type: ignore

    def __init__(self, nscf_outdir: FilePointer, hr_files: Dict[str, FilePointer], wannier90_calculations, wannier90_pp_calculations, **kwargs):
        super().__init__(**kwargs)
        self._nscf_outdir = nscf_outdir
        self._hr_files = hr_files
        self._wannier90_calculations = wannier90_calculations
        self._wannier90_pp_calculations = wannier90_pp_calculations

    def _run(self) -> None:
        '''

        Wrapper for folding Wannier or Kohn-Sham functions from the primitive cell
        to the supercell and convert them to a kcp.x friendly format.

        '''

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs']:
            # Loop over the various subblocks that we have wannierized separately
            converted_files = {}
            w2k_calcs = []
            for w90_calc, w90_pp_calc, block in zip(self._wannier90_calculations, self._wannier90_pp_calculations, self.projections):
                # Create the calculator
                calc_w2k = self.new_calculator('wann2kcp', spin_component=block.spin, wan_mode='wannier2kcp')
                calc_w2k.prefix = f'convert_{block.name}_to_supercell'

                # Checking that gamma_trick is consistent with gamma_only
                if calc_w2k.parameters.gamma_trick and not self.kpoints.gamma_only:
                    calc_w2k.parameters.gamma_trick = False
                elif not calc_w2k.parameters.gamma_trick and self.kpoints.gamma_only:
                    calc_w2k.parameters.gamma_trick = True
                else:
                    pass

                # Link the input files
                self.link(*self._hr_files[block.id], calc_w2k, self._hr_files[block.id].name, symlink=True)
                self.link(*self._nscf_outdir, calc_w2k, calc_w2k.parameters.outdir, recursive_symlink=True)
                self.link(w90_pp_calc, w90_pp_calc.prefix + '.nnkp', calc_w2k,
                          calc_w2k.parameters.seedname + '.nnkp', symlink=True)
                self.link(w90_calc, w90_calc.prefix + '.chk', calc_w2k,
                          calc_w2k.parameters.seedname + '.chk', symlink=True)

                w2k_calcs.append(calc_w2k)

            # Run the calculators (possibly in parallel)
            status = self.run_steps(w2k_calcs)
            if status != Status.COMPLETED:
                return

            for block, calc_w2k in zip(self.projections, w2k_calcs):

                if self.parameters.spin_polarized:
                    converted_files[block.id] = [FilePointer(calc_w2k, Path("evcw.dat"))]
                else:
                    converted_files[block.id] = [FilePointer(calc_w2k, Path("evcw1.dat")),
                                                 FilePointer(calc_w2k, Path("evcw2.dat"))]

            # Merging evcw files
            merged_files = {}
            for merged_id, subset in self.projections.to_merge.items():
                if len(subset) == 1:
                    for f in converted_files[subset[0].id]:
                        dest_file = _construct_dest_filename(f.name, merged_id)
                        merged_files[dest_file] = f
                else:
                    if self.parameters.spin_polarized:
                        evc_fnames = ['evcw.dat']
                    else:
                        evc_fnames = ['evcw1.dat', 'evcw2.dat']

                    for evc_fname in evc_fnames:
                        src_files = [f for s in subset for f in converted_files[s.id] if str(f.name) == evc_fname]
                        merge_proc = MergeEVCProcess(kgrid=self.kpoints.grid,
                                                     src_files=src_files, dest_filename=evc_fname)
                        if merged_id.filled:
                            tidy_label = 'occupied'
                        else:
                            tidy_label = 'empty'
                        if merged_id.spin != None:
                            tidy_label += f'_spin_{subset[0].spin}'
                        merge_proc.name = 'merge_wavefunctions_for_' + tidy_label
                        status = self.run_steps(merge_proc)
                        if status != Status.COMPLETED:
                            return

                        dest_file = _construct_dest_filename(evc_fname, merged_id)
                        merged_files[dest_file] = merge_proc.outputs.merged_file

        else:
            # Create the calculator
            calc_w2k = self.new_calculator('wann2kcp', directory='ks2kcp', wan_mode='ks2kcp', seedname=None)
            calc_w2k.prefix = 'ks2kcp'

            # Run the calculator
            status = self.run_steps(calc_w2k)
            if status != Status.COMPLETED:
                return

            raise NotImplementedError("Need to populate converted and merged file dictionaries")

        self.outputs = self.output_model(kcp_files=merged_files)

        self.status = Status.COMPLETED

        return


def _construct_dest_filename(fname: str | Path, merged_id: BlockID) -> str:
    fname = str(fname)
    if fname[4] in ['1', '2']:
        spin = int(fname[4])
    elif merged_id.spin == 'up':
        spin = 1
    elif merged_id.spin == 'down':
        spin = 2
    else:
        raise ValueError('Should not arrive here')
    if merged_id.filled:
        return f'evc_occupied{spin}.dat'
    else:
        return f'evc0_empty{spin}.dat'
