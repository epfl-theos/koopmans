"""Processes for modifying Wannier90 files."""

import math
from pathlib import Path
from typing import Callable, List

import numpy as np
from ase_koopmans import Atoms
from pydantic import ConfigDict

from koopmans import utils
from koopmans.files import File

from ._process import IOModel, Process


def merge_wannier_hr_file_contents(filecontents: List[str]) -> str:
    """Merge the contents of multiple `Wannier90` hr files into a single file."""
    # Reading in each hr file in turn
    hr_list = []
    weights_out = None
    rvect_out = None
    for filecontent in filecontents:
        # Parsing the hr file contents
        hr, rvect, weights, nrpts = utils.parse_wannier_hr_file_contents(filecontent)

        # Sanity checking
        if weights_out is None:
            weights_out = weights
        elif weights != weights_out:
            raise ValueError('Cannot merge HR file contents that have differing weights.')
        if rvect_out is None:
            rvect_out = rvect
        elif np.all(rvect != rvect_out):
            raise ValueError('Cannot merge HR file contents with differing sets of R-vectors.')

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

    rvect_out_list = rvect_out.tolist()
    assert isinstance(rvect_out_list, list) and all([isinstance(sublist, list) for sublist in rvect_out_list]) and all(
        [isinstance(x, int) for sublist in rvect_out_list for x in sublist])

    return utils.generate_wannier_hr_file_contents(hr_out, rvect_out_list, weights_out)


def merge_wannier_u_file_contents(filecontents: List[str]) -> str:
    """Merge the contents of multiple `Wannier90` U matrix files into a single file."""
    u_list = []
    kpts_master = None
    for filecontent in filecontents:
        # Parsing the U file contents
        umat, kpts, nkpts = utils.parse_wannier_u_file_contents(filecontent)

        if kpts_master is None:
            kpts_master = kpts
        elif nkpts == len(kpts_master) and np.allclose(kpts, kpts_master):
            pass
        else:
            raise ValueError('Cannot merge U matrix file contents with differing sets of k-points.')

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
    return utils.generate_wannier_u_file_contents(u_merged, kpts)


def merge_wannier_centers_file_contents(filecontents: list[str], atoms: Atoms) -> str:
    """Merge the contents of multiple `Wannier90` centers files into a single file."""
    centers_list = []
    for filecontent in filecontents:
        centers, _ = utils.parse_wannier_centers_file_contents(filecontent)

        centers_list += centers

    # Writing the centers file
    return utils.generate_wannier_centers_file_contents(centers_list, atoms)


def extend_wannier_u_dis_file_content(filecontent: str, nbnd: int, nwann: int) -> str:
    """Extend the content of a `Wannier90` u_dis file to include the identity matrix for lower bands."""
    # Parse the file content
    udis_mat, kpts, _ = utils.parse_wannier_u_file_contents(filecontent)

    # Build up the larger U_dis matrix, which is a nkpts x nwann_emp x nbnd_emp matrix...
    udis_mat_large = np.zeros((len(kpts), nwann, nbnd), dtype=complex)
    # ... with the diagonal entries equal to 1...
    udis_mat_large[:, :nwann, :nwann] = np.identity(nwann)
    # ... except for the last block, where we insert the contents of the corresponding u_dis file
    udis_mat_large[:, -udis_mat.shape[1]:, -udis_mat.shape[2]:] = udis_mat

    # Generate the contents of the extended U_dis file
    return utils.generate_wannier_u_file_contents(udis_mat_large, kpts)


class MergeInputModel(IOModel):
    """Input model for a `MergeProcess`."""

    src_files: List[File]
    dst_file: Path
    model_config = ConfigDict(arbitrary_types_allowed=True)


class MergeOutputModel(IOModel):
    """Output model for a `MergeProcess`."""

    dst_file: File
    model_config = ConfigDict(arbitrary_types_allowed=True)


class MergeProcess(Process[MergeInputModel, MergeOutputModel]):
    """Generic process for merging the contents of multiple Wannier files (e.g. hr, u, centers)."""

    input_model = MergeInputModel
    output_model = MergeOutputModel

    def __init__(self, merge_function: Callable[[List[str]], str], **kwargs):
        self.merge_function = merge_function
        super().__init__(**kwargs)

    def _run(self):
        if len(self.inputs.src_files) == 0:
            raise ValueError('No input files provided to merge.')

        filecontents = [f.read_text() for f in self.inputs.src_files]

        merged_filecontents = self.merge_function(filecontents)

        dst_file = File(self, self.inputs.dst_file)
        dst_file.write_text(merged_filecontents)

        self.outputs = self.output_model(dst_file=dst_file)


class ExtendInputModel(IOModel):
    """Input model for an `ExtendProcess`."""

    src_file: File
    dst_file: Path
    model_config = ConfigDict(arbitrary_types_allowed=True)


class ExtendOutputModel(IOModel):
    """Output model for an `ExtendProcess`."""

    dst_file: File
    model_config = ConfigDict(arbitrary_types_allowed=True)


class ExtendProcess(Process):
    """Generic process for extending the contents of a Wannier file (usually a u_dis file)."""

    input_model = ExtendInputModel
    output_model = ExtendOutputModel

    def __init__(self, extend_function: Callable[[str], str], **kwargs):
        self.extend_function = extend_function
        super().__init__(**kwargs)

    def _run(self):

        filecontent = self.inputs.src_file.read_text()

        extended_filecontent = self.extend_function(filecontent)

        dst_file = File(self, self.inputs.dst_file)

        dst_file.write_text(extended_filecontent)

        self.outputs = self.output_model(dst_file=dst_file)
