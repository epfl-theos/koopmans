"""Processes for calling WannierJL."""

import juliapkg
from juliacall import Main as jl
from pydantic import ConfigDict

from koopmans.files import File
from koopmans.process_io import IOModel

from ._process import Process


class WannierJLSplitVCInput(IOModel):
    pass


class WannierJLSplitVCOutput(IOModel):
    pass


class WannierJLSplitVCProcess(Process[WannierJLSplitVCInput, WannierJLSplitVCOutput]):

    def _run(self):
        pass


class WannierJLCheckNNKPInput(IOModel):
    pass


class WannierJLCheckNNKPOutput(IOModel):

    missing_bvectors: bool


class WannierJLCheckNNKPProcess(Process[WannierJLCheckNNKPInput, WannierJLCheckNNKPOutput]):

    input_model = WannierJLCheckNNKPInput
    output_model = WannierJLCheckNNKPOutput

    def _run(self):
        self.outputs = self.output_model(missing_bvectors=True)


class WannierJLGenerateCubicNNKPInput(IOModel):
    wannier_input_file: File

    model_config = ConfigDict(arbitrary_types_allowed=True)


class WannierJLGenerateCubicNNKPOutput(IOModel):
    nnkp_file: File

    model_config = ConfigDict(arbitrary_types_allowed=True)


class WannierJLGenerateCubicNNKPProcess(Process[WannierJLGenerateCubicNNKPInput, WannierJLGenerateCubicNNKPOutput]):

    input_model = WannierJLGenerateCubicNNKPInput
    output_model = WannierJLGenerateCubicNNKPOutput

    def _run(self):
        # Load the Wannier julia module
        juliapkg.add("Wannier", uuid="2b19380a-1f7e-4d7d-b1b8-8aa60b3321c9", rev="bvec_cubic")
        juliapkg.resolve()
        jl.seval("using Wannier")

        # Parse the Wannier input file
        win = jl.read_win(str(self.inputs.wannier_input_file))

        # Construct the reciprocal lattice
        reciprocal_lattice = jl.Wannier.get_recip_lattice(win.unit_cell_cart)

        # Generate the k-point stencil
        kstencil = jl.Wannier.get_bvectors_nearest(win.kpoints, win.mp_grid, reciprocal_lattice)

        # Write the nnkp file
        nnkp_file = File(self, 'cubic.nnkp')
        jl.write_nnkp(str(nnkp_file), kstencil, win.num_wann)

        self.outputs = self.output_model(nnkp_file=nnkp_file)
