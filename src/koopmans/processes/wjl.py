"""Processes for calling WannierJL."""

from functools import lru_cache

from pydantic import field_validator

from koopmans.files import File, FileViaSymlink
from koopmans.pydantic import BaseModel

from .pythonprocess import PythonProcess

WANNIER_JL_UUID = "2b19380a-1f7e-4d7d-b1b8-8aa60b3321c9"
WANNIER_JL_REV = "65245c59"


@lru_cache
def import_julia():
    """Attempt to import the julia package and the Main module from juliacall."""
    try:
        import juliapkg
        from juliacall import Main as jl
    except ImportError:
        raise ImportError("The installed version of `koopmans` does not contain support for `WannierJL`. "
                          "Run `pip install koopmans[julia]` to enable this functionality.")

    # Capturing stdout
    jl.seval("""
    function capture_stdio(f, output_file, args...; kwargs...)
        result = nothing
        redirect_stdio(stdout=output_file, stderr=output_file) do
            result = f(args...; kwargs...)
        end
        return result
    end
    """)
    return juliapkg, jl


@lru_cache
def load_wannierjl():
    """Load the WannierJL julia module."""
    juliapkg, jl = import_julia()
    juliapkg.add("Wannier", uuid=WANNIER_JL_UUID, rev=WANNIER_JL_REV)
    jl.seval("using Wannier")
    return jl


class WannierJLCheckNeighborsInput(BaseModel):
    """Input model for the WannierJLCheckNeighborsProcess."""

    wannier90_input_file: FileViaSymlink
    chk_file: FileViaSymlink
    mmn_file: FileViaSymlink
    amn_file: FileViaSymlink
    eig_file: FileViaSymlink


class WannierJLCheckNeighborsOutput(BaseModel):
    """Output model for the WannierJLCheckNeighborsProcess."""

    has_cubic_neighbors: bool
    julia_stdio: File


class WannierJLCheckNeighborsProcess(PythonProcess[WannierJLCheckNeighborsInput, WannierJLCheckNeighborsOutput]):
    """Process for checking if the default set of b-vectors contains six cubic nearest neighbors."""

    input_model = WannierJLCheckNeighborsInput
    output_model = WannierJLCheckNeighborsOutput

    def _run(self):
        # Load the Wannier julia module
        jl = load_wannierjl()

        # Read the Wannier files
        model = jl.capture_stdio(jl.read_w90_with_chk,
                                 str(self.directory / "julia.out"),
                                 str(self.directory / self.inputs.wannier90_input_file.dest.with_suffix("")))

        # Check if the default set of b-vectors contains six cubic nearest neighbors
        has_cubic_neighbors = jl.Wannier.has_cubic_neighbors(model.kstencil)

        # Save the output
        self.outputs = self.output_model(has_cubic_neighbors=has_cubic_neighbors,
                                         julia_stdio=File(parent_process_uid=self.uid, path="julia.out"))


class WannierJLGenerateNeighborsInput(BaseModel):
    """Input model for the WannierJLGenerateNeighborsProcess."""

    wannier90_input_file: File


class WannierJLGenerateNeighborsOutput(BaseModel):
    """Output model for the WannierJLGenerateNeighborsProcess."""

    nnkp_file: File


class WannierJLGenerateNeighborsProcess(PythonProcess[WannierJLGenerateNeighborsInput, WannierJLGenerateNeighborsOutput]):
    """Process for generating a nnkp file with six cubic nearest neighbors."""

    input_model = WannierJLGenerateNeighborsInput
    output_model = WannierJLGenerateNeighborsOutput

    def _run(self):
        # Load the Wannier julia module
        jl = load_wannierjl()

        # Define the output nnkp file
        nnkp_file = File(self, 'cubic.nnkp')

        # Generate the nnkp file
        jl.Wannier.write_nnkp_cubic(str(nnkp_file), str(self.inputs.wannier90_input_file))

        # Save the output
        self.outputs = self.output_model(nnkp_file=nnkp_file)


class WannierJLSplitInput(BaseModel):
    """Input model for the WannierJLSplitProcess."""

    indices: list[list[int]]
    outdirs: list[str]
    wannier90_input_file: FileViaSymlink
    chk_file: FileViaSymlink
    mmn_file: FileViaSymlink
    amn_file: FileViaSymlink
    eig_file: FileViaSymlink
    cubic_mmn_file: FileViaSymlink | None = None


class WannierJLSplitBlockOutput(BaseModel):
    """Output for one block within the WannierJLSplitProcess."""

    amn_file: File
    eig_file: File
    mmn_file: File
    u_file: File
    win_file: File


class WannierJLSplitOutput(BaseModel):
    """Output model for the WannierJLSplitProcess."""

    blocks: list[WannierJLSplitBlockOutput]
    julia_stdio: File


class WannierJLSplitProcess(PythonProcess[WannierJLSplitInput, WannierJLSplitOutput]):
    """Process for splitting a Wannier manifold into multiple blocks using WannierJL."""

    input_model = WannierJLSplitInput
    output_model = WannierJLSplitOutput

    def _run(self):
        # Load the Wannier julia module
        jl = load_wannierjl()

        # Construct the julia indices (need to do this so that the following mrwf interface finds a match)
        julia_indices = jl.seval("[" + ', '.join([str(i) for i in self.inputs.indices]) + "]")
        julia_outdirs = jl.seval("[" + ', '.join([f'"{self.directory / d}"' for d in self.inputs.outdirs]) + "]")

        # Perform the parallel transport algorithm to split the manifolds
        cubic_mmn_file = str(self.inputs.cubic_mmn_file) if self.inputs.cubic_mmn_file else None
        jl.capture_stdio(jl.Wannier.Tools.mrwf,
                         str(self.directory / "julia.out"),
                         str(self.directory / self.inputs.wannier90_input_file.dest.with_suffix("")),
                         julia_indices,
                         julia_outdirs,
                         cubic_mmn_file)

        # Save the output
        prefix = str(self.inputs.wannier90_input_file.dest.with_suffix(""))
        self.outputs = WannierJLSplitOutput(
            blocks=[{'amn_file': File(parent_process_uid=self.uid, path=d) / (prefix + ".amn"),
                     'eig_file': File(parent_process_uid=self.uid, path=d) / (prefix + ".eig"),
                     'mmn_file': File(parent_process_uid=self.uid, path=d) / (prefix + ".mmn"),
                     'u_file': File(parent_process_uid=self.uid, path=d) / (prefix + "_split.amn"),
                     'win_file': File(parent_process_uid=self.uid, path=d) / (prefix + ".win")}
                    for d in self.inputs.outdirs],
            julia_stdio=File(parent_process_uid=self.uid, path="julia.out"),
        )
