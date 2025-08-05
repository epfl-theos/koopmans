"""Processes for calling WannierJL."""

from functools import lru_cache

from pydantic import field_validator

from koopmans.files import File
from koopmans.process_io import IOModel

from ._process import Process

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
    return juliapkg, jl


@lru_cache
def load_wannierjl():
    """Load the WannierJL julia module."""
    juliapkg, jl = import_julia()
    juliapkg.add("Wannier", uuid=WANNIER_JL_UUID, rev=WANNIER_JL_REV)
    jl.seval("using Wannier")
    return jl


class WannierJLCheckNeighborsInput(IOModel):
    """Input model for the WannierJLCheckNeighborsProcess."""

    wannier90_input_file: File
    chk_file: File
    mmn_file: File
    amn_file: File
    eig_file: File


class WannierJLCheckNeighborsOutput(IOModel):
    """Output model for the WannierJLCheckNeighborsProcess."""

    has_cubic_neighbors: bool


class WannierJLCheckNeighborsProcess(Process[WannierJLCheckNeighborsInput, WannierJLCheckNeighborsOutput]):
    """Process for checking if the default set of b-vectors contains six cubic nearest neighbors."""

    input_model = WannierJLCheckNeighborsInput
    output_model = WannierJLCheckNeighborsOutput

    def _run(self):
        # Load the Wannier julia module
        jl = load_wannierjl()

        # Link the files locally
        File(self, 'wannier90.win').symlink_to(self.inputs.wannier90_input_file)
        File(self, 'wannier90.chk').symlink_to(self.inputs.chk_file)
        File(self, 'wannier90.mmn').symlink_to(self.inputs.mmn_file)
        File(self, 'wannier90.amn').symlink_to(self.inputs.amn_file)
        File(self, 'wannier90.eig').symlink_to(self.inputs.eig_file)

        print([f for f in self.directory.rglob('*')])  # Debugging line to see the files in the directory

        # Read the Wannier files
        model = jl.read_w90_with_chk(str(self.inputs.wannier90_input_file.with_suffix("")), str(self.inputs.chk_file))

        # Check if the default set of b-vectors contains six cubic nearest neighbors
        has_cubic_neighbors = jl.Wannier.has_cubic_neighbors(model.kstencil)

        # Save the output
        self.outputs = self.output_model(has_cubic_neighbors=has_cubic_neighbors)


class WannierJLGenerateNeighborsInput(IOModel):
    """Input model for the WannierJLGenerateNeighborsProcess."""

    wannier90_input_file: File


class WannierJLGenerateNeighborsOutput(IOModel):
    """Output model for the WannierJLGenerateNeighborsProcess."""

    nnkp_file: File


class WannierJLGenerateNeighborsProcess(Process[WannierJLGenerateNeighborsInput, WannierJLGenerateNeighborsOutput]):
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


class WannierJLSplitInput(IOModel):
    """Input model for the WannierJLSplitProcess."""

    indices: list[list[int]]
    outdirs: list[str]
    wannier90_input_file: File
    chk_file: File
    mmn_file: File
    amn_file: File
    eig_file: File
    cubic_nnkp_file: File | None = None
    cubic_mmn_file: File | None = None

    @field_validator("wannier90_input_file", "chk_file", "mmn_file", "amn_file", "eig_file", "cubic_nnkp_file",
                     "cubic_mmn_file", mode="after")
    @classmethod
    def exists(cls, value: File | None) -> File | None:
        """Ensure the Wannier90 input file exists."""
        if value is not None and not value.exists():
            raise FileNotFoundError(f"{value} does not exist.")
        return value


class WannierJLSplitBlockOutput(IOModel):
    """Output for one block within the WannierJLSplitProcess."""

    amn_file: File
    eig_file: File
    mmn_file: File
    u_file: File
    win_file: File

    @field_validator("amn_file", "eig_file", "mmn_file", "u_file", mode="after")
    @classmethod
    def exists(cls, value: File | None) -> File | None:
        """Ensure the Wannier90 input file exists."""
        if value is not None and not value.exists():
            raise FileNotFoundError(f"{value} does not exist.")
        return value


class WannierJLSplitOutput(IOModel):
    """Output model for the WannierJLSplitProcess."""

    blocks: list[WannierJLSplitBlockOutput]


class WannierJLSplitProcess(Process[WannierJLSplitInput, WannierJLSplitOutput]):
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
        jl.Wannier.Tools.mrwf(str(self.inputs.wannier90_input_file.with_suffix("")),
                              julia_indices,
                              julia_outdirs,
                              cubic_mmn_file)

        # Save the output
        prefix = str(self.inputs.wannier90_input_file.with_suffix("").name)
        self.outputs = WannierJLSplitOutput(
            blocks=[{'amn_file': File(self, d) / (prefix + ".amn"),
                     'eig_file': File(self, d) / (prefix + ".eig"),
                     'mmn_file': File(self, d) / (prefix + ".mmn"),
                     'u_file': File(self, d) / (prefix + "_split.amn"),
                     'win_file': File(self, d) / (prefix + ".win")}
                    for d in self.inputs.outdirs]
        )
