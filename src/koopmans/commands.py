"""Pydantic models that define the configurations of commands."""

import os
from abc import ABC, abstractmethod
from typing import Any, ClassVar

from pydantic import Field, field_validator, model_validator
from typing_extensions import Self

from koopmans.base import BaseModel


def _redirection_string(stdin: str, stdout: str | None, stderr: None | str = None,
                        stdin_redirection_string: str = "<") -> str:
    redirection = f"{stdin_redirection_string} {stdin}"
    if stdout is not None:
        redirection += f" > {stdout}"
        if stderr is not None:
            redirection += f" 2> {stderr}"
        else:
            redirection += " 2>&1"
    return redirection


def para_prefix_to_dict(para_prefix: str) -> dict[str, str | None | int]:
    """Extract the MPI executable and number of tasks from the PARA_PREFIX string."""
    parts = para_prefix.split(" ")
    if len(parts) == 3:
        [mpi_executable, _, num_tasks_str] = parts
        num_tasks = int(num_tasks_str)
    elif len(parts) == 1:
        mpi_executable = parts[0]
        num_tasks = None
    else:
        raise ValueError(
            f"Invalid PARA_PREFIX string `{para_prefix}`\nExpected something of the format "
            "`<mpi_executable> -n <num_tasks>` or `<mpi_executable>`")

    return {
        "mpi_executable": mpi_executable,
        "num_tasks": num_tasks,
    }


def command_string_to_dict(command_string: str) -> dict[str, Any]:
    """Convert a command string into a dictionary of its constituent parts."""
    # Split the string into parts
    parts = command_string.split(" ")[::-1]
    stdout = None
    stderr = None
    options: dict[str, Any] = {}
    try:
        executable = parts.pop()

        # Command-line options
        while parts[-1].startswith("-"):
            if parts[-1] == "-in":
                break
            key = parts.pop().strip('-')
            value: Any
            if parts[-1].startswith("-") or parts[-1] == "<" or len(parts) < 2 or parts[-2] == ">":
                # This must be a boolean flag
                value = True
            else:
                value = parts.pop().strip('.')
            options[key] = value
        if parts[-1] in ["-in", "<"]:
            stdin_redirection_string = parts.pop()
        else:
            stdin_redirection_string = ""
        stdin = parts.pop()
        if parts and parts.pop() == ">":
            stdout = parts.pop()
            if parts and parts[-1] == "2>&1":
                stderr = None
            elif parts.pop() == "2>":
                stderr = parts.pop()
    except IndexError:
        raise ValueError(
            f"Invalid string `{command_string}`\nExpected something of the format `<executable> "
            "(-<flag> <value>) <stdin_redirection> <stdin> > <stdout> 2>&1` or `2> <stderr>`")

    return {
        "executable": executable,
        "stdin": stdin,
        "stdout": stdout,
        "stderr": stderr,
        "stdin_redirection_string": stdin_redirection_string} | options


class CommandConfig(BaseModel, ABC):
    """Abstract base class for a command configuration."""

    executable: str = Field(..., description="Path to the command executable.")
    stdin: str = Field(..., description="Input file for the command.")
    stdout: str | None = Field(None, description="Output file for the command.")
    stderr: str | None = Field(None, description="Error file for the command.")
    stdin_redirection_string: ClassVar[str] = "<"
    ase_env_var: ClassVar[str] = "ASE_COMMAND"

    @property
    @abstractmethod
    def options_str(self) -> str:
        """Return the string of options to provide to the executable."""

    @property
    def redirection_string(self) -> str:
        """Return the string of redirection to provide to the executable."""
        return _redirection_string(
            self.stdin, self.stdout, self.stderr, self.stdin_redirection_string
        )

    def command(self, additional_flags: list[str] = [], additional_options: dict[str, str] = {}) -> str:
        """Return the string of the command to run."""
        additional_options_str = ' '.join([f"{k} {v}" for k, v in additional_options.items()])

        command_str = (f"{self.executable} {self.options_str} {' '.join(additional_flags)} "
                       f"{additional_options_str} {self.redirection_string}")

        while '  ' in command_str:
            command_str = command_str.replace('  ', ' ')

        return command_str.strip()

    @classmethod
    def from_string(cls, value: str, **kwargs: Any) -> Self:
        """Create a CommandConfig object from a string."""
        kwargs.update(command_string_to_dict(value))
        stdin_redirect = kwargs.pop("stdin_redirection_string")
        if stdin_redirect != cls.stdin_redirection_string:
            raise ValueError(
                f"Redirection string `{stdin_redirect}` does not match expected `{cls.stdin_redirection_string}`")

        return cls(**kwargs)

    @classmethod
    def from_env(cls) -> Self:
        """Create a CommandConfig based on the corresponding environment variable."""
        env_var = os.environ.get(cls.ase_env_var)
        if env_var is None:
            raise OSError(f"Environment variable `{cls.ase_env_var}` not set")
        return cls.from_string(env_var)


class ParallelCommandConfig(CommandConfig, ABC):
    """Abstract base class for a command configuration that can be run in parallel."""

    mpi_executable: str = Field("mpirun", description="MPI executable to use for parallel execution.")
    num_tasks: int | None = Field(None, description="Number of tasks to run in parallel.")

    def command(self, additional_flags: list[str] = [], additional_options: dict[str, str] = {}) -> str:
        """Return the string of the command to run."""
        if self.num_tasks == 1:
            prefix = ""
        else:
            ntasks_str = f" -n {self.num_tasks}" if self.num_tasks is not None else ""
            prefix = f"{self.mpi_executable}{ntasks_str} "
        return prefix + super().command(additional_flags, additional_options)

    @classmethod
    def from_string(cls, value: str, **kwargs: Any) -> Self:
        """Create a ParallelCommandConfig object from a string."""
        if value.split(" ")[0] not in ["srun", "mpirun"]:
            mpi_executable = "mpirun"
            num_tasks = 1
            rest = value
        elif value.split(" ")[1] not in ["-n", "--np"]:
            mpi_executable, rest = value.split(" ", 1)
            num_tasks = None
        else:
            mpi_executable, _, num_tasks_str, rest = value.split(" ", 3)
            num_tasks = int(num_tasks_str)

        kwargs["mpi_executable"] = mpi_executable
        kwargs["num_tasks"] = num_tasks

        return super(ParallelCommandConfig, cls).from_string(rest, **kwargs)


class KCWHamConfig(ParallelCommandConfig):
    """Configuration for the kcw.x command."""

    executable: str = "kcw.x"
    stdin: str = "PREFIX.khi"
    stdout: str | None = "PREFIX.kho"
    stderr: str | None = None
    stdin_redirection_string: ClassVar[str] = "-in"
    ase_env_var: ClassVar[str] = "ASE_KCW_HAM_COMMAND"

    @property
    def options_str(self) -> str:
        """Return the string of options to provide to the executable."""
        return ''


class KCWScreenConfig(ParallelCommandConfig):
    """Configuration for the kcw.x command."""

    executable: str = "kcw.x"
    stdin: str = "PREFIX.ksi"
    stdout: str | None = "PREFIX.kso"
    stderr: str | None = None
    stdin_redirection_string: ClassVar[str] = "-in"
    ase_env_var: ClassVar[str] = "ASE_KCW_SCREEN_COMMAND"

    @property
    def options_str(self) -> str:
        """Return the string of options to provide to the executable."""
        return ''


class KCWWannierConfig(ParallelCommandConfig):
    """Configuration for the kcw.x command."""

    executable: str = "kcw.x"
    stdin: str = "PREFIX.w2ki"
    stdout: str | None = "PREFIX.w2ko"
    stderr: str | None = None
    stdin_redirection_string: ClassVar[str] = "-in"
    ase_env_var: ClassVar[str] = "ASE_KCW_WANNIER_COMMAND"

    @property
    def options_str(self) -> str:
        """Return the string of options to provide to the executable."""
        return ''


class KCPConfig(ParallelCommandConfig):
    """Configuration for the kcp.x command."""

    executable: str = "kcp.x"
    stdin: str = "PREFIX.cpi"
    stdout: str | None = "PREFIX.cpo"
    stderr: str | None = None
    stdin_redirection_string: ClassVar[str] = "-in"
    ase_env_var: ClassVar[str] = "ASE_KCP_COMMAND"

    @property
    def options_str(self) -> str:
        """Return the string of options to provide to the executable."""
        return ''


class PhConfig(ParallelCommandConfig):
    """Configuration for the ph.x command."""

    executable: str = "ph.x"
    stdin: str = "PREFIX.phi"
    stdout: str | None = "PREFIX.pho"
    stderr: str | None = None
    stdin_redirection_string: ClassVar[str] = "-in"
    ase_env_var: ClassVar[str] = "ASE_PH_COMMAND"

    @property
    def options_str(self) -> str:
        """Return the string of options to provide to the executable."""
        return ''


class ProjwfcConfig(ParallelCommandConfig):
    """Configuration for the projwfc.x command."""

    executable: str = "projwfc.x"
    stdin: str = "PREFIX.pri"
    stdout: str | None = "PREFIX.pro"
    stderr: str | None = None
    stdin_redirection_string: ClassVar[str] = "-in"
    ase_env_var: ClassVar[str] = "ASE_PROJWFC_COMMAND"

    @property
    def options_str(self) -> str:
        """Return the string of options to provide to the executable."""
        return ''


class Wann2KCPConfig(ParallelCommandConfig):
    """Configuration for the wann2kcp.x command."""

    executable: str = "wann2kcp.x"
    stdin: str = "PREFIX.wki"
    stdout: str | None = "PREFIX.wko"
    stderr: str | None = None
    stdin_redirection_string: ClassVar[str] = "-in"
    ase_env_var: ClassVar[str] = "ASE_WANN2KCP_COMMAND"

    @property
    def options_str(self) -> str:
        """Return the string of options to provide to the executable."""
        return ''


class PWConfig(ParallelCommandConfig):
    """Configuration for the pw.x command."""

    executable: str = "pw.x"
    stdin: str = "PREFIX.pwi"
    stdout: str | None = "PREFIX.pwo"
    stderr: str | None = None
    stdin_redirection_string: ClassVar[str] = "-in"
    npool: int | None = Field(None, description="Number of pools to use for parallel execution.")
    pd: bool = Field(False, description="Use pencil decomposition.")
    ase_env_var: ClassVar[str] = "ASE_PW_COMMAND"

    @property
    def options_str(self) -> str:
        """Return the string of options to provide to the executable."""
        options: list[str] = []
        if self.npool is not None:
            options.append(f"-npool {self.npool}")
        if self.pd:
            options.append(f"-pd {self.pd}")
        return ' '.join(options)


class PW2Wannier90Config(ParallelCommandConfig):
    """Configuration for the pw2wannier90.x command."""

    executable: str = "pw2wannier90.x"
    stdin: str = "PREFIX.p2wi"
    stdout: str | None = "PREFIX.p2wo"
    stderr: str | None = None
    stdin_redirection_string: ClassVar[str] = "-in"
    pd: bool = Field(False, description="Use pencil decomposition.")
    ase_env_var: ClassVar[str] = "ASE_PW2WANNIER90_COMMAND"

    @property
    def options_str(self) -> str:
        """Return the string of options to provide to the executable."""
        options: list[str] = []
        if self.pd:
            options.append(f"-pd {self.pd}")
        return ' '.join(options)


class Wannier90Config(CommandConfig):
    """Configuration for the wannier90.x command."""

    executable: str = "wannier90.x"
    stdin: str = "PREFIX.win"
    stdin_redirection_string: ClassVar[str] = ""
    pp: bool = Field(False, description="Perform preprocessing.")
    ase_env_var: ClassVar[str] = "ASE_WANNIER90_COMMAND"

    @field_validator("stdout", mode="after")
    @classmethod
    def check_stdout_is_none(cls, value: str | None) -> None:
        """Check that stdout is None."""
        if value is not None:
            raise ValueError("stdout must be None for wannier90.x as it automatically performs redirection.")
        return value

    @property
    def options_str(self) -> str:
        """Return the string of options to provide to the executable."""
        return '-pp' if self.pp else ''

    @property
    def redirection_string(self) -> str:
        """Return the string of redirection to provide to the executable.

        Here we override the default redirection string because Wannier90 automatically
        performs redirection.
        """
        return f"{self.stdin}"


class WannierJLConfig(CommandConfig):
    """Configuration for the wannierjl command."""

    executable: str = "wjl"
    stdin: str = "PREFIX.wjli"
    stdout: str | None = "PREFIX.wjlo"
    stderr: str | None = None
    stdin_redirection_string: ClassVar[str] = ">"
    ase_env_var: ClassVar[str] = "ASE_WANNIERJL_COMMAND"

    @property
    def options_str(self) -> str:
        """Return the string of options to provide to the executable."""
        return ''


class CommandConfigs(BaseModel):
    """Collection of command configurations.

    Will consult the environment variables `ASE_<CALC_NAME>_COMMAND` to set the default values for each command.
    Note that the various "type: ignore"s exist because the type checker isn't aware that we populate missing
    fields during the model validators
    """

    mpi_executable: str = Field(
        default="mpirun",
        description="MPI executable to use for parallel execution."
    )
    kcw_ham: KCWHamConfig = Field(
        default_factory=lambda: KCWHamConfig(),  # type: ignore
        description="Configuration for the kcw.x code when running `calculation='ham'`; "
        "defaults to `ASE_KCW_HAMILTONIAN_COMMAND` if set."
    )
    kcw_screen: KCWScreenConfig = Field(
        default_factory=lambda: KCWScreenConfig(),  # type: ignore
        description="Configuration for the kcw.x code when running `calculation='screen'`; "
        "defaults to `ASE_KCW_SCREENING_COMMAND` if set."
    )
    kcw_wannier: KCWWannierConfig = Field(
        default_factory=lambda: KCWWannierConfig(),   # type: ignore
        description="Configuration for the kcw.x code when running `calculation='wann2kcw'`; "
        "defaults to `ASE_KCW_WANNIER_COMMAND` if set."
    )
    kcp: KCPConfig = Field(
        default_factory=lambda: KCPConfig(),  # type: ignore
        description="Configuration for the kcp.x code; defaults to `ASE_KCP_COMMAND` if set."
    )
    ph: PhConfig = Field(
        default_factory=lambda: PhConfig(),  # type: ignore
        description="Configuration for the ph.x code; defaults to `ASE_PH_COMMAND` if set."
    )
    projwfc: ProjwfcConfig = Field(
        default_factory=lambda: ProjwfcConfig(),  # type: ignore
        description="Configuration for the projwfc.x code; defaults to `ASE_PROJWFC_COMMAND` if set."
    )
    pw2wannier90: PW2Wannier90Config = Field(
        default_factory=lambda: PW2Wannier90Config(),  # type: ignore
        description="Configuration for the pw2wannier90.x code; defaults to `ASE_PW2WANNIER90_COMMAND` if set."
    )
    pw: PWConfig = Field(
        default_factory=lambda: PWConfig(),  # type: ignore
        description="Configuration for the pw.x code; defaults to `ASE_PW_COMMAND` if set."
    )
    wann2kcp: Wann2KCPConfig = Field(
        default_factory=lambda: Wann2KCPConfig(),  # type: ignore
        description="Configuration for the wann2kcp.x code; defaults to `ASE_WANN2KCP_COMMAND` if set."
    )
    wannier90: Wannier90Config = Field(
        default_factory=lambda: Wannier90Config(),  # type: ignore
        description="Configuration for the wannier90.x code; defaults to `ASE_WANNIER90_COMMAND` if set."
    )
    wannierjl: WannierJLConfig = Field(
        default_factory=lambda: WannierJLConfig(),  # type: ignore
        description="Configuration for the WannierJL code; defaults to `ASE_WANNIERJL_COMMAND` if set."
    )

    @model_validator(mode="before")
    def default_mpi_executable(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Use the provided `mpi_executable` for all codes."""
        # The default value
        mpi_executable = cls.__pydantic_fields__['mpi_executable'].default

        # PARA_PREFIX takes precedence over the default value
        para_prefix = os.environ.get("PARA_PREFIX", None)
        if para_prefix is not None:
            mpi_settings = para_prefix_to_dict(para_prefix)
            mpi_executable = mpi_settings["mpi_executable"]
            num_tasks = mpi_settings["num_tasks"]
        else:
            num_tasks = None

        # Value provided explicitly takes precedence over both
        provided_mpi_executable = data.get("mpi_executable", None)
        if provided_mpi_executable is not None:
            mpi_executable = provided_mpi_executable

        data["mpi_executable"] = mpi_executable

        for name, field in cls.__pydantic_fields__.items():
            # Check if the field corresponds to a ParallelCommandConfig
            subclass = field.annotation
            if subclass is None or not issubclass(subclass, ParallelCommandConfig):
                continue

            # Propagate mpi_executable and num_tasks to the configuration
            code_settings = data.get(name, {})
            if isinstance(code_settings, dict):
                # Only propagate to configurations if they were not provided explicitly or were
                # provided as a dictionary (in which case the mpi_executable propagation is convenient)
                data[name] = {"mpi_executable": mpi_executable, "num_tasks": num_tasks} | code_settings

        return data

    @model_validator(mode="before")
    def get_env_vars(cls, data: dict[str, Any]) -> dict[str, Any]:
        """Set the configurations based on the corresponding environment variables, if present."""
        for name, field in cls.__pydantic_fields__.items():
            # Check if the field corresponds to a CommandConfig
            subclass = field.annotation
            if subclass is None or not issubclass(subclass, CommandConfig):
                continue

            # If the field is not provided, use the environment variable
            env_var = os.environ.get(subclass.ase_env_var)
            if env_var is not None and name not in data:
                data[name] = subclass.from_string(env_var)

        return data
