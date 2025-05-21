"""Tests for the CommandConfig subclasses."""

from koopmans.commands import (CommandConfigs, KCPConfig, KCWHamConfig,
                               PWConfig, Wannier90Config)


def test_kcw_config():
    """Test the KCWHamConfig class."""
    config = KCWHamConfig(mpi_executable="srun", num_tasks=8, stdin="PREFIX.khi", stdout="PREFIX.kho")

    assert config.executable == "kcw.x"
    assert config.stdin == "PREFIX.khi"
    assert config.stdout == "PREFIX.kho"

    assert config.command() == "srun -n 8 kcw.x -in PREFIX.khi > PREFIX.kho 2>&1"

    second_config = KCWHamConfig.from_string(config.command())
    assert second_config == config


def test_kcp_config():
    """Test the KCPConfig class."""
    config = KCPConfig(mpi_executable="mpirun", num_tasks=4)

    assert config.executable == "kcp.x"
    assert config.stdin == "PREFIX.cpi"
    assert config.stdout == "PREFIX.cpo"
    assert config.command() == "mpirun -n 4 kcp.x -in PREFIX.cpi > PREFIX.cpo 2>&1"

    second_config = KCPConfig.from_string(config.command())
    assert second_config == config


def test_calculator_configs_ase_command(monkeypatch):
    """Test the CalculatorConfigs class when ASE_KCW_HAM_COMMAND is set."""
    custom_command = "mpirun -n 3 kcw.x -in PREFIX.khi > PREFIX.kho 2> PREFIX.err"
    monkeypatch.delenv("PARA_PREFIX", False)
    monkeypatch.setenv("ASE_KCW_HAM_COMMAND", custom_command)

    config = CommandConfigs()
    assert config.kcw_ham.command() == custom_command


def test_calculator_configs_para_prefix(monkeypatch):
    """Test the CalculatorConfigs class when PARA_PREFIX is set."""
    # Make sure we don't have any environment variables
    monkeypatch.delenv("ASE_KCW_HAM_COMMAND", False)
    monkeypatch.delenv("ASE_PW_COMMAND", False)
    monkeypatch.delenv("PARA_PREFIX", False)

    # Check by default that mpirun is used
    config = CommandConfigs()
    assert config.kcw_ham.mpi_executable == "mpirun"

    # Check that when PARA_PREFIX is set it propagates to all subclasses
    monkeypatch.setenv("PARA_PREFIX", "srun")
    config = CommandConfigs()
    assert config.kcw_ham.mpi_executable == "srun"
    assert config.kcw_ham.num_tasks is None

    # Check that explicitly providing `mpi_executable` to a calculator takes precedence over PARA_PREFIX
    config = CommandConfigs(pw={"mpi_executable": "mpirun"})
    assert config.kcw_ham.mpi_executable == "srun"
    assert config.pw.mpi_executable == "mpirun"

    # Check that explicitly providing `mpi_executable` to the entire config class takes precedence over PARA_PREFIX
    config = CommandConfigs(mpi_executable="srun")
    assert config.kcw_ham.mpi_executable == "srun"
    assert config.pw.mpi_executable == "srun"


def test_pw_config_from_string():
    """Test PWConfig.from_string()."""
    command = "mpirun -n 4 pw.x -npool 4 -pd .true. -in PREFIX.pwi > PREFIX.pwo 2>&1"
    config = PWConfig.from_string(command)
    assert config.pd
    assert config.npool == 4


def test_wannier90_preproc_from_string():
    """Test Wannier90PreprocConfig.from_string()."""
    command = "wannier90.x -pp PREFIX.w90"
    config = Wannier90Config.from_string(command)
    assert config.command() == command
    assert config.executable == "wannier90.x"
    assert config.stdin == "PREFIX.w90"
    assert config.pp
    assert config.command() == command
