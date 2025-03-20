import shutil

import pytest
from pytest_console_scripts import script_runner

from koopmans import __version__
from koopmans.utils import chdir


def test_cli_version(script_runner):
    result = script_runner.run(["koopmans", "--version"])
    assert result.success
    assert result.stdout.strip() == __version__


def test_cli_help(script_runner):
    result = script_runner.run(["koopmans", "--help"])
    assert result.success
    assert result.stdout.startswith("usage: koopmans")


def test_cli_no_args(script_runner):
    result = script_runner.run(["koopmans"])
    assert result.success
    assert result.stdout.startswith("usage: koopmans")


def test_cli_run_no_args(script_runner):
    result = script_runner.run(["koopmans", "run"])
    assert not result.success


def test_cli_run_si(script_runner, workflow_patch, tmpdir, datadir):
    with chdir(tmpdir):
        shutil.copy(datadir / "json" / "si.json", tmpdir)
        result = script_runner.run(["koopmans", "run", "si.json"])
    assert result.success


def test_cli_run_si_with_logs(script_runner, workflow_patch, tmpdir, datadir):
    with chdir(tmpdir):
        shutil.copy(datadir / "json" / "si.json", tmpdir)
        result = script_runner.run(["koopmans", "run", "si.json", "--log"])
    assert result.success
    assert (tmpdir / "koopmans.log").exists()


class TestPseudos:
    def test_cli_pseudos_list(self, script_runner):
        result = script_runner.run(["koopmans", "pseudos", "list"])
        assert result.success
        assert "PseudoDojo/0.4/LDA/SR/standard/upf" in result.stdout
        assert result.stderr == ''

    @pytest.mark.dependency()
    def test_cli_pseudos_install(self, script_runner, datadir):
        install_result = script_runner.run(
            ["koopmans", "pseudos", "install", f"{datadir}/pseudos/H_ONCV_PBE-1.2.upf", "--library", "TestPseudos"])
        assert install_result.success

    @pytest.mark.dependency(depends=["TestPseudos::test_cli_pseudos_install[inprocess]"])
    def test_cli_pseudos_list_after_install(self, script_runner):
        list_result = script_runner.run(["koopmans", "pseudos", "list"])
        assert list_result.success
        assert "TestPseudos" in list_result.stdout

    @pytest.mark.dependency(depends=["TestPseudos::test_cli_pseudos_list_after_install[inprocess]"])
    def test_cli_pseudos_uninstall(self, script_runner):
        uninstall_result = script_runner.run(["koopmans", "pseudos", "uninstall", "TestPseudos"])
        assert uninstall_result.success
        assert "TestPseudos" not in uninstall_result.stdout

    @pytest.mark.dependency(depends=["TestPseudos::test_cli_pseudos_uninstall[inprocess]"])
    def test_cli_pseudos_list_after_uninstall(self, script_runner):
        second_list_result = script_runner.run(["koopmans", "pseudos", "list"])
        assert second_list_result.success
        assert "TestPseudos" not in second_list_result.stdout
