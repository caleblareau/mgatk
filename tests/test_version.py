"""
Regression tests for CLI version handling.

PR #107 moved version resolution off of the removed ``pkg_resources`` API.
A later merge left ``@click.version_option(version=__version__, ...)`` at module
import time while the module-level ``__version__`` binding was dropped, so every
console script crashed with ``NameError`` on import. These tests lock in that the
version is resolvable at import time and surfaced through ``--version``.
"""
from importlib.metadata import version

import mgatk
from mgatk import cli
from mgatk.deletioncalling import clidel, clifind

from click.testing import CliRunner


EXPECTED_VERSION = version("mgatk")


def test_package_version_matches_installed_metadata():
    assert mgatk.__version__ == EXPECTED_VERSION


def test_cli_version_option():
    result = CliRunner().invoke(cli.main, ["--version"])
    assert result.exit_code == 0, result.output
    assert "mgatk" in result.output
    assert EXPECTED_VERSION in result.output


def test_del_version_option():
    result = CliRunner().invoke(clidel.main, ["--version"])
    assert result.exit_code == 0, result.output
    assert "mgatk-del" in result.output
    assert EXPECTED_VERSION in result.output


def test_del_find_version_option():
    result = CliRunner().invoke(clifind.main, ["--version"])
    assert result.exit_code == 0, result.output
    assert "mgatk-del-find" in result.output
    assert EXPECTED_VERSION in result.output
