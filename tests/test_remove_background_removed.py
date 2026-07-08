"""
Validates that the CellBender / mitobender "remove-background" functionality has
been removed from mgatk and only survives as a graceful, no-op deprecation stub.

The heavy background-removal implementation used to live in-tree; it now points
users at the external mitobender project. These tests guard against the
implementation (or its imports) sneaking back in.
"""
import importlib.util

import pytest
from click.testing import CliRunner

from mgatk import cli


@pytest.fixture
def runner():
    return CliRunner()


def test_remove_background_is_a_deprecation_stub(runner, tmp_path):
    result = runner.invoke(cli.main, ["remove-background", "-i", "."])
    assert result.exit_code == 0, result.output
    assert "deprecated" in result.output
    assert "mitobender" in result.output


def test_remove_background_does_no_processing(runner, tmp_path):
    out_dir = tmp_path / "out"
    result = runner.invoke(
        cli.main, ["remove-background", "-i", ".", "-o", str(out_dir)]
    )
    assert result.exit_code == 0, result.output
    # The stub must not spin up any pipeline output directories.
    assert not out_dir.exists()


@pytest.mark.parametrize(
    "module_name",
    [
        "mgatk.cellbender",
        "mgatk.mitobender",
        "mgatk.processing.cellbender",
        "mgatk.processing.remove_background",
        "mgatk.remove_background",
    ],
)
def test_no_background_removal_implementation_modules(module_name):
    assert importlib.util.find_spec(module_name) is None


def test_deprecated_cellbender_flags_are_hidden_noops(runner):
    """The old --ncells_fg/--ncells_bg CellBender knobs remain accepted (so old
    scripts don't break) but are hidden and do nothing."""
    help_result = runner.invoke(cli.main, ["--help"])
    assert help_result.exit_code == 0
    assert "ncells_fg" not in help_result.output
    assert "ncells_bg" not in help_result.output
    assert "CellBender" not in help_result.output

    # Still parsed without error when explicitly supplied.
    result = runner.invoke(
        cli.main, ["remove-background", "-i", ".", "-nfg", "5", "-nbg", "9"]
    )
    assert result.exit_code == 0, result.output
