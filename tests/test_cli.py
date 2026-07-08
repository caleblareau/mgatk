import gzip

import pytest
from click.testing import CliRunner

from mgatk import cli
from conftest import requires_mgatk_r_packages, requires_snakemake


@pytest.fixture
def runner():
    return CliRunner()


# ---------------------------------------------------------------------------
# Argument validation (fast, no external tools required)
# ---------------------------------------------------------------------------

def test_check_without_barcode_tag_fails(runner, tmp_path):
    result = runner.invoke(
        cli.main,
        ["check", "-i", "input", "-o", str(tmp_path / "out"), "-n", "name"],
    )
    assert result.exit_code != 0
    assert "must specify a valid read tag ID" in str(result.exception)


def test_check_rejects_non_bam_input(runner, tmp_path):
    not_a_bam = tmp_path / "notabam.txt"
    not_a_bam.write_text("hello")
    result = runner.invoke(
        cli.main,
        ["check", "-i", str(not_a_bam), "-o", str(tmp_path / "out"), "-n", "name", "-bt", "CB"],
    )
    assert result.exit_code != 0
    assert "should be an individual .bam file" in str(result.exception)


def test_check_missing_bam_file_fails(runner, tmp_path):
    missing_bam = tmp_path / "missing.bam"
    result = runner.invoke(
        cli.main,
        ["check", "-i", str(missing_bam), "-o", str(tmp_path / "out"), "-n", "name", "-bt", "CB"],
    )
    assert result.exit_code != 0
    assert "No file found" in str(result.exception)


def test_support_lists_builtin_genomes(runner):
    result = runner.invoke(cli.main, ["support", "-i", "."])
    assert result.exit_code != 0  # `support` mode always exits after printing
    assert "rCRS" in result.output
    assert "hg19" in result.output


def test_check_passes_with_valid_bam(runner, tmp_path, humanbam_dir):
    # `check` mode always exits via sys.exit(message) once validation succeeds,
    # which Python reports as exit code 1 even though it's the "happy path".
    bam = humanbam_dir / "MGH60-P6-A11.mito.bam"
    result = runner.invoke(
        cli.main,
        [
            "check", "-i", str(bam), "-o", str(tmp_path / "out"), "-n", "name",
            "-bt", "CB", "-g", "hg19",
        ],
    )
    assert isinstance(result.exception, SystemExit)
    assert "mgatk check passed" in str(result.exception)


# ---------------------------------------------------------------------------
# End-to-end smoke tests against the bundled sample data.
# These require snakemake and the R/Bioconductor packages mgatk uses to build
# final output; skipped automatically if absent.
# ---------------------------------------------------------------------------

@requires_snakemake
@requires_mgatk_r_packages
def test_call_end_to_end_on_humanbam(runner, tmp_path, humanbam_dir):
    out_dir = tmp_path / "out"
    result = runner.invoke(
        cli.main,
        [
            "call", "-i", str(humanbam_dir), "-o", str(out_dir), "-n", "glio",
            "-g", "hg19", "--snake-stdout",
        ],
    )
    assert result.exit_code == 0, result.output

    final = out_dir / "final"
    assert (final / "glio.rds").exists()
    assert (final / "glio.depthTable.txt").exists()
    for base in "ACGT":
        assert (final / f"glio.{base}.txt.gz").exists()

    with gzip.open(final / "glio.coverage.txt.gz", "rt") as fh:
        header_line = fh.readline()
    assert header_line.strip() != ""


@requires_snakemake
@requires_mgatk_r_packages
def test_bcall_end_to_end_on_known_barcodes(runner, tmp_path, barcode_dir):
    out_dir = tmp_path / "out"
    result = runner.invoke(
        cli.main,
        [
            "bcall", "-i", str(barcode_dir / "test_barcode.bam"), "-n", "bc1",
            "-o", str(out_dir), "-bt", "CB", "-b", str(barcode_dir / "test_barcodes.txt"),
            "-z", "--snake-stdout",
        ],
    )
    assert result.exit_code == 0, result.output
    assert (out_dir / "final" / "bc1.rds").exists()


@requires_snakemake
@requires_mgatk_r_packages
def test_tenx_end_to_end_on_known_barcodes(runner, tmp_path, barcode_dir):
    out_dir = tmp_path / "out"
    result = runner.invoke(
        cli.main,
        [
            "tenx", "-i", str(barcode_dir / "test_barcode.bam"), "-n", "bc1",
            "-o", str(out_dir), "-bt", "CB", "-b", str(barcode_dir / "test_barcodes.txt"),
            "-c", "2", "--snake-stdout",
        ],
    )
    assert result.exit_code == 0, result.output

    final = out_dir / "final"
    assert (final / "bc1.rds").exists()
    with gzip.open(final / "bc1.variant_stats.tsv.gz", "rt") as fh:
        header = fh.readline().strip().split("\t")
    assert header == [
        "position", "nucleotide", "variant", "vmr", "mean", "variance",
        "n_cells_conf_detected", "n_cells_over_5", "n_cells_over_10", "n_cells_over_20",
        "n_cells_over_95", "max_heteroplasmy", "strand_correlation", "mean_coverage",
    ]
